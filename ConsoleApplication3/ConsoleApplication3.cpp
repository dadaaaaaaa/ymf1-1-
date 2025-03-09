#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple> 

double lamda = 1.0;
double gamma = 1.0;


double _yu(double h) {
    return (-lamda / (h * h));
}

// Функция для вычисления центрального коэффициента
double _cu(double h, double hy) {
    return (2 * lamda * (1 / (h * h) + 1 / (hy * hy)) + gamma);
}

// Структура для хранения точки
struct Point {
    double x, y;
};

// Структура для хранения граничного условия
struct BoundaryCondition {
    Point p1, p2; // Точки, образующие прямую
    int type; // Тип граничного условия
};

// Структура для хранения узла сетки
struct GridNode {
    double x, y;
    std::string type; // Тип узла: "фиктивный", "граничный", "внутренний"
    int boundaryType; // Тип граничного условия (если узел граничный)
    int index; // Индекс узла в матрице
};

bool isPointOnLineSegment(const Point& p1, const Point& p2, const Point& p) {
    // Проверка, лежит ли точка на прямой
    double cross = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
    if (std::abs(cross) > 1e-12) return false; // Точка не на прямой

    // Проверка, лежит ли точка внутри отрезка
    if (p1.x != p2.x) { // Горизонтальная линия
        if ((p.x < std::min(p1.x, p2.x)) || (p.x > std::max(p1.x, p2.x))) return false;
    }
    else { // Вертикальная линия
        if ((p.y < std::min(p1.y, p2.y)) || (p.y > std::max(p1.y, p2.y))) return false;
    }
    // Проверка, лежит ли точка внутри отрезка
    double dot = (p.x - p1.x) * (p2.x - p1.x) + (p.y - p1.y) * (p2.y - p1.y);
    if (dot < 0) return false; // Точка за пределами отрезка

    double squaredLength = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
    if (dot > squaredLength) return false; // Точка за пределами отрезка
    return true;
}

// Функция для проверки, находится ли точка внутри многоугольника
bool isPointInsidePolygon(const std::vector<Point>& polygon, const Point& p) {
    int n = polygon.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (((polygon[i].y > p.y) != (polygon[j].y > p.y)) &&
            (p.x < (polygon[j].x - polygon[i].x) * (p.y - polygon[i].y) / (polygon[j].y - polygon[i].y) + polygon[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

// Функция для вычисления значения функции на границе (пример)
double boundaryFunction(double x, double y) {
    return x + y; // Пример функции
}
void readInputData(const std::string& filename, double& x_min, double& x_max, int& x_steps,
    double& y_min, double& y_max, int& y_steps,
    std::vector<BoundaryCondition>& boundaryConditions) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        exit(1);
    }

    file >> x_min >> x_max >> x_steps;
    file >> y_min >> y_max >> y_steps;

    BoundaryCondition bc;
    while (file >> bc.p1.x >> bc.p1.y >> bc.p2.x >> bc.p2.y >> bc.type) {
        boundaryConditions.push_back(bc);
    }
    file.close();
}


std::tuple<std::vector<std::vector<GridNode>>, std::vector<double>, std::vector<double>>
createGrid(double x_min, double x_max, int x_steps,
    double y_min, double y_max, int y_steps,
    const std::vector<BoundaryCondition>& boundaryConditions,
    const std::vector<Point>& polygon) {
    std::vector<double> x_nodes(x_steps + 1);
    std::vector<double> y_nodes(y_steps + 1);
    for (int i = 0; i <= x_steps; ++i) {
        x_nodes[i] = x_min + i * (x_max - x_min) / x_steps;
    }
    for (int i = 0; i <= y_steps; ++i) {
        y_nodes[i] = y_min + i * (y_max - y_min) / y_steps;
    }

    std::vector<std::vector<GridNode>> grid(y_steps + 1, std::vector<GridNode>(x_steps + 1));
    int nodeIndex = 0;
    for (int i = 0; i <= y_steps; ++i) {
        for (int j = 0; j <= x_steps; ++j) {
            grid[i][j].x = x_nodes[j];
            grid[i][j].y = y_nodes[i];
            grid[i][j].boundaryType = -1;
            grid[i][j].index = nodeIndex++;

            bool isBoundary = false;
            for (const auto& condition : boundaryConditions) {
                if (isPointOnLineSegment(condition.p1, condition.p2, { grid[i][j].x, grid[i][j].y })) {
                    isBoundary = true;
                    grid[i][j].boundaryType = condition.type;
                    break;
                }
            }

            if (isBoundary) {
                grid[i][j].type = "граничный";
            }
            else if (isPointInsidePolygon(polygon, { grid[i][j].x, grid[i][j].y })) {
                grid[i][j].type = "внутренний";
            }
            else {
                grid[i][j].type = "фиктивный";
            }
        }
    }

    return std::make_tuple(grid, x_nodes, y_nodes);
}

void buildMatrixAndRightHandSide(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int N = grid.size() * grid[0].size();
    A.assign(N, std::vector<double>(N, 0.0));
    b.assign(N, 0.0);

    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            if (grid[i][j].type == "фиктивный") continue;

            int idx = grid[i][j].index;

            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                A[idx][idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = x_nodes[j] - x_nodes[j - 1];
                double hy = y_nodes[i] - y_nodes[i - 1];

                A[idx][idx] = _cu(hx, hy);

                if (j > 0 && grid[i][j - 1].type != "фиктивный") {
                    A[idx][grid[i][j - 1].index] = _yu(hx);
                }
                if (j < grid[i].size() - 1 && grid[i][j + 1].type != "фиктивный") {
                    A[idx][grid[i][j + 1].index] = _yu(hx);
                }
                if (i > 0 && grid[i - 1][j].type != "фиктивный") {
                    A[idx][grid[i - 1][j].index] = _yu(hy);
                }
                if (i < grid.size() - 1 && grid[i + 1][j].type != "фиктивный") {
                    A[idx][grid[i + 1][j].index] = _yu(hy);
                }
            }
        }
    }
}

void printResults(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    std::cout << "Номера узлов и их координаты:" << std::endl;
    for (const auto& row : grid) {
        for (const auto& node : row) {
            std::cout << "Узел " << node.index << ": (" << node.x << ", " << node.y << "), тип: " << node.type;
            if (node.type == "граничный") {
                std::cout << ", тип граничного условия: " << node.boundaryType;
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Матрица A:" << std::endl;
    for (const auto& row : A) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Правая часть b:" << std::endl;
    for (double val : b) {
        std::cout << val << std::endl;
    }
}


int main() {
    std::setlocale(LC_ALL, "Russian");
    double x_min, x_max, y_min, y_max;
    int x_steps, y_steps;
    std::vector<BoundaryCondition> boundaryConditions;
    readInputData("input.txt", x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions);
    std::vector<Point> polygon;
    for (const auto& condition : boundaryConditions) {
        polygon.push_back(condition.p1);
        polygon.push_back(condition.p2);
    }

    std::vector<std::vector<GridNode>> grid;
    std::vector<double> x_nodes, y_nodes;
    std::tie(grid, x_nodes, y_nodes) = createGrid(x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions, polygon);
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    buildMatrixAndRightHandSide(grid, x_nodes, y_nodes, A, b);

    // Вывод результатов
    printResults(grid, A, b);

    return 0;
}