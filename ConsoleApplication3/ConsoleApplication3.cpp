#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple>
#include <map>

double lamda = 1.0;
double gamma = 1.0;

double _yu(double h) {
    return (-lamda / (h * h));
}

double _cu(double h, double hy) {
    return (2 * lamda * (1 / (h * h) + 1 / (hy * hy)) + gamma);
}

double boundaryFunction(double x, double y) {
    return x * x * x + y * y * y; // Пример функции
}

struct Point {
    double x, y;
};

struct BoundaryCondition {
    Point p1, p2;
    int type;
};

struct GridNode {
    double x, y;
    std::string type;
    int boundaryType;
    int index;
};

bool isPointOnLineSegment(const Point& p1, const Point& p2, const Point& p) {
    double cross = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
    if (std::abs(cross) > 1e-12) return false;

    if (p1.x != p2.x) {
        if ((p.x < std::min(p1.x, p2.x)) || (p.x > std::max(p1.x, p2.x))) return false;
    }
    else {
        if ((p.y < std::min(p1.y, p2.y)) || (p.y > std::max(p1.y, p2.y))) return false;
    }

    double dot = (p.x - p1.x) * (p2.x - p1.x) + (p.y - p1.y) * (p2.y - p1.y);
    if (dot < 0) return false;

    double squaredLength = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
    if (dot > squaredLength) return false;
    return true;
}

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

void printResults(const std::vector<std::vector<GridNode>>& grid, const std::vector<double>& b) {
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
    std::cout << "Правая часть b:" << std::endl;
    for (int i = 0; i < b.size(); i++) {
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
    }
}

std::vector<double> solveWithGaussSeidel(const std::vector<double>& di, const std::vector<double>& u1, const std::vector<double>& u2,
    const std::vector<double>& l1, const std::vector<double>& l2, const std::vector<double>& b,
    int maxIterations = 1000, double tolerance = 1e-12) {
    int N = di.size();
    std::vector<double> x(N, 0.0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        std::vector<double> x_old = x;
        double maxError = 0.0;

        for (int i = 0; i < N; ++i) {
            double sum = b[i];

            if (i > 0) sum -= l2[i - 1] * x[i - 1];
            if (i > 1) sum -= l1[i - 2] * x[i - 2];
            if (i < N - 1) sum -= u2[i] * x[i + 1];
            if (i < N - 2) sum -= u1[i] * x[i + 2];

            x[i] = sum / di[i];
            maxError = std::max(maxError, std::abs(x[i] - x_old[i]));
        }

        if (maxError < tolerance) {
            std::cout << "Метод Гаусса-Зейделя сошелся за " << iter + 1 << " итераций." << std::endl;
            return x;
        }
    }

    std::cout << "Метод Гаусса-Зейделя не сошелся за " << maxIterations << " итераций." << std::endl;
    return x;
}

void buildMatrixAndRightHandSides(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    std::vector<double>& di, std::vector<double>& u1, std::vector<double>& u2,
    std::vector<double>& l1, std::vector<double>& l2, std::vector<double>& b, int x_steps, int y_steps) {
    int N = grid.size() * grid[0].size();
    di.assign(N, 0.0);
    u1.assign(N, 0.0);
    u2.assign(N, 0.0);
    l1.assign(N, 0.0);
    l2.assign(N, 0.0);
    b.assign(N, 0.0);

    int rows = grid.size();
    int cols = grid[0].size();
    // Заполнение матрицы
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int idx = i * cols + j; // Индекс текущего узла

            if (grid[i][j].type == "фиктивный") {
                di[idx] = 1.0;
                b[idx] = 0;
                continue;
            }

            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                di[idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = (j > 0) ? (x_nodes[j] - x_nodes[j - 1]) : 1.0; // Шаг по X
                double hy = (i > 0) ? (y_nodes[i] - y_nodes[i - 1]) : 1.0; // Шаг по Y

                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
                di[idx] = _cu(hx, hy);

                // Заполнение нижних диагоналей
                if (j > 0 && grid[i][j - 1].type != "фиктивный") {
                    l1[idx] = _yu(hx); // Первая нижняя диагональ
                }
                if (i > 0 && grid[i - 1][j].type != "фиктивный") {
                    l2[idx] = _yu(hy); // Вторая нижняя диагональ
                }

                // Заполнение верхних диагоналей
                if (j < cols - 1 && grid[i][j + 1].type != "фиктивный") {
                    u1[idx] = _yu(hx); // Первая верхняя диагональ
                }
                if (i < rows - 1 && grid[i + 1][j].type != "фиктивный") {
                    u2[idx] = _yu(hy); // Вторая верхняя диагональ
                }
            }
        }
    }

    // Для проверки корректности заполнения
    for (int i = 0; i < N; ++i) {
        std::cout << "di[" << i << "] = " << di[i] << ", ";
        if (i < N - 1) std::cout << "u1[" << i << "] = " << u1[i] << ", ";
        if (i < N - cols - 1) std::cout << "u2[" << i << "] = " << u2[i] << ", ";
        if (i > 0) std::cout << "l1[" << i << "] = " << l1[i] << ", ";
        if (i > cols - 1) std::cout << "l2[" << i << "] = " << l2[i] << ", ";
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
    }

}

void writeResultsToFile(const std::string& filename, const std::vector<std::vector<GridNode>>& grid, const std::vector<double>& solutions, const std::vector<double>& b) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Ошибка открытия файла для записи!" << std::endl;
        return;
    }

    outFile << "N \tx\ty\tu\tu*\t|u - u*|\n";
    for (const auto& row : grid) {
        for (const auto& node : row) {
            int idx = node.index;
            double difference = std::abs(solutions[idx] - b[idx]);
            outFile << node.index << "\t" << node.x << "\t" << node.y << "\t"
                << b[idx] << "\t" << solutions[idx] << "\t" << difference << "\n";
        }
    }

    outFile.close();
    std::cout << "Результаты успешно сохранены в файл: " << filename << std::endl;
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

    std::vector<double> di, u1, u2, l1, l2, b;
    buildMatrixAndRightHandSides(grid, x_nodes, y_nodes, di, u1, u2, l1, l2, b, x_steps, y_steps);

    std::vector<double> x = solveWithGaussSeidel(di, u1, u2, l1, l2, b);

    writeResultsToFile("results.txt", grid, x, b);

    printResults(grid, b);

    std::cout << "Решение СЛАУ:" << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    std::cout << "Решение Ошибка:" << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        std::cout << "x[" << i << "] = " << x[i] - b[i] << std::endl;
    }

    return 0;
}