#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple>
#include <map>

struct FiveDiagonalMatrix {
    std::vector<double> mainDiagonal;      // Главная диагональ
    std::vector<double> upperDiagonal1;     // Первая верхняя диагональ
    std::vector<double> upperDiagonal2;     // Вторая верхняя диагональ
    std::vector<double> lowerDiagonal1;     // Первая нижняя диагональ
    std::vector<double> lowerDiagonal2;     // Вторая нижняя диагональ
};
// Параметры задачи
double lamda = 1.0;  // Коэффициент теплопроводности
double gamma = 1.0;  // Коэффициент источника тепла

// Функция для вычисления коэффициента вне диагонали
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
    int type;     // Тип граничного условия
};

// Структура для хранения узла сетки
struct GridNode {
    double x, y;
    std::string type; // Тип узла: "фиктивный", "граничный", "внутренний"
    int boundaryType; // Тип граничного условия (если узел граничный)
    int index;        // Индекс узла в матрице
};

// Проверка, лежит ли точка на отрезке
bool isPointOnLineSegment(const Point& p1, const Point& p2, const Point& p) {
    double cross = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
    if (std::abs(cross) > 1e-12) return false; // Точка не на прямой

    // Проверка, лежит ли точка внутри отрезка
    if (p1.x != p2.x) { // Горизонтальная линия
        if ((p.x < std::min(p1.x, p2.x)) || (p.x > std::max(p1.x, p2.x))) return false;
    }
    else { // Вертикальная линия
        if ((p.y < std::min(p1.y, p2.y)) || (p.y > std::max(p1.y, p2.y))) return false;
    }
    return true;
}

// Проверка, находится ли точка внутри многоугольника
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

// Функция для вычисления значения функции на границе
double boundaryFunction(double x, double y) {
    return x + y; // Пример функции
}

// Чтение входных данных из файла
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

// Создание сетки
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

// Построение 5-диагональной матрицы и правой части
void buildFiveDiagonalMatrix(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    FiveDiagonalMatrix& matrix, std::vector<double>& b) {
    int N = grid.size() * grid[0].size();
    matrix.mainDiagonal.resize(N, 0.0);
    matrix.upperDiagonal1.resize(N - 1, 0.0);
    matrix.upperDiagonal2.resize(N - 2, 0.0);
    matrix.lowerDiagonal1.resize(N - 1, 0.0);
    matrix.lowerDiagonal2.resize(N - 2, 0.0);
    b.assign(N, 0.0);

    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            int idx = grid[i][j].index;
            if (grid[i][j].type == "фиктивный") {
                matrix.mainDiagonal[idx] = 1.0;
                b[idx] = 0.0;
                continue;
            }

            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                matrix.mainDiagonal[idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = x_nodes[j] - x_nodes[j - 1];
                double hy = y_nodes[i] - y_nodes[i - 1];
                matrix.mainDiagonal[idx] = _cu(hx, hy);

                if (j < grid[i].size() - 1 && grid[i][j + 1].type != "фиктивный") {
                    matrix.upperDiagonal1[idx] = _yu(hx);
                }
                if (j < grid[i].size() - 2 && grid[i][j + 2].type != "фиктивный") {
                    matrix.upperDiagonal2[idx] = _yu(hx);
                }
                if (j > 0 && grid[i][j - 1].type != "фиктивный") {
                    matrix.lowerDiagonal1[idx] = _yu(hy);
                }
                if (j > 1 && grid[i][j - 2].type != "фиктивный") {
                    matrix.lowerDiagonal2[idx] = _yu(hy);
                }
            }
        }
    }
}

// Решение СЛАУ методом блочной релаксации
void blockRelaxation(const FiveDiagonalMatrix& matrix, const std::vector<double>& b, std::vector<double>& x, int maxIterations = 10000, double tolerance = 1e-12) {
    int N = b.size();
    x.assign(N, 0.0); // Начальное приближение (нулевой вектор)

    for (int iter = 0; iter < maxIterations; ++iter) {
        double maxError = 0.0;

        for (int i = 0; i < N; ++i) {
            double sum = 0.0;

            // Учитываем элементы из нижних диагоналей
            if (i >= 1) {
                sum += matrix.lowerDiagonal1[i - 1] * x[i - 1];
            }
            if (i >= 2) {
                sum += matrix.lowerDiagonal2[i - 2] * x[i - 2];
            }

            // Учитываем элементы из верхних диагоналей
            if (i < N - 1) {
                sum += matrix.upperDiagonal1[i] * x[i + 1];
            }
            if (i < N - 2) {
                sum += matrix.upperDiagonal2[i] * x[i + 2];
            }

            // Вычисляем новое значение x[i]
            double newX = (b[i] - sum) / matrix.mainDiagonal[i];

            // Обновляем максимальную ошибку
            maxError = std::max(maxError, std::abs(newX - x[i]));

            // Обновляем значение x[i]
            x[i] = newX;
        }

        // Проверка на сходимость
        if (maxError < tolerance) {
            std::cout << "Метод сошелся за " << iter + 1 << " итераций." << std::endl;
            return;
        }
    }

    std::cout << "Метод не сошелся за " << maxIterations << " итераций." << std::endl;
}

// Запись результатов в файл
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

    // Чтение данных
    double x_min, x_max, y_min, y_max;
    int x_steps, y_steps;
    std::vector<BoundaryCondition> boundaryConditions;
    readInputData("input.txt", x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions);
    std::vector<Point> polygon;
    for (const auto& condition : boundaryConditions) {
        polygon.push_back(condition.p1);
        polygon.push_back(condition.p2);
    }

    // Создание сетки
    std::vector<std::vector<GridNode>> grid;
    std::vector<double> x_nodes, y_nodes;
    std::tie(grid, x_nodes, y_nodes) = createGrid(x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions, polygon);

    // Построение 5-диагональной матрицы и правой части
    FiveDiagonalMatrix matrix;
    std::vector<double> b;
    buildFiveDiagonalMatrix(grid, x_nodes, y_nodes, matrix, b);

    // Решение СЛАУ методом блочной релаксации
    std::vector<double> x;
    blockRelaxation(matrix, b, x);

    // Запись результатов в файл
    writeResultsToFile("results.txt", grid, x, b);

    return 0;
}