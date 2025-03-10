#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple> 
#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple> 
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple> 
#include <map>
// Структура для хранения 5-ти диагональной матрицы
struct FiveDiagonalMatrix {
    std::vector<double> mainDiagonal;      // Главная диагональ
    std::vector<double> upperDiagonal1;     // Первая верхняя диагональ
    std::vector<double> upperDiagonal2;     // Вторая верхняя диагональ
    std::vector<double> lowerDiagonal1;     // Первая нижняя диагональ
    std::vector<double> lowerDiagonal2;     // Вторая нижняя диагональ
};

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
    const FiveDiagonalMatrix& matrix, const std::vector<double>& b) {
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

    std::cout << "5-ти диагональная матрица:" << std::endl;

    // Вывод главной диагонали
    std::cout << "Главная диагональ: ";
    for (double val : matrix.mainDiagonal) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Вывод первой верхней диагонали
    std::cout << "Первая верхняя диагональ: ";
    for (double val : matrix.upperDiagonal1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Вывод второй верхней диагонали
    std::cout << "Вторая верхняя диагональ: ";
    for (double val : matrix.upperDiagonal2) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Вывод первой нижней диагонали
    std::cout << "Первая нижняя диагональ: ";
    for (double val : matrix.lowerDiagonal1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Вывод второй нижней диагонали
    std::cout << "Вторая нижняя диагональ: ";
    for (double val : matrix.lowerDiagonal2) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Правая часть b:" << std::endl;
    for (int i = 0; i < b.size(); i++) {
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
    }
}


// Функция для сохранения 5-ти диагональной матрицы в файл
void saveFiveDiagonalMatrixToFile(const FiveDiagonalMatrix& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла для записи!" << std::endl;
        return;
    }

    // Сохраняем главную диагональ
    for (double val : matrix.mainDiagonal) {
        file << val << " ";
    }
    file << std::endl;

    // Сохраняем первую верхнюю диагональ
    for (double val : matrix.upperDiagonal1) {
        file << val << " ";
    }
    file << std::endl;

    // Сохраняем вторую верхнюю диагональ
    for (double val : matrix.upperDiagonal2) {
        file << val << " ";
    }
    file << std::endl;

    // Сохраняем первую нижнюю диагональ
    for (double val : matrix.lowerDiagonal1) {
        file << val << " ";
    }
    file << std::endl;

    // Сохраняем вторую нижнюю диагональ
    for (double val : matrix.lowerDiagonal2) {
        file << val << " ";
    }
    file << std::endl;

    file.close();
    std::cout << "5-ти диагональная матрица успешно сохранена в файл: " << filename << std::endl;
}

// Функция для построения 5-ти диагональной матрицы и правой части
void buildFiveDiagonalMatrix(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    FiveDiagonalMatrix& matrix, std::vector<double>& b) {
    int N = grid.size() * grid[0].size();
    matrix.mainDiagonal.resize(N, 0.0);
    matrix.upperDiagonal1.resize(N - 1, 0.0); // Для первой верхней диагонали
    matrix.upperDiagonal2.resize(N - 2, 0.0); // Для второй верхней диагонали
    matrix.lowerDiagonal1.resize(N - 1, 0.0); // Для первой нижней диагонали
    matrix.lowerDiagonal2.resize(N - 2, 0.0); // Для второй нижней диагонали
    b.assign(N, 0.0);

    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            int idx = grid[i][j].index;
            if (grid[i][j].type == "фиктивный") { matrix.mainDiagonal[idx] = 1.0; b[idx] = 0.0; continue; }



            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                matrix.mainDiagonal[idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = x_nodes[j] - x_nodes[j - 1];
                double hy = y_nodes[i] - y_nodes[i - 1];
                // Заполнение главной диагонали
                matrix.mainDiagonal[idx] = _cu(hx, hy);

                // Заполнение верхних диагоналей
                if (j < grid[i].size() - 1 && grid[i][j + 1].type != "фиктивный") {
                    matrix.upperDiagonal1[idx] = _yu(hx); // Первая верхняя диагональ
                }
                if (j < grid[i].size() - 2 && grid[i][j + 2].type != "фиктивный") {
                    matrix.upperDiagonal2[idx] = _yu(hx); // Вторая верхняя диагональ
                }

                // Заполнение нижних диагоналей
                if (j > 0 && grid[i][j - 1].type != "фиктивный") {
                    matrix.lowerDiagonal1[idx] = _yu(hy); // Первая нижняя диагональ
                }
                if (j > 1 && grid[i][j - 2].type != "фиктивный") {
                    matrix.lowerDiagonal2[idx] = _yu(hy); // Вторая нижняя диагональ
                }
            }
        }
    }
}

// Функция для решения СЛАУ методом блочной релаксации
void blockRelaxation(const FiveDiagonalMatrix& matrix, const std::vector<double>& b, std::vector<double>& x, int maxIterations = 1000, double tolerance = 1e-6) {
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
// Функция для создания полной матрицы из 5-ти диагональной
std::vector<std::vector<double>> createFullMatrix(const FiveDiagonalMatrix& matrix, int N) {
    std::vector<std::vector<double>> fullMatrix(N, std::vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i) {
        fullMatrix[i][i] = matrix.mainDiagonal[i];
        if (i < N - 1) {
            fullMatrix[i][i + 1] = matrix.upperDiagonal1[i];
            fullMatrix[i + 1][i] = matrix.lowerDiagonal1[i]; // Симметрично
        }
        if (i < N - 2) {
            fullMatrix[i][i + 2] = matrix.upperDiagonal2[i];
            fullMatrix[i + 2][i] = matrix.lowerDiagonal2[i]; // Симметрично
        }
    }

    return fullMatrix;
}

// Функция для сохранения полной матрицы в файл
void saveFullMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла для записи полной матрицы!" << std::endl;
        return;
    }

    for (const auto& row : matrix) {
        for (double val : row) {
            file << val << " ";
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Полная матрица успешно сохранена в файл: " << filename << std::endl;
}

// Функция для решения СЛАУ методом Гаусса
std::vector<double> solveWithGaussianElimination(const std::vector<std::vector<double>>& matrix, const std::vector<double>& b) {
    int N = matrix.size();
    std::vector<double> x(N, 0.0);
    std::vector<std::vector<double>> augmentedMatrix(N, std::vector<double>(N + 1));

    // Создаем расширенную матрицу
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            augmentedMatrix[i][j] = matrix[i][j];
        }
        augmentedMatrix[i][N] = b[i];
    }
    saveFullMatrixToFile(augmentedMatrix, "full_matrix_output.txt");
    // Прямой ход Гаусса
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double factor = augmentedMatrix[j][i] / augmentedMatrix[i][i];
            for (int k = i; k <= N; ++k) {
                augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
            }
        }
    }

    // Обратный ход
    for (int i = N - 1; i >= 0; --i) {
        x[i] = augmentedMatrix[i][N] / augmentedMatrix[i][i];
        for (int j = i - 1; j >= 0; --j) {
            augmentedMatrix[j][N] -= augmentedMatrix[j][i] * x[i];
        }
    }

    return x;
}
void buildMatrixAndRightHandSides(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int N = grid.size() * grid[0].size();
    A.assign(N, std::vector<double>(N, 0.0));
    b.assign(N, 0.0);

    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            int idx = grid[i][j].index;
            if (grid[i][j].type == "фиктивный") { A[idx][idx] = 1.0;  b[idx] = 0; continue; }



            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                A[idx][idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = x_nodes[j] - x_nodes[j - 1];
                double hy = y_nodes[i] - y_nodes[i - 1];
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
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

int main() {
    std::setlocale(LC_ALL, "Russian");

    // Чтение данных
    double x_min, x_max, y_min;
    double y_max;
    int x_steps, y_steps;
    std::vector<BoundaryCondition> boundaryConditions;
    readInputData("input.txt", x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions);

    // Создание полигона
    std::vector<Point> polygon;
    for (const auto& condition : boundaryConditions) {
        polygon.push_back(condition.p1);
        polygon.push_back(condition.p2);
    }

    // Создание сетки и получение x_nodes и y_nodes
    std::vector<std::vector<GridNode>> grid;
    std::vector<double> x_nodes, y_nodes;
    std::tie(grid, x_nodes, y_nodes) = createGrid(x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions, polygon);

    // Построение 5-ти диагональной матрицы и правой части
    FiveDiagonalMatrix matrix;
    std::vector<double> b;
    buildFiveDiagonalMatrix(grid, x_nodes, y_nodes, matrix, b);

    // Построение матрицы и правой части
    std::vector<std::vector<double>> A;
    buildMatrixAndRightHandSides(grid, x_nodes, y_nodes, A, b);

    // Решение СЛАУ методом Гаусса
    std::vector<double> x = solveWithGaussianElimination(A, b);


    printResults(grid, matrix, b);
    // Вывод решения
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