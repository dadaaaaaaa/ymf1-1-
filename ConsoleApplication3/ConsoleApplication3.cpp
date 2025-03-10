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
double hx, hy, h; // Шаг сетки
double left = 0, right = 12, bottom = 0, top = 12; // Границы области
enum NodeType {
    DUMMY = 0,    // Фиктивные узлы
    BOUNDARY = 1, // Граничные узлы
    INTERNAL = 2  // Внутренние узлы
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

constexpr double EPSILON = 1e-6;
int shotfile() {
    std::ofstream file(filename);
    file >> left >> right>>hx;
    file >> bottom >> top>>hy;
    hx =(abs(left) + abs(right))/ hx;
    hy = (abs(bottom) + abs(top)) / hy;
}
bool isOnBoundary(double x, double y) {
    const double EPSILON = 1e-9; // Погрешность для проверки

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

void createMatrix(const std::vector<Node>& nodes, std::vector<std::vector<int>>& matrix) {
    size_t size = nodes.size();

    // Инициализация матрицы нулями
    matrix.resize(size, std::vector<int>(size, 0));

    for (size_t i = 0; i < size; ++i) {
        if (nodes[i].type != DUMMY) {
            // Узел связан сам с собой
            matrix[i][i] = 1;

            // Проверка соседей в пределах 2 единиц по x и y
            for (size_t j = 0; j < size; ++j) {
                if (i != j && nodes[j].type != DUMMY) {
                    // Проверка на горизонтальную связь
                    if (nodes[i].y == nodes[j].y && std::abs(nodes[i].x - nodes[j].x) <= 2.0) {
                        matrix[i][j] = -1/4; // Связь с соседом
                    }
                    // Проверка на вертикальную связь
                    else if (nodes[i].x == nodes[j].x && std::abs(nodes[i].y - nodes[j].y) <= 2.0) {
                        matrix[i][j] = -1/4; // Связь с соседом
                    }
                }
            }
        }
    }
}

void saveMatrixToFile(const std::vector<std::vector<int>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : matrix) {
            for (const auto& value : row) {
                file << value << " ";
            }
            file << "\n";
        }
        file.close();
    }
    else {
        std::cerr << "Ошибка открытия файла!" << std::endl;
    }
}

int countOnes(const std::vector<std::vector<int>>& matrix) {
    int count = 0;
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            if (value == 1) {
                count++;
            }
        }
    }
    return count;
}

void printNodes(const std::vector<Node>& nodes) {
    std::cout << "Список узлов с координатами:" << std::endl;
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "Узел " << i + 1 << " (" << nodes[i].x << ", " << nodes[i].y << "), тип: ";
        switch (nodes[i].type) {
        case DUMMY:
            std::cout << "Фиктивный" << std::endl;
            break;
        case BOUNDARY:
            std::cout << "Граничный" << std::endl;
            break;
        case INTERNAL:
            std::cout << "Внутренний" << std::endl;
            break;
        }
    }
}

    for (const auto& row : matrix) {
        for (double val : row) {
            file << val << " ";
        }
        file << std::endl;
    }

    shotfile();

// Функция для решения СЛАУ методом Гаусса
std::vector<double> solveWithGaussianElimination(const std::vector<std::vector<double>>& matrix, const std::vector<double>& b) {
    int N = matrix.size();
    std::vector<double> x(N, 0.0);
    std::vector<std::vector<double>> augmentedMatrix(N, std::vector<double>(N + 1));

    // Генерация узлов
    for (double x = left; x <= right; x += hx) {
        for (double y = bottom; y <= top; y += hy) {
            Node node;
            node.x = x;
            node.y = y;

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

    // Создание матрицы смежности
    std::vector<std::vector<int>> adjacencyMatrix;
    createMatrix(nodes, adjacencyMatrix);

    // Сохранение матрицы в файл
    saveMatrixToFile(adjacencyMatrix, "matrix.txt");

    // Подсчет единиц в матрице
    int onesCount = countOnes(adjacencyMatrix);
    std::cout << "Количество единиц в матрице смежности: " << onesCount << std::endl;

    // Вывод результатов
    printResults(grid, matrix, b);

    std::cout << "Матрица смежности успешно сохранена в файл 'matrix.txt'." << std::endl;

    return 0;
}
