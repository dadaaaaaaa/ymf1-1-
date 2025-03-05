#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

enum NodeType {
    DUMMY = 0,    // Фиктивные узлы
    BOUNDARY = 1, // Граничные узлы
    INTERNAL = 2  // Внутренние узлы
};

struct Node {
    double x, y;
    NodeType type;
};

constexpr double EPSILON = 1e-6;

bool isOnBoundary(double x, double y) {
    const double EPSILON = 1e-9; // Погрешность для проверки

    // Проверка на горизонтальную и вертикальную прямые
    if (std::abs(y - 0.0) < EPSILON && x >= 4 && x <= 8) return true;
    if (std::abs(x - 8.0) < EPSILON && y >= 0 && y <= 4) return true;
    if (std::abs(y - 4.0) < EPSILON && x >= 8 && x <= 12) return true;
    if (std::abs(x - 12.0) < EPSILON && y >= 4 && y <= 8) return true;
    if (std::abs(y - 8.0) < EPSILON && x >= 0 && x <= 12) return true;
    if (std::abs(x - 0.0) < EPSILON && y >= 4 && y <= 8) return true;
    if (std::abs(y - 4.0) < EPSILON && x >= 0 && x <= 4) return true;
    if (std::abs(x - 4.0) < EPSILON && y >= 0 && y <= 4) return true;

    return false; // Если ни одна проверка не прошла
}

bool isInsideTShape(double x, double y) {
    if (x > 4 && x < 8 && y > 0 && y < 8) return true; // Ножка
    if (x > 0 && x < 12 && y > 4 && y < 8) return true; // Шляпка
    return false;
}

// Функция для расчета значения функции u(x, y)
double _u(double x, double y) {
    return (x + y); // Пример функции
}

// Функция для поиска точки в списке узлов
const Node* findNode(double x, double y, const std::vector<Node>& nodes) {
    for (const auto& node : nodes) {
        if (std::abs(node.x - x) < EPSILON && std::abs(node.y - y) < EPSILON) {
            return &node;
        }
    }
    return nullptr; // Точка не найдена
}

// Функции для аппроксимации производных
double _uLx(const Node& node, double h, const std::vector<Node>& nodes) {
    if (node.type == DUMMY) return 0.0; // Если узел фиктивный, возвращаем 0
    double x = node.x + h;
    double y = node.y;
    const Node* neighbor = findNode(x, y, nodes);
    if (!neighbor || neighbor->type == DUMMY) return (_u(node.x, node.y) / h); // Если точка не найдена или фиктивная
    return (_u(x, y) - _u(node.x, node.y)) / h;
}

double _uRx(const Node& node, double h, const std::vector<Node>& nodes) {
    if (node.type == DUMMY) return 0.0; // Если узел фиктивный, возвращаем 0
    double x = node.x - h;
    double y = node.y;
    const Node* neighbor = findNode(x, y, nodes);
    if (!neighbor || neighbor->type == DUMMY) return (_u(node.x, node.y) / h); // Если точка не найдена или фиктивная
    return (_u(node.x, node.y) - _u(x, y)) / h;
}

double _uCx(const Node& node, double h, const std::vector<Node>& nodes) {
    if (node.type == DUMMY) return 0.0; // Если узел фиктивный, возвращаем 0
    double x1 = node.x + h;
    double x2 = node.x - h;
    double y = node.y;
    const Node* neighbor1 = findNode(x1, y, nodes);
    const Node* neighbor2 = findNode(x2, y, nodes);
    if (!neighbor1 || neighbor1->type == DUMMY) {
        return (_u(x1, y) / (2 * h));
    }
    if (!neighbor2 || neighbor2->type == DUMMY) {
        return (_u(x2, y)) / (2 * h);
    }
    return (_u(x1, y) - _u(x2, y)) / (2 * h);
}

double _uLy(const Node& node, double h, const std::vector<Node>& nodes) {
    if (node.type == DUMMY) return 0.0; // Если узел фиктивный, возвращаем 0
    double x = node.x;
    double y = node.y + h;
    const Node* neighbor = findNode(x, y, nodes);
    if (!neighbor || neighbor->type == DUMMY) return (_u(node.x, node.y) / h); // Если точка не найдена или фиктивная
    return (_u(x, y) - _u(node.x, node.y)) / h;
}

double _uRy(const Node& node, double h, const std::vector<Node>& nodes) {
    if (node.type == DUMMY) return 0.0; // Если узел фиктивный, возвращаем 0
    double x = node.x;
    double y = node.y - h;
    const Node* neighbor = findNode(x, y, nodes);
    if (!neighbor || neighbor->type == DUMMY) return (_u(node.x, node.y) / h); // Если точка не найдена или фиктивная
    return (_u(node.x, node.y) - _u(x, y)) / h;
}

void createMatrix(const std::vector<Node>& nodes, std::vector<std::vector<double>>& matrix) {
    size_t size = nodes.size();
    matrix.resize(size, std::vector<double>(size, 0.0));

    for (size_t i = 0; i < size; ++i) {
        if (nodes[i].type != DUMMY) {
            matrix[i][i] = _uCx(nodes[i], 2, nodes);

            for (size_t j = 0; j < size; ++j) {
                if (i != j && nodes[j].type != DUMMY) {
                    // Проверка на горизонтальную связь
                    if (nodes[i].y == nodes[j].y && std::abs(nodes[i].x - nodes[j].x) <= 2.0) {
                        if (nodes[i].x < nodes[j].x) {
                            matrix[i][j] = _uRx(nodes[i], 2, nodes); // Правый сосед
                        }
                        else {
                            matrix[i][j] = _uLx(nodes[i], 2, nodes); // Левый сосед
                        }
                    }
                    // Проверка на вертикальную связь
                    else if (nodes[i].x == nodes[j].x && std::abs(nodes[i].y - nodes[j].y) <= 2.0) {
                        if (nodes[i].y < nodes[j].y) {
                            matrix[i][j] = _uLy(nodes[i], 2, nodes); // Верхний сосед
                        }
                        else {
                            matrix[i][j] = _uRy(nodes[i], 2, nodes); // Нижний сосед
                        }
                    }
                }
            }
        }
    }
}

void convertTo5Diagonal(const std::vector<std::vector<double>>& matrix, std::vector<double>& mainDiagonal, std::vector<double>& upperDiagonal, std::vector<double>& lowerDiagonal) {
    size_t size = matrix.size();
    mainDiagonal.resize(size);
    upperDiagonal.resize(size - 1);
    lowerDiagonal.resize(size - 1);

    for (size_t i = 0; i < size; ++i) {
        mainDiagonal[i] = matrix[i][i];
        if (i < size - 1) {
            upperDiagonal[i] = matrix[i][i + 1];
            lowerDiagonal[i] = matrix[i + 1][i];
        }
    }
}

void save5DiagonalToFile(const std::vector<double>& mainDiagonal, const std::vector<double>& upperDiagonal, const std::vector<double>& lowerDiagonal, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        // Сохраняем главную диагональ
        for (const auto& value : mainDiagonal) {
            file << std::fixed << std::setprecision(6) << value << " ";
        }
        file << "\n";

        // Сохраняем верхнюю диагональ
        for (const auto& value : upperDiagonal) {
            file << std::fixed << std::setprecision(6) << value << " ";
        }
        file << "\n";

        // Сохраняем нижнюю диагональ
        for (const auto& value : lowerDiagonal) {
            file << std::fixed << std::setprecision(6) << value << " ";
        }
        file << "\n";

        file.close();
    }
    else {
        std::cerr << "Ошибка открытия файла!" << std::endl;
    }
}
void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : matrix) {
            for (const auto& value : row) {
                file << std::fixed << std::setprecision(6) << value << " ";
            }
            file << "\n";
        }
        file.close();
    }
    else {
        std::cerr << "Ошибка открытия файла для сохранения матрицы!" << std::endl;
    }
}
void printNodes(const std::vector<Node>& nodes) {
    for (const auto& node : nodes) {
        std::cout << "Node: (" << node.x << ", " << node.y << "), Type: " << node.type << std::endl;
    }
}
int countNonZero(const std::vector<std::vector<double>>& matrix) {
    int count = 0;
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            if (std::abs(value) > EPSILON) { // Проверяем, что элемент ненулевой
                count++;
            }
        }
    }
    return count;
}

int main() {
    setlocale(LC_ALL, "Russian");

    double h = 2.0; // Шаг сетки
    double left = 0, right = 12, bottom = 0, top = 12; // Границы области

    std::vector<Node> nodes;

    // Генерация узлов
    for (double x = left; x <= right; x += h) {
        for (double y = bottom; y <= top; y += h) {
            Node node;
            node.x = x;
            node.y = y;

            if (isOnBoundary(x, y)) {
                node.type = BOUNDARY; // Узел на границе
            }
            else if (isInsideTShape(x, y)) {
                node.type = INTERNAL; // Внутренний узел
            }
            else {
                node.type = DUMMY; // Фиктивный узел
            }

            nodes.push_back(node);
        }
    }

    // Создание матрицы смежности
    std::vector<std::vector<double>> adjacencyMatrix;
    createMatrix(nodes, adjacencyMatrix);

    // Сохранение полной матрицы в файл
    saveMatrixToFile(adjacencyMatrix, "matrix.txt");

    // Конвертация в 5-ти диагональную матрицу
    std::vector<double> mainDiagonal, upperDiagonal, lowerDiagonal;
    convertTo5Diagonal(adjacencyMatrix, mainDiagonal, upperDiagonal, lowerDiagonal);

    // Сохранение 5-ти диагональной матрицы в файл
    save5DiagonalToFile(mainDiagonal, upperDiagonal, lowerDiagonal, "5diagonal_matrix.txt");

    // Подсчет ненулевых элементов в матрице
    int nonZeroCount = countNonZero(adjacencyMatrix);
    std::cout << "Количество ненулевых элементов в матрице смежности: " << nonZeroCount << std::endl;

    // Вывод списка узлов
    printNodes(nodes);

    std::cout << "Матрица смежности успешно сохранена в файл 'matrix.txt'." << std::endl;
    std::cout << "5-ти диагональная матрица успешно сохранена в файл '5diagonal_matrix.txt'." << std::endl;

    return 0;
}

