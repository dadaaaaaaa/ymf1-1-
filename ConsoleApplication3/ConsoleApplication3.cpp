#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

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
                        matrix[i][j] = 1; // Связь с соседом
                    }
                    // Проверка на вертикальную связь
                    else if (nodes[i].x == nodes[j].x && std::abs(nodes[i].y - nodes[j].y) <= 2.0) {
                        matrix[i][j] = 1; // Связь с соседом
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
    std::vector<std::vector<int>> adjacencyMatrix;
    createMatrix(nodes, adjacencyMatrix);

    // Сохранение матрицы в файл
    saveMatrixToFile(adjacencyMatrix, "matrix.txt");

    // Подсчет единиц в матрице
    int onesCount = countOnes(adjacencyMatrix);
    std::cout << "Количество единиц в матрице смежности: " << onesCount << std::endl;

    // Вывод списка узлов
    printNodes(nodes);

    std::cout << "Матрица смежности успешно сохранена в файл 'matrix.txt'." << std::endl;

    return 0;
}
