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

int _u(float x, float y) {
    return (x + y);
}

void createMatrix(const std::vector<Node>& nodes, int minX, int maxX, int minY, int maxY) {
    for (size_t i = 0; i < nodes.size(); ++i) {
        int connections = 0;
        std::vector<int> neighbors;

        // Проверка соседей в пределах 2 единиц по x и y
        for (size_t j = 0; j < nodes.size(); ++j) {
            if (i != j) { // Пропустить сам узел
                // Проверка, что узел не DUMMY
                if (nodes[i].type != DUMMY && nodes[j].type != DUMMY) {
                    // Проверка соседей по горизонтали
                    if (nodes[i].y == nodes[j].y && std::abs(nodes[i].x - nodes[j].x) <= 2.0) {
                        // Проверка на границы
                        if (nodes[j].x >= minX && nodes[j].x <= maxX && nodes[j].y >= minY && nodes[j].y <= maxY) {
                            connections++;
                            neighbors.push_back(j);
                        }
                    }
                    // Проверка соседей по вертикали
                    else if (nodes[i].x == nodes[j].x && std::abs(nodes[i].y - nodes[j].y) <= 2.0) {
                        // Проверка на границы
                        if (nodes[j].x >= minX && nodes[j].x <= maxX && nodes[j].y >= minY && nodes[j].y <= maxY) {
                            connections++;
                            neighbors.push_back(j);
                        }
                    }
                }
            }
        }

        // Вывод информации о текущем узле, его координатах и соседях
        std::cout << "Узел " << i
            << " (тип: " << (nodes[i].type == INTERNAL ? "INTERNAL" : (nodes[i].type == BOUNDARY ? "BOUNDARY" : "DUMMY"))
            << ", координаты: (" << nodes[i].x << ", " << nodes[i].y << ")) "
            << "имеет " << connections << " соседей: ";
        for (int neighbor : neighbors) {
            std::cout << neighbor << " (координаты: (" << nodes[neighbor].x << ", " << nodes[neighbor].y << ")) ";
        }
        std::cout << std::endl;
    }
}


void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
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
}int main() {
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

    std::cout << "Узлы сетки в T-образной области и фиктивные узлы:\n";
    for (const auto& node : nodes) {
        std::cout << "(" << node.x << ", " << node.y << ")-" << static_cast<int>(node.type) << "\n";
    }

    // Создание 5-ти диагональной матрицы
    std::vector<std::vector<double>> matrix;
    int size;
    createMatrix(nodes,0,12,0,8);

    // Сохранение матрицы в файл
    saveMatrixToFile(matrix, "matrix.txt");

    std::cout << "5-ти диагональная матрица сохранена в файл matrix.txt." << std::endl;

    return 0;
}
