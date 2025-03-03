#include <vector>
#include <iostream>
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

    // Проверка на горизонтальную прямую (4, 0) - (8, 0)
    if (std::abs(y - 0.0) < EPSILON && x >= 4 && x <= 8) return true;

    // Проверка на вертикальную прямую (8, 0) - (8, 4)
    if (std::abs(x - 8.0) < EPSILON && y >= 0 && y <= 4) return true;

    // Проверка на горизонтальную прямую (8, 4) - (12, 4)
    if (std::abs(y - 4.0) < EPSILON && x >= 8 && x <= 12) return true;

    // Проверка на вертикальную прямую (12, 4) - (12, 8)
    if (std::abs(x - 12.0) < EPSILON && y >= 4 && y <= 8) return true;

    // Проверка на горизонтальную прямую (12, 8) - (0, 8)
    if (std::abs(y - 8.0) < EPSILON && x >= 0 && x <= 12) return true;

    // Проверка на вертикальную прямую (0, 8) - (0, 4)
    if (std::abs(x - 0.0) < EPSILON && y >= 4 && y <= 8) return true;

    // Проверка на горизонтальную прямую (0, 4) - (4, 4)
    if (std::abs(y - 4.0) < EPSILON && x >= 0 && x <= 4) return true;

    // Проверка на вертикальную прямую (4, 4) - (4, 0)
    if (std::abs(x - 4.0) < EPSILON && y >= 0 && y <= 4) return true;

    return false; // Если ни одна проверка не прошла
}

// Проверка, находится ли точка внутри T-образной области
bool isInsideTShape(double x, double y) {
    // Внутренние узлы ножки (x от 4 до 8, y от 0 до 4)
    if (x > 4 && x < 8 && y > 0 && y < 8) return true;

    // Внутренние узлы шляпки (x от 0 до 12, y от 4 до 8, исключая границы)
    if (x > 0 && x < 12 && y > 4 && y < 8) return true;

    return false;
}

int main() {
    setlocale(LC_ALL, "Russian");

    double h = 2.0; // Шаг сетки
    double left = 0, right = 12, bottom = 0, top = 12; // Границы области

    std::vector<Node> nodes;

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

    return 0;
}