#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <cmath>
#include <locale>
#include <algorithm>
#include <tuple> 
#include <iomanip>

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
    // Проверка на горизонтальную и вертикальную прямые
    if (std::abs(y - 0.0) < EPSILON && x >= 4 && x <= 8) return true;
    if (std::abs(x - 8.0) < EPSILON && y >= 0 && y <= 4) return true;
    if (std::abs(y - 4.0) < EPSILON && x >= 8 && x <= 12) return true;
    if (std::abs(x - 12.0) < EPSILON && y >= 4 && y <= 8) return true;
    if (std::abs(y - 8.0) < EPSILON && x >= 0 && x <= 12) return true;
    if (std::abs(x - 0.0) < EPSILON && y >= 4 && y <= 8) return true;
    if (std::abs(y - 4.0) < EPSILON && x >= 0 && x <= 4) return true;
    if (std::abs(x - 4.0) < EPSILON && y >= 0 && y <= 4) return true;

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

            int idx = grid[i][j].index;

            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                A[idx][idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = x_nodes[j] - x_nodes[j - 1];
                double hy = y_nodes[i] - y_nodes[i - 1];

                A[idx][idx] = _cu(hx, hy);

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

    // Создание сетки и получение x_nodes и y_nodes
    std::vector<std::vector<GridNode>> grid;
    std::vector<double> x_nodes, y_nodes;
    std::tie(grid, x_nodes, y_nodes) = createGrid(x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions, polygon);

    // Построение матрицы и правой части
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    buildMatrixAndRightHandSide(grid, x_nodes, y_nodes, A, b);

    // Вывод результатов
    printResults(grid, A, b);

    return 0;
}

