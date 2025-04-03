#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <clocale>
#include <iomanip>

using namespace std;

const int PRECISION = 6;
const double GAMMA = 1.0; // Коэффициент для матрицы массы

// Базовые функции для МКЭ
double default_lambda(double x, double t) { return 1.0 + x * t; }
double default_f(double x, double t) { return x * t; }

struct BoundaryCondition {
    double point1;
    double point2;
    double value;
    vector<double> additionalPoints;
};

struct GridData {
    vector<double> xPoints;
    vector<double> tPoints;
    vector<double> hSteps;
    vector<double> tauSteps;
    vector<double> hAvgSteps;
    vector<double> tauAvgSteps;
    vector<BoundaryCondition> boundaryConditions;
    bool isUniformX;
    bool isUniformT;

    double (*lambda)(double, double) = default_lambda;
    double (*f)(double, double) = default_f;
};

struct FEMatrix {
    vector<vector<double>> M; // Матрица массы
    vector<vector<double>> G; // Матрица жесткости
    vector<double> b;         // Вектор нагрузки
    vector<double> solution;  // Вектор решения
};

void generateUniformGrid(vector<double>& points, vector<double>& steps,
    vector<double>& avgSteps, double min, double max,
    int pointsCount) {
    points.clear(); steps.clear(); avgSteps.clear();

    double step = (max - min) / (pointsCount - 1);
    for (int i = 0; i < pointsCount; ++i) {
        points.push_back(min + i * step);
    }

    for (size_t i = 1; i < points.size(); ++i) {
        steps.push_back(points[i] - points[i - 1]);
    }

    for (size_t i = 1; i < steps.size(); ++i) {
        avgSteps.push_back((steps[i - 1] + steps[i]) / 2.0);
    }
}

void generateNonUniformGrid(vector<double>& points, vector<double>& steps,
    vector<double>& avgSteps, double min, double max,
    int pointsCount) {
    points.clear(); steps.clear(); avgSteps.clear();

    double total = 0;
    for (int i = 1; i < pointsCount; ++i) total += 1.0 / i;

    double current = min;
    points.push_back(current);
    for (int i = 1; i < pointsCount; ++i) {
        current += (max - min) / (i * total);
        points.push_back(current);
    }

    for (size_t i = 1; i < points.size(); ++i) {
        steps.push_back(points[i] - points[i - 1]);
    }

    for (size_t i = 1; i < steps.size(); ++i) {
        avgSteps.push_back((steps[i - 1] + steps[i]) / 2.0);
    }
}

GridData readAndGenerateGrid(const string& filename, bool uniformX, bool uniformT) {
    GridData data;
    data.isUniformX = uniformX;
    data.isUniformT = uniformT;

    ifstream file(filename);
    if (!file) {
        cerr << "Ошибка открытия файла " << filename << endl;
        exit(1);
    }

    string line;

    // Чтение параметров сетки
    getline(file, line);
    istringstream xStream(line);
    double xMin, xMax;
    int xPoints;
    xStream >> xMin >> xMax >> xPoints;

    if (uniformX) {
        generateUniformGrid(data.xPoints, data.hSteps, data.hAvgSteps, xMin, xMax, xPoints);
    }
    else {
        generateNonUniformGrid(data.xPoints, data.hSteps, data.hAvgSteps, xMin, xMax, xPoints);
    }

    getline(file, line);
    istringstream tStream(line);
    double tMin, tMax;
    int tPoints;
    tStream >> tMin >> tMax >> tPoints;

    if (uniformT) {
        generateUniformGrid(data.tPoints, data.tauSteps, data.tauAvgSteps, tMin, tMax, tPoints);
    }
    else {
        generateNonUniformGrid(data.tPoints, data.tauSteps, data.tauAvgSteps, tMin, tMax, tPoints);
    }

    // Чтение граничных условий
    while (getline(file, line)) {
        istringstream bcStream(line);
        BoundaryCondition bc;
        bcStream >> bc.point1 >> bc.point2 >> bc.value;

        double point;
        while (bcStream >> point) bc.additionalPoints.push_back(point);

        data.boundaryConditions.push_back(bc);
    }

    return data;
}

vector<vector<double>> calculateMassMatrix(double h) {
    vector<vector<double>> M(2, vector<double>(2));
    double coef = GAMMA * h / 6.0;

    M[0][0] = 2.0 * coef;
    M[0][1] = 1.0 * coef;
    M[1][0] = 1.0 * coef;
    M[1][1] = 2.0 * coef;

    return M;
}

vector<vector<double>> calculateStiffnessMatrix(double lambda1, double lambda2, double h) {
    vector<vector<double>> G(2, vector<double>(2));
    double coef = (lambda1 + lambda2) / (2.0 * h);

    G[0][0] = 1.0 * coef;
    G[0][1] = -1.0 * coef;
    G[1][0] = -1.0 * coef;
    G[1][1] = 1.0 * coef;

    return G;
}

vector<double> calculateLoadVector(double f1, double f2, double h) {
    vector<double> b(2);

    b[0] = h * (2.0 * f1 + f2) / 6.0;
    b[1] = h * (f1 + 2.0 * f2) / 6.0;

    return b;
}

// Простое решение системы (для демонстрации)
vector<double> solveSystem(const vector<vector<double>>& K, const vector<double>& F) {
    int n = F.size();
    vector<double> solution(n, 0.0);

    // Здесь должна быть реальная процедура решения
    // Для демонстрации просто копируем вектор нагрузки
    solution = F;

    return solution;
}

FEMatrix assembleGlobalSystem(const GridData& grid) {
    FEMatrix global;
    int n = grid.xPoints.size();

    global.M.resize(n, vector<double>(n, 0.0));
    global.G.resize(n, vector<double>(n, 0.0));
    global.b.resize(n, 0.0);

    // Сборка глобальных матриц
    for (int i = 0; i < n - 1; i++) {
        double h = grid.hSteps[i];
        double x1 = grid.xPoints[i];
        double x2 = grid.xPoints[i + 1];

        double f1 = grid.f(x1, 0);
        double f2 = grid.f(x2, 0);
        double lambda1 = grid.lambda(x1, 0);
        double lambda2 = grid.lambda(x2, 0);

        auto localM = calculateMassMatrix(h);
        auto localG = calculateStiffnessMatrix(lambda1, lambda2, h);
        auto localb = calculateLoadVector(f1, f2, h);

        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                global.M[i + j][i + k] += localM[j][k];
                global.G[i + j][i + k] += localG[j][k];
            }
            global.b[i + j] += localb[j];
        }
    }

    // Решение системы (K = G + M)
    vector<vector<double>> K(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            K[i][j] = global.G[i][j] + global.M[i][j];
        }
    }

    global.solution = solveSystem(K, global.b);

    return global;
}

void printMatrix(const string& name, const vector<vector<double>>& matrix) {
    cout << "\n" << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << setw(10) << fixed << setprecision(PRECISION) << val << " ";
        }
        cout << endl;
    }
}

void printVector(const string& name, const vector<double>& vec) {
    cout << "\n" << name << ":\n";
    for (double val : vec) {
        cout << fixed << setprecision(PRECISION) << val << endl;
    }
}

void printLocalMatrices(const GridData& grid, const FEMatrix& fem) {
    cout << "\n=== ЛОКАЛЬНЫЕ МАТРИЦЫ ===";

    for (int i = 0; i < grid.xPoints.size() - 1; i++) {
        double h = grid.hSteps[i];
        double x1 = grid.xPoints[i];
        double x2 = grid.xPoints[i + 1];

        cout << "\n\nЭлемент " << i + 1 << " (узлы " << i << ";" << i + 1 << "):";
        cout << "\nДлина: " << h;
        cout << "\nX1 = " << x1 << ", X2 = " << x2;

        double f1 = grid.f(x1, 0);
        double f2 = grid.f(x2, 0);
        double lambda1 = grid.lambda(x1, 0);
        double lambda2 = grid.lambda(x2, 0);

        auto M = calculateMassMatrix(h);
        auto G = calculateStiffnessMatrix(lambda1, lambda2, h);
        auto b = calculateLoadVector(f1, f2, h);

        printMatrix("\nЛокальная матрица массы", M);
        printMatrix("Локальная матрица жесткости", G);
        printVector("Локальный вектор нагрузки", b);
    }
}

void printGlobalSystem(const FEMatrix& fem) {
    cout << "\n=== ГЛОБАЛЬНАЯ СИСТЕМА ===";

    // Глобальная матрица системы (K = G + M)
    vector<vector<double>> K(fem.M.size(), vector<double>(fem.M.size()));
    for (size_t i = 0; i < fem.M.size(); i++) {
        for (size_t j = 0; j < fem.M.size(); j++) {
            K[i][j] = fem.G[i][j] + fem.M[i][j];
        }
    }

    printMatrix("\nГлобальная матрица массы", fem.M);
    printMatrix("Глобальная матрица жесткости", fem.G);
    printMatrix("Глобальная матрица системы (K = G + M)", K);
    printVector("Глобальный вектор нагрузки", fem.b);
    printVector("Вектор решения", fem.solution);
}

int main() {
    setlocale(LC_ALL, "Russian");

    string filename = "input.txt";
    bool uniformX = false;
    bool uniformT = true;

    cout << "Чтение данных и генерация сетки...\n";
    GridData grid = readAndGenerateGrid(filename, uniformX, uniformT);

    cout << "\n=== ПАРАМЕТРЫ СЕТКИ ===";
    cout << "\nКоличество точек по X: " << grid.xPoints.size();
    cout << "\nКоличество точек по T: " << grid.tPoints.size();
    cout << "\nТип сетки X: " << (grid.isUniformX ? "РАВНОМЕРНАЯ" : "НЕРАВНОМЕРНАЯ");
    cout << "\nТип сетки T: " << (grid.isUniformT ? "РАВНОМЕРНАЯ" : "НЕРАВНОМЕРНАЯ") << endl;

    cout << "\nРасчет МКЭ системы...";
    FEMatrix fem = assembleGlobalSystem(grid);

    // Вывод всех результатов
    printLocalMatrices(grid, fem);
    printGlobalSystem(fem);

    return 0;
}