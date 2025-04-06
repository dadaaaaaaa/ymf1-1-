#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <clocale>
#include <iomanip>
#include <algorithm>

using namespace std;

const int PRECISION = 6;

double default_lambda(double x, double t) { return 1+ x * t; }
double default_f(double x, double t) { return 1 + x * t; }
double default_u(double x, double t) { return 1 +  x * t; }
double default_gamma(double x, double t) { return 1 +  x * t; }

struct BoundaryCondition {
    int node;       // Номер узла
    double value;   // Значение условия
    double time;    // Время действия условия (-1 для всех времен)
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
};

struct FEMatrix {
    vector<vector<double>> M; // Матрица массы
    vector<vector<double>> G; // Матрица жесткости
    vector<double> solution;  // Вектор решения
    vector<vector<double>> allF; // Все векторы правой части
    vector<vector<double>> allSolutions; // Все решения
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
        bcStream >> bc.node >> bc.value >> bc.time;
        data.boundaryConditions.push_back(bc);
    }

    return data;
}

vector<vector<double>> calculateMassMatrix(double h, double x, double t) {
    vector<vector<double>> M(2, vector<double>(2));
    double gamma_val = default_gamma(x, t);
    double coef = gamma_val * h / 6.0;

    M[0][0] = 2.0 * coef;
    M[0][1] = 1.0 * coef;
    M[1][0] = 1.0 * coef;
    M[1][1] = 2.0 * coef;

    return M;
}

vector<vector<double>> calculateStiffnessMatrix(double x1, double x2, double t, double h) {
    vector<vector<double>> G(2, vector<double>(2));
    double lambda1 = default_lambda(x1, t);
    double lambda2 = default_lambda(x2, t);
    double coef = (lambda1 + lambda2) / (2.0 * h);

    G[0][0] = 1.0 * coef;
    G[0][1] = -1.0 * coef;
    G[1][0] = -1.0 * coef;
    G[1][1] = 1.0 * coef;

    return G;
}


vector<double> calculateLoadVector(double x1, double x2, double t, double h) {
    vector<double> b(2);
    double f1 = default_f(x1, t);
    double f2 = default_f(x2, t);

    b[0] = default_u(x1, t);
    b[1] = default_u(x2, t);

    return b;
}


vector<double> solveSystem(vector<vector<double>>& K, vector<double>& F) {
    int n = F.size();
    vector<double> solution(n, 0.0);

    // Простейший решатель (заглушка)
    for (int i = 0; i < n; i++) {
        if (K[i][i] != 0) {
            solution[i] = F[i] / K[i][i];
        }
    }
    return solution;
}

void applyBoundaryConditions(vector<vector<double>>& K, vector<double>& F,
    const vector<BoundaryCondition>& bcs, double currentTime) {
    for (const auto& bc : bcs) {
        int node = bc.node;
        if (bc.time < 0 || fabs(bc.time - currentTime) < 1e-6) {
            // Обнуляем строку и столбец
            for (int k = 0; k < K.size(); k++) {
                K[node][k] = 0.0;
                K[k][node] = 0.0;
            }
            K[node][node] = 1.0;

            F[node] = bc.value;

            cout << "  Применено граничное условие: q[" << node << "] = "
                << bc.value << " (время " << currentTime << ")" << endl;
        }
    }
}

void printLocalMatrix(const string& name, const vector<vector<double>>& matrix) {
    cout << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << setw(10) << fixed << setprecision(PRECISION) << val << " ";
        }
        cout << endl;
    }
}

void printLocalVector(const string& name, const vector<double>& vec) {
    cout << name << ":\n";
    for (double val : vec) {
        cout << fixed << setprecision(PRECISION) << val << endl;
    }
}

void printGlobalMatrix(const string& name, const vector<vector<double>>& matrix, int showSize = 5) {
    cout << "\n" << name << " (первые " << showSize << "x" << showSize << " элементов):\n";
    int size = min(showSize, (int)matrix.size());
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << setw(10) << fixed << setprecision(3) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printGlobalVector(const string& name, const vector<double>& vec, int showSize = 5) {
    cout << name << " (первые " << showSize << " элементов):\n";
    int size = min(showSize, (int)vec.size());
    for (int i = 0; i < size; i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}

void printFullVector(const string& name, const vector<double>& vec) {
    cout << name << " (полный вектор):\n";
    for (size_t i = 0; i < vec.size(); i++) {
        cout << "[" << i << "] = " << vec[i] << endl;
    }
}
FEMatrix solveTimeDependentSystem(const GridData& grid, const vector<double>& q0) {
    FEMatrix global;
    int n = grid.xPoints.size();
    int timeSteps = grid.tPoints.size();

    global.M.resize(n, vector<double>(n, 0.0));
    global.G.resize(n, vector<double>(n, 0.0));

    vector<double> q_prev = q0;

    // Основной временной цикл
    for (int j = 1; j < timeSteps; j++) {
        double dt = grid.tauSteps[j - 1];
        double t = grid.tPoints[j-1]; // Текущее время из сетки
        double q1 = q_prev[j - 1];
        double q12 = q_prev[j];
        cout << "\n\n=== ВРЕМЕННОЙ ШАГ " << j << " ===";
        cout << "\nВремя t = " << t << ", Δt = " << dt << endl;

        // Обнуляем глобальные матрицы и вектора для текущего шага
        vector<vector<double>> globalM(n, vector<double>(n, 0.0));
        vector<vector<double>> globalG(n, vector<double>(n, 0.0));
        vector<double> b(n, 0.0);

        // Сборка глобальных матриц для текущего времени
        for (int i = 0; i < n - 1; i++) {
            double h = grid.hSteps[i];
            double x1 = grid.xPoints[i];
            double x2 = grid.xPoints[i + 1];

            // Вычисляем локальные матрицы для текущего времени t
            auto localM = calculateMassMatrix(h, x1, t);
            auto localG = calculateStiffnessMatrix(x1, x2, t, h);
            auto localb = calculateLoadVector(x1, x2, t, h);

            cout << "\nЭлемент " << i + 1 << " (узлы " << i << "-" << i + 1 << "):";
            cout << "\nДлина h = " << h << ", X1 = " << x1 << ", X2 = " << x2;
            printLocalMatrix("\nЛокальная матрица массы:", localM);
            printLocalMatrix("Локальная матрица жесткости:", localG);
            printLocalVector("Локальный вектор нагрузки:", localb);

            // Добавляем в глобальные матрицы
            for (int k = 0; k < 2; k++) {
                for (int m = 0; m < 2; m++) {
                    globalM[i + k][i + m] += localM[k][m];
                    globalG[i + k][i + m] += localG[k][m];
                }
                b[i + k] += localb[k];
            }
        }

        printGlobalMatrix("\nГлобальная матрица массы:", globalM);
        printGlobalMatrix("Глобальная матрица жесткости:", globalG);
        printGlobalVector("\nГлобальный вектор нагрузки:", b);

        // Формируем систему уравнений
        vector<vector<double>> K(n, vector<double>(n));
        vector<double> F(n);

        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                K[i][k] = (1.0 / dt) * globalM[i][k] + globalG[i][k];
            }
            F[i] = b[i];
            for (int k = 0; k < n; k++) {
                F[i] += (1.0 / dt) * globalM[i][k] * q_prev[k];
            }
        }

        global.allF.push_back(F);

        printGlobalMatrix("\nМатрица системы K:", K);
        printFullVector("Полный вектор правой части F:", F);

        applyBoundaryConditions(K, F, grid.boundaryConditions, t);

        vector<double> q = solveSystem(K, F);
        q_prev = q;
        global.allSolutions.push_back(q);

        cout << "\nРешение на шаге " << j << ":";
        printFullVector("Вектор решения q:", q);
    }

    global.solution = q_prev;
    return global;
}
int main() {
    setlocale(LC_ALL, "Russian");

    string filename = "input.txt";
    bool uniformX = true;
    bool uniformT = true;

    cout << "=== МЕТОД КОНЕЧНЫХ ЭЛЕМЕНТОВ ДЛЯ ВРЕМЕННОЙ ЗАДАЧИ ===" << endl;
    cout << "Чтение данных и генерация сетки...\n";
    GridData grid = readAndGenerateGrid(filename, uniformX, uniformT);

    cout << "\n=== ПАРАМЕТРЫ СЕТКИ ===" << endl;
    cout << "Пространственная сетка (X): " << grid.xPoints.size() << " узлов, "
        << (grid.isUniformX ? "равномерная" : "неравномерная") << endl;
    cout << "Временная сетка (T): " << grid.tPoints.size() << " шагов, "
        << (grid.isUniformT ? "равномерная" : "неравномерная") << endl;
    cout << "Граничные условий: " << grid.boundaryConditions.size() << endl;

    // Задание начального условия q0
    vector<double> q0(grid.xPoints.size(), 0.0);
    cout << "\n=== ЗАДАНИЕ НАЧАЛЬНОГО УСЛОВИЯ ===" << endl;
    cout << "Выберите способ задания q0:\n"
        << "1. Константное значение для всех узлов\n"
        << "2. Линейное распределение\n"
        << "3. Вручную для каждого узла\n"
        << "Ваш выбор: ";

    int choice;
    cin >> choice;

    switch (choice) {
    case 1: {
        double value;
        cout << "Введите значение для всех узлов: ";
        cin >> value;
        fill(q0.begin(), q0.end(), value);
        break;
    }
    case 2: {
        double q_start, q_end;
        cout << "Введите значение в первом узле: ";
        cin >> q_start;
        cout << "Введите значение в последнем узле: ";
        cin >> q_end;
        for (size_t i = 0; i < q0.size(); i++) {
            double t = static_cast<double>(i) / (q0.size() - 1);
            q0[i] = q_start + t * (q_end - q_start);
        }
        break;
    }
    case 3: {
        cout << "Введите значения для каждого узла:" << endl;
        for (size_t i = 0; i < q0.size(); i++) {
            cout << "Узел " << i << " (x=" << grid.xPoints[i] << "): ";
            cin >> q0[i];
        }
        break;
    }
    default:
        cout << "Неверный выбор, используется q0 = 0" << endl;
    }

    cout << "\n=== НАЧАЛЬНОЕ УСЛОВИЕ q0 ===" << endl;
    printFullVector("Вектор q0:", q0);

    cout << "\n=== НАЧАЛО РАСЧЕТА ===" << endl;
    FEMatrix fem = solveTimeDependentSystem(grid, q0);

    cout << "\n=== ИТОГОВОЕ РЕШЕНИЕ ===" << endl;
    printFullVector("Финальное решение:", fem.solution);

    // Вывод всех правых частей
    cout << "\n=== ВЕКТОРЫ ПРАВОЙ ЧАСТИ НА КАЖДОМ ШАГЕ ===" << endl;
    for (size_t i = 0; i < fem.allF.size(); i++) {
        cout << "\nШаг " << i + 1 << " (t = " << grid.tPoints[i + 1] << "):";
        printFullVector("Вектор F:", fem.allF[i]);
    }

    // Вывод всех решений
    cout << "\n=== РЕШЕНИЯ НА КАЖДОМ ШАГЕ ===" << endl;
    for (size_t i = 0; i < fem.allSolutions.size(); i++) {
        cout << "\nШаг " << i + 1 << " (t = " << grid.tPoints[i + 1] << "):";
        printFullVector("Решение q:", fem.allSolutions[i]);
    }

    cout << "\n=== РАСЧЕТ ЗАВЕРШЕН ===" << endl;
    return 0;
}