#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <clocale>
#include <iomanip> // для точного вывода чисел

using namespace std;

// Точность вывода
const int PRECISION = 6;

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
};

// Функция для генерации равномерной сетки с правильным расчетом шагов
void generateUniformGrid(vector<double>& points, vector<double>& steps, vector<double>& avgSteps,
    double min, double max, int pointsCount) {
    points.clear();
    steps.clear();
    avgSteps.clear();

    if (pointsCount < 2) {
        cerr << "Ошибка: количество точек должно быть не менее 2" << endl;
        return;
    }

    double step = (max - min) / (pointsCount - 1);
    for (int i = 0; i < pointsCount; ++i) {
        points.push_back(min + i * step);
    }

    // Расчет шагов
    for (size_t i = 1; i < points.size(); ++i) {
        steps.push_back(points[i] - points[i - 1]);
    }

    // Расчет средних шагов (для внутренних точек)
    for (size_t i = 1; i < steps.size(); ++i) {
        avgSteps.push_back((steps[i - 1] + steps[i]) / 2.0);
    }
}

// Функция для генерации неравномерной сетки
void generateNonUniformGrid(vector<double>& points, vector<double>& steps, vector<double>& avgSteps,
    double min, double max, int pointsCount) {
    points.clear();
    steps.clear();
    avgSteps.clear();

    if (pointsCount < 2) {
        cerr << "Ошибка: количество точек должно быть не менее 2" << endl;
        return;
    }

    // Неравномерная сетка сгущается к началу
    double total = 0;
    for (int i = 1; i < pointsCount; ++i) {
        total += 1.0 / i;
    }

    double current = min;
    points.push_back(current);
    for (int i = 1; i < pointsCount; ++i) {
        current += (max - min) / (i * total);
        points.push_back(current);
    }

    // Расчет шагов
    for (size_t i = 1; i < points.size(); ++i) {
        steps.push_back(points[i] - points[i - 1]);
    }

    // Расчет средних шагов
    for (size_t i = 1; i < steps.size(); ++i) {
        avgSteps.push_back((steps[i - 1] + steps[i]) / 2.0);
    }
}

GridData readAndGenerateGrid(const string& filename, bool uniformX, bool uniformT) {
    GridData data;
    data.isUniformX = uniformX;
    data.isUniformT = uniformT;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось открыть файл " << filename << endl;
        exit(1);
    }

    string line;

    // Чтение параметров для X
    getline(file, line);
    double xMin, xMax;
    int xPoints;
    istringstream xStream(line);
    xStream >> xMin >> xMax >> xPoints;

    if (uniformX) {
        generateUniformGrid(data.xPoints, data.hSteps, data.hAvgSteps, xMin, xMax, xPoints);
    }
    else {
        generateNonUniformGrid(data.xPoints, data.hSteps, data.hAvgSteps, xMin, xMax, xPoints);
    }

    // Чтение параметров для T
    getline(file, line);
    double tMin, tMax;
    int tPoints;
    istringstream tStream(line);
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

        double additionalPoint;
        while (bcStream >> additionalPoint) {
            bc.additionalPoints.push_back(additionalPoint);
        }
        data.boundaryConditions.push_back(bc);
    }

    return data;
}

void printGridInfo(const GridData& grid) {
    cout << fixed << setprecision(PRECISION);

    cout << "\nТекущие настройки:";
    cout << "\nСетка по X: " << (grid.isUniformX ? "РАВНОМЕРНАЯ" : "НЕРАВНОМЕРНАЯ");
    cout << "\nСетка по T: " << (grid.isUniformT ? "РАВНОМЕРНАЯ" : "НЕРАВНОМЕРНАЯ") << endl;

    cout << "\nСетка по X:";
    for (size_t i = 0; i < grid.xPoints.size(); ++i) {
        cout << "\nx[" << i << "] = " << grid.xPoints[i];
        if (i > 0) cout << ", h_" << i << " = " << grid.hSteps[i - 1];
        if (i > 0 && i - 1 < grid.hAvgSteps.size())
            cout << ", h̄_" << i << " = " << grid.hAvgSteps[i - 1];
    }

    cout << "\n\nСетка по T:";
    for (size_t i = 0; i < grid.tPoints.size(); ++i) {
        cout << "\nt[" << i << "] = " << grid.tPoints[i];
        if (i > 0) cout << ", τ_" << i << " = " << grid.tauSteps[i - 1];
        if (i > 0 && i - 1 < grid.tauAvgSteps.size())
            cout << ", τ̄_" << i << " = " << grid.tauAvgSteps[i - 1];
    }
    cout << endl;
}

int main() {
    setlocale(LC_ALL, "Russian");

    string filename = "input.txt";
    bool uniformX = false;  // По умолчанию равномерная сетка по X
    bool uniformT = true;  // По умолчанию равномерная сетка по T

    cout << "Чтение данных и генерация сетки...";
    GridData grid = readAndGenerateGrid(filename, uniformX, uniformT);

    printGridInfo(grid);

    // Тестирование функций
    double testX = 2.5, testT = 1.5;
    cout << "\nТестирование функций в точке (" << testX << ", " << testT << "):";
    cout << "\nu(x,t) = " << testX * testT;
    cout << "\nf(x,t) = " << testX * testT;
    cout << "\nlambda(x,t) = " << testX * testT;
    cout << "\ngamma(x,t) = " << testX * testT << endl;

    // Вывод граничных условий
    if (!grid.boundaryConditions.empty()) {
        cout << "\nГраничные условия:";
        for (const auto& bc : grid.boundaryConditions) {
            cout << "\nТочки " << bc.point1 << " и " << bc.point2 << " → значение " << bc.value;
            if (!bc.additionalPoints.empty()) {
                cout << " (Доп. точки:";
                for (double pt : bc.additionalPoints) cout << " " << pt;
                cout << ")";
            }
        }
        cout << endl;
    }

    return 0;
}