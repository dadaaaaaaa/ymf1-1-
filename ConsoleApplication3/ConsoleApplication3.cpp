
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
    Point p1, p2;
    int type;
};

// Структура для хранения узла сетки
struct GridNode {
    double x, y;
    std::string type;
    int boundaryType;
    int index;
};
using namespace std;

// Метод сопряженных градиентов для пятидиагональной матрицы
void MCG(const vector<double>& di, const vector<double>& u1, const vector<double>& u2,
	const vector<double>& l1, const vector<double>& l2, const vector<double>& b,
	vector<double>& x, int max_iter = 1000, double eps = 1e-6) {

	int n = di.size();
	vector<double> r(n), p(n), Ap(n);
	double rr, rr_new, alpha, beta;

	// Инициализация
	x.assign(n, 0.0);
	r = b;

	// Учет граничных условий
	for (int i = 0; i < n; ++i) {
		if (abs(di[i] - 1.0) < 1e-9) { // Граничный узел
			x[i] = b[i];
			r[i] = 0.0;
		}
	}

	// Первое вычисление невязки
	rr = 0.0;
	for (int i = 0; i < n; ++i) rr += r[i] * r[i];

	// Главный цикл
	for (int iter = 0; iter < max_iter; ++iter) {
		// Умножение матрицы на вектор
		for (int i = 0; i < n; ++i) {
			Ap[i] = di[i] * p[i];
			if (i > 0) Ap[i] += l1[i - 1] * p[i - 1];
			if (i < n - 1) Ap[i] += u1[i] * p[i + 1];
			if (i >= 5) Ap[i] += l2[i - 5] * p[i - 5];
			if (i < n - 5) Ap[i] += u2[i] * p[i + 5];
		}

		// Вычисление коэффициентов
		double pAp = 0.0;
		for (int i = 0; i < n; ++i) pAp += p[i] * Ap[i];
		alpha = rr / pAp;

		// Обновление решения и невязки
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		// Проверка сходимости
		rr_new = 0.0;
		for (int i = 0; i < n; ++i) rr_new += r[i] * r[i];
		if (sqrt(rr_new) < eps) {
			cout << "Сходимость достигнута за " << iter << " итераций" << endl;
			return;
		}

		// Обновление направления
		beta = rr_new / rr;
		for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
		rr = rr_new;
	}
	cout << "Достигнуто максимальное число итераций" << endl;
}
bool isPointOnLineSegment(const Point& p1, const Point& p2, const Point& p) {
    double cross = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
    if (std::abs(cross) > 1e-12) return false;

    if (p1.x != p2.x) {
        if ((p.x < std::min(p1.x, p2.x)) || (p.x > std::max(p1.x, p2.x))) return false;
    }
    else {
        if ((p.y < std::min(p1.y, p2.y)) || (p.y > std::max(p1.y, p2.y))) return false;
    }

    double dot = (p.x - p1.x) * (p2.x - p1.x) + (p.y - p1.y) * (p2.y - p1.y);
    if (dot < 0) return false;

    double squaredLength = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
    if (dot > squaredLength) return false;
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

void printResults(const std::vector<std::vector<GridNode>>& grid, const std::vector<double>& b) {
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
    std::cout << "Правая часть b:" << std::endl;
    for (int i = 0; i < b.size(); i++) {
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
    }
}

// Функции для метода Гаусса-Зейделя
double norm(int n, vector<double> mas) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += mas[i] * mas[i];
    return sqrt(sum);
}

double norm2(int n, vector<double> v, vector<double> mas) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += (v[i] - mas[i]) * (v[i] - mas[i]);
    return sqrt(sum);
}

double iter(int i, int n, int m, double w, vector<double> di, vector<double> u1, vector<double> u2, vector<double> l1, vector<double> l2, vector<double> f, vector<double> x0, vector<double>& xnext) {
    double sum = 0, a = 0;
    for (int j = 0; j < n; j++) {
        if (i == j) a = di[i]; // главная диагональ
        else {
            int r = j - i;
            if (r == m) a = u1[i]; // верхние диагонали
            else if (r == 1) a = u2[i];
            else if (r == -m) a = l1[j]; // нижние диагонали
            else if (r == -1) a = l2[j];
            else a = 0;
        }
        sum += a * x0[j];
    }
    xnext[i] = x0[i] + w / di[i] * (f[i] - sum);
    return sum;
}

void zeidel(int n, int m, int maxit, double e, double w, vector<double> di, vector<double> u1, vector<double> u2, vector<double> l1, vector<double> l2, vector<double> f, vector<double>& x0) {
    ofstream out;
    out.open("output.txt");
    int flag = 1, k = 0;
    vector<double> Ax(n, 0.0);
    while (flag) {
        for (int i = 0; i < n; i++)
            Ax[i] = iter(i, n, m, w, di, u1, u2, l1, l2, f, x0, x0);
        k++;
        if ((norm2(n, f, Ax) / norm(n, f) > e && k < maxit))
            flag = 1;
        else flag = 0;
    }
    out << "Метод Гаусса-Зейделя:\n";
        for (int i = 0; i < n; i++)
            out << scientific << setprecision(15) << x0[i] << "\n";
        out.close();
}

void buildMatrixAndRightHandSides(const std::vector<std::vector<GridNode>>& grid,
    const std::vector<double>& x_nodes, const std::vector<double>& y_nodes,
    std::vector<double>& di, std::vector<double>& u1, std::vector<double>& u2,
    std::vector<double>& l1, std::vector<double>& l2, std::vector<double>& b,
    std::vector<double>& x0) {
    int N = grid.size() * grid[0].size();
    di.assign(N, 0.0);
    u1.assign(N, 0.0);
    u2.assign(N, 0.0);
    l1.assign(N, 0.0);
    l2.assign(N, 0.0);
    b.assign(N, 0.0);
    x0.assign(N, 0.0);

    int rows = grid.size();
    int cols = grid[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int idx = i * cols + j;

            if (grid[i][j].type == "фиктивный") {
                di[idx] = 1.0;
                b[idx] = 0;
                x0[idx] = 0;
                continue;
            }

            if (grid[i][j].type == "граничный" && grid[i][j].boundaryType == 1) {
                di[idx] = 1.0;
                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
                x0[idx] = b[idx]; // Устанавливаем начальное приближение для граничных узлов
            }
            else if (grid[i][j].type == "внутренний") {
                double hx = (j > 0) ? (x_nodes[j] - x_nodes[j - 1]) : 1.0;
                double hy = (i > 0) ? (y_nodes[i] - y_nodes[i - 1]) : 1.0;

                b[idx] = boundaryFunction(grid[i][j].x, grid[i][j].y);
                di[idx] = _cu(hx, hy);

                if (j > 0 && grid[i][j - 1].type != "фиктивный") {
                    l1[idx - 1] = _yu(hx);
                }
                if (i > 0 && grid[i - 1][j].type != "фиктивный") {
                    l2[idx - cols] = _yu(hy);
                }
                if (j < cols - 1 && grid[i][j + 1].type != "фиктивный") {
                    u1[idx] = _yu(hx);
                }
                if (i < rows - 1 && grid[i + 1][j].type != "фиктивный") {
                    u2[idx] = _yu(hy);
                }
            }
        }
    }
}


void writeResultsToFile(const std::string& filename, const std::vector<std::vector<GridNode>>& grid, const std::vector<double>& solutions, const std::vector<double>& b) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Ошибка открытия файла для записи!" << std::endl;
        return;
    }

    outFile << "N \tx\ty\tu\tu*\t|u - u*|\n";
    for (const auto& row : grid) {
        for (const auto& node : row) {
            int idx = node.index;
            double difference = std::abs(solutions[idx] - b[idx]);
            outFile << node.index << "\t" << node.x << "\t" << node.y << "\t"
                << b[idx] << "\t" << solutions[idx] << "\t" << difference << "\n";
        }
    }

    outFile.close();
    std::cout << "Результаты успешно сохранены в файл: " << filename << std::endl;
}
int main() {
	std::setlocale(LC_ALL, "Russian");

    double x_min, x_max, y_min, y_max;
    int x_steps, y_steps;
    std::vector<BoundaryCondition> boundaryConditions;
    readInputData("input.txt", x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions);
    std::vector<Point> polygon;
    for (const auto& condition : boundaryConditions) {
        polygon.push_back(condition.p1);
        polygon.push_back(condition.p2);
    }

    std::vector<std::vector<GridNode>> grid;
    std::vector<double> x_nodes, y_nodes;
    std::tie(grid, x_nodes, y_nodes) = createGrid(x_min, x_max, x_steps, y_min, y_max, y_steps, boundaryConditions, polygon);

    std::vector<double> di, u1, u2, l1, l2, b, x0;
    buildMatrixAndRightHandSides(grid, x_nodes, y_nodes, di, u1, u2, l1, l2, b, x0);
    int maxit = 1000;
    double e = 1e-15, w = 1.1;
    zeidel(di.size(), grid[0].size(), maxit, e, w, di, u1, u2, l1, l2, b, x0);

    writeResultsToFile("results.txt", grid, x0, b);

    printResults(grid, b);

    std::cout << "Решение СЛАУ:" << std::endl;
    for (int i = 0; i < x0.size(); ++i) {
        std::cout << "x[" << i << "] = " << x0[i] << std::endl;
    }
    std::cout << "Решение Ошибка:" << std::endl;
    for (int i = 0; i < x0.size(); ++i) {
        std::cout << "x[" << i << "] = " << x0[i] - b[i] << std::endl;
    }

	return 0;
}