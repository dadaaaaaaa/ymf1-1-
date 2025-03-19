#include<stdio.h>
#include<vector>
#include<iostream>
#include<fstream>
#include <iomanip>
using namespace std;
double lambda = 1, gamma = 1;
int o = 8, p = 8;

struct area
{
	vector <double> knotes;	// Координаты оси
	vector <int> counts;		// Количество узлов между соответствующими координатами
	vector <double> coefs;	// Коэфициент q
	vector <double> h;		// Шаги
	vector <double> coords; 	// Координаты
	int kol_uz;			// Количество всех узлов
} areaX, areaY;

double func(double x, double y)
{
	return x + y;
		//x * x + 2 * y * y - 6;
		//0.5 * pow(x, 3) + pow(y, 3) + 3 * x + 6 * y;
		//-12*pow(x,2) - 12*pow(y, 2) + pow(x, 4) + pow(y, 4);		
		//-12*x*x + pow(x, 4);
		//2 * cos(x);

}

double funcU(double x, double y)
{
	return x + y;
		//x * x + 2 * y * y;
		//0.5 * pow(x, 3) + pow(y, 3);
		//pow(x, 4) + pow(y, 4);
		//pow(x, 4);
		//cos(x);
}

void input(area& area, ifstream& name)
{
	int n;
	name >> n;
	area.knotes.resize(n);
	area.counts.resize(n - 1);
	area.coefs.resize(n - 1);

	for (int i = 0; i < n; i++)
		name >> area.knotes[i];

	for (int i = 0; i < n - 1; i++)
		name >> area.counts[i];

	for (int i = 0; i < n - 1; i++)
		name >> area.coefs[i];

	double sum = 0;
	for (int i = 0; i < n - 1; i++)
		sum += area.counts[i];
	area.h.resize(sum);
	area.coords.resize(sum + 1);
	name.close();
}

void h(area& ar)
{
	double h0;
	int k = 0;
	for (int i = 0; i < ar.counts.size(); i++)
	{
		if (ar.coefs[i] != 1)
			h0 = (ar.knotes[i + 1] - ar.knotes[i]) * (1 - ar.coefs[i]) / (1 - pow(ar.coefs[i], ar.counts[i])); // Неравномерный шаг
		else h0 = (ar.knotes[i + 1] - ar.knotes[i]) / ar.counts[i]; // Равномерный шаг

		ar.h[k] = h0;
		k++;
		double h = h0;
		for (int j = 0; j < ar.counts[i] - 1; j++, k++)
		{
			h *= ar.coefs[i];
			ar.h[k] = h;
		}
	}
}

void coords(area& ar)
{
	ar.coords[0] = ar.knotes[0];
	for (int i = 1; i < ar.h.size() + 1; i++)
		ar.coords[i] = ar.coords[i - 1] + ar.h[i - 1];
	ar.kol_uz = ar.coords.size();
}

void null(vector <double>& di, int n)
{
	di.resize(n);
	for (int i = 0; i < di.size(); i++)
		di[i] = 0;
}

void search_neighbors(area arX, area arY, int t, int i, int j, vector <double>& u1, vector <double>& u2, vector <double>& l1, vector <double>& l2)
{
	//Для равномерного шага
	l2[t - 1] = -lambda / pow(arX.h[j - 1], 2);	// Левая точка
	u2[t] = -lambda / pow(arX.h[j], 2);	// Правая точка
	l1[t - arX.kol_uz] = -lambda / pow(arY.h[i - 1], 2);	// Нижняя точка
	u1[t] = -lambda / pow(arY.h[i], 2);		// Верхняя точка

//	//Для неравномерного шага
//	double h2 = arX.h[j - 1] + arX.h[j];
//	l2[t - 1] = -2 * lambda / (arX.h[j - 1] * h2);	// Левая точка
//	u2[t] = -2 * lambda / (arX.h[j] * h2);	// Правая точка
//	h2 = arY.h[i - 1] + arY.h[i];
//	l1[t - arX.kol_uz] = -2 * lambda / (arY.h[i - 1] * h2);	// Нижняя точка
//	u1[t] = -2 * lambda / (arY.h[i] * h2);		// Верхняя точка
}

void fictiv(area arX, area arY, vector <bool>& fict, vector <double>& di, vector <double>& f)
{
	double x0, x1, y0, y1;
	fict.resize(arX.kol_uz * arY.kol_uz);
	ifstream in;
	in.open("Fict.txt");
	in >> x0 >> x1 >> y0 >> y1;
	in.close();

	// Центр сетки по оси X
	double centerX = (arX.knotes[0] + arX.knotes[arX.knotes.size() - 1]) / 2.0;

	for (int i = 0; i < arY.kol_uz; i++)
	{
		for (int j = 0; j < arX.kol_uz; j++)
		{
			int t = i * arX.kol_uz + j;

			// Горизонтальная линия (верхняя часть)
			bool isHorizontal = (arY.coords[i] >= y0 && arY.coords[i] <= y1);

			// Вертикальная линия (центр по оси X)
			bool isVertical = (arX.coords[j] >= centerX - x0 && arX.coords[j] <= centerX + x0);

			// Если узел принадлежит букве "Т", он не фиктивный
			if ((isHorizontal && arX.coords[j] >= centerX - x1 && arX.coords[j] <= centerX + x1) ||
				(isVertical && arY.coords[i] <= y1))
			{
				fict[t] = false;
				di[t] = 1;
				f[t] = 0;
			}
			else
			{
				fict[t] = true;
			}
		}
	}
}

void build_matrix(area arX, area arY, vector <bool> fict, vector <double>& di, vector <double>& u1, vector <double>& u2, vector <double>& l1, vector <double>& l2, vector <double>& f)
{
	for (int i = 0; i < arY.kol_uz; i++)
		for (int j = 0; j < arX.kol_uz; j++)
		{
			int t = i * arX.kol_uz + j;
			if (fict[t]) //если нефиктивная точка
			{
				if (i == 0 || j == 0 || (i == arY.kol_uz - 1) || (j == arX.kol_uz - 1) || (j == o && i <= p) || (i == p && j >= o)) //если граничная точка
				{
					di[t] = 1;
					f[t] = funcU(arX.coords[j], arY.coords[i]);
				}
				else /*if (di[t] == 0)*/
				{
					di[t] = 2 * lambda / (arX.h[j - 1] * arX.h[j]) + 2 * lambda / (arY.h[i - 1] * arY.h[i]) + gamma;
					f[t] = func(arX.coords[j], arY.coords[i]);
					search_neighbors(arX, arY, t, i, j, u1, u2, l1, l2);
				}
			}
			else //если фиктивная точка
			{
				di[t] = 1;
				f[t] = 0;
			}
		}
}

double norm(int n, vector <double> mas)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += mas[i] * mas[i];
	return sqrt(sum);
}

double norm2(int n, vector <double> v, vector <double> mas)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += (v[i] - mas[i]) * (v[i] - mas[i]);
	return sqrt(sum);
}

double iter(int i, int n, int m, double w, vector <double> di, vector <double> u1, vector <double> u2, vector <double> l1, vector <double> l2, vector <double> f, vector <double> x0, vector <double>& xnext)
{
	double sum = 0, a = 0;
	for (int j = 0; j < n; j++)
	{
		if (i == j) a = di[i]; // главная диагональ
		else
		{
			int r = j - i;
			if (r == m) a = u1[i];				//верхние диагонали
			else if (r == 1) a = u2[i];
			else if (r == -m) a = l1[j];			//нижние диагонали
			else if (r == -1) a = l2[j];
			else a = 0;
		}
		sum += a * x0[j];
	}
	xnext[i] = x0[i] + w / di[i] * (f[i] - sum);
	return sum;
}

void zeidel(int n, int m, int maxit, double e, double w, vector <double> di, vector <double> u1, vector <double> u2, vector <double> l1, vector <double> l2, vector <double> f, vector <double>& x0, area arX, area arY)
{
	ofstream out;
	out.open("output.txt");
	int flag = 1, k = 0;
	vector <double> Ax;
	null(Ax, n);
	while (flag)
	{
		for (int i = 0; i < n; i++)
			Ax[i] = iter(i, n, m, w, di, u1, u2, l1, l2, f, x0, x0);
		k++;
		if ((norm2(n, f, Ax) / norm(n, f) > e) && k < maxit)
			flag = 1;
		else flag = 0;
	}

	// Заголовок таблицы
	out << "N\tx\ty\tu\tu*\t|u - u*|\n";

	for (int i = 0; i < n; i++)
	{
		// Вычисляем координаты узла
		int row = i / arX.kol_uz; // Строка (ось Y)
		int col = i % arX.kol_uz; // Столбец (ось X)
		double x = arX.coords[col];
		double y = arY.coords[row];

		// Вычисляем точное значение функции (u*)
		double u_star = funcU(x, y);

		// Вычисляем разницу между вычисленным и точным значением
		double diff = abs(x0[i] - u_star);

		// Вывод в формате N  x  y  u  u*  |u - u*|
		out << i << "\t" << x << "\t" << y << "\t"
			<< scientific << setprecision(15) << x0[i] << "\t"
			<< u_star << "\t" << diff << "\n";
	}
	out.close();
}
void print_full_matrix(int n, int m, vector<double>& di, vector<double>& u1, vector<double>& u2, vector<double>& l1, vector<double>& l2)
{
	ofstream out;
	out.open("matrix.txt");
	out << "Полная матрица системы:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double a = 0;
			if (i == j) a = di[i]; // Главная диагональ
			else if (j == i + 1) a = u2[i]; // Первая верхняя диагональ
			else if (j == i + m) a = u1[i]; // Вторая верхняя диагональ
			else if (j == i - 1) a = l2[j]; // Первая нижняя диагональ
			else if (j == i - m) a = l1[j]; // Вторая нижняя диагональ
			out << scientific << setprecision(6) << a << "\t";
		}
		out << "\n";
	}
	out.close();
}
void print_non_fictitious_points(area arX, area arY, vector<bool>& fict) {
	cout << "Не фиктивные точки:" << endl;
	for (int i = 0; i < arY.kol_uz; i++) {
		for (int j = 0; j < arX.kol_uz; j++) {
			int t = i * arX.kol_uz + j; // Номер узла
			if (!fict[t]) { // Если точка не фиктивная
				double x = arX.coords[j];
				double y = arY.coords[i];

				// Определение типа узла
				string nodeType;
				if (i == 0 || j == 0 || (i == arY.kol_uz - 1) || (j == arX.kol_uz - 1)) {
					nodeType = "0"; // Граничный узел
				}
				else {
					nodeType = "1"; // Внутренний узел
				}

				cout << "Номер узла: " << t << " Координаты: (" << x << ", " << y << ") Тип: " << nodeType << endl;
			}
		}
	}
}
int main()
{
	int maxit = 1000;
	double e = 1e-15, w = 1.1;
	vector <double> di, u1, u2, l1, l2, f, x0;
	vector <bool> fict, kraev;
	ifstream inX, inY;
	inX.open("AreaX.txt");
	input(areaX, inX);
	inY.open("AreaY.txt");
	input(areaY, inY);
	h(areaX);
	h(areaY);
	coords(areaX);
	coords(areaY);

	null(di, areaX.kol_uz * areaY.kol_uz);
	null(f, areaX.kol_uz * areaY.kol_uz);
	null(u2, areaX.kol_uz * areaY.kol_uz - 1);	// Правые точки
	null(l2, areaX.kol_uz * areaY.kol_uz - 1);	// Левые точки
	null(u1, areaX.kol_uz * areaY.kol_uz - areaX.kol_uz);	// Верхние точки
	null(l1, areaX.kol_uz * areaY.kol_uz - areaX.kol_uz);	// Нижние точки
	//kraev_us(areaX, areaY, kraev, di, f);
	fictiv(areaX, areaY, fict, di, f);
	//inside_area(areaX, areaY, fict, di, u1, u2, l1, l2, f);
	build_matrix(areaX, areaY, fict, di, u1, u2, l1, l2, f);

	// Вывод полной матрицы системы
	print_full_matrix(areaX.kol_uz * areaY.kol_uz, areaX.kol_uz, di, u1, u2, l1, l2);
	print_non_fictitious_points(areaX, areaY, fict);
	null(x0, areaX.kol_uz * areaY.kol_uz);
	zeidel(areaX.kol_uz * areaY.kol_uz, areaX.kol_uz, maxit, e, w, di, u1, u2, l1, l2, f, x0, areaX, areaY);
}