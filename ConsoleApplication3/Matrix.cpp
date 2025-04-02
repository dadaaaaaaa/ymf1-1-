#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

struct Matrix {
    double* di;    // диагональные элементы
    double* gg;    // элементы нижней(верхней) диагонали
    double* F;     // вектор правой части
    int* jg;       // массив с индексами начала строк(столбцов) в массиве gg 
    int* ig;       // номера столбцов(строк) элементов в массиве gg
    int N;         // размерность матрицы
    double normF;  // норма вектора F
};

extern Matrix A;
extern int Maxiter;
extern double Eps;

void readFromFile();
double norm(double* x, int size);
void multMV(double* x, double* r);
void MCG();
void MCGdiag();
void printToFile(double* x, int iter, double nevR);
double scMult(double* x, double* y);
void Hilbert();

#endif // MATRIX_H
