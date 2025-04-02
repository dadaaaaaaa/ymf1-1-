#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

struct Matrix {
    double* di;    // ������������ ��������
    double* gg;    // �������� ������(�������) ���������
    double* F;     // ������ ������ �����
    int* jg;       // ������ � ��������� ������ �����(��������) � ������� gg 
    int* ig;       // ������ ��������(�����) ��������� � ������� gg
    int N;         // ����������� �������
    double normF;  // ����� ������� F
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
