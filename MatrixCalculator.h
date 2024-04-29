#ifndef MATRIX_CALCULATOR_H
#define MATRIX_CALCULATOR_H

#include "Matrix.h"

using namespace std;

class MatrixCalculator {
public:
    vector<vector<int>> addMatrices(const Matrix& A, const Matrix& B);
    vector<vector<int>> subtractMatrices(const Matrix& A, const Matrix& B);
    vector<vector<int>> multiplyMatrices(const Matrix& A, const Matrix& B);
    double determinant(const Matrix& A);
};

#endif