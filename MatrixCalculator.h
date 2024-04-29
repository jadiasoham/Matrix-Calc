#ifndef MATRIX_CALCULATOR_H
#define MATRIX_CALCULATOR_H

#include "Matrix.h"

class MatrixCalculator {
public:
    Matrix addMatrices(const Matrix& A, const Matrix& B);
    Matrix subtractMatrices(const Matrix& A, const Matrix& B);
    Matrix elementWiseProduct(const Matrix& A, const Matrix& B);
    Matrix dotProduct(const Matrix& A, const Matrix& B);
};

#endif