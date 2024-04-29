#ifndef MATRIX_CALCULATOR_H
#define MATRIX_CALCULATOR_H

#include "Matrix.h"

class MatrixCalculator {
public:
    std::vector<std::vector<int>> addMatrices(const Matrix& A, const Matrix& B);
    std::vector<std::vector<int>> subtractMatrices(const Matrix& A, const Matrix& B);
    std::vector<std::vector<int>> matrixProduct(const Matrix& A, const Matrix& B);
    std::vector<std::vector<int>> dotProduct(const Matrix& A, const Matrix& B);
    double determinant(const Matrix& A);
};

#endif