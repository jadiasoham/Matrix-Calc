#ifndef MATRIX_CALCULATOR_H
#define MATRIX_CALCULATOR_H

#include "Matrix.h"

class MatrixCalculator {
public:
    Matrix addMatrices(const Matrix& A, const Matrix& B);
    Matrix subtractMatrices(const Matrix& A, const Matrix& B);
    Matrix elementWiseProduct(const Matrix& A, const Matrix& B);
    Matrix dotProduct(const Matrix& A, const Matrix& B);
    std::pair<Matrix, Matrix> eigen(const Matrix& A);

private:
    template<typename T>
    void findMaxOffDiagonal(const std::vector<std::vector<T>>& matrix, int& p, int& q, T& maxOffDiagonal);

    template<typename T>
    std::pair<std::vector<T>, std::vector<std::vector<T>>> jacobiEigen(const std::vector<std::vector<T>>& matrix);
};

#endif
