#include <iostream>
#include "MatrixCalculator.h"

std::vector<std::vector<int>> addMatrices(const Matrix& A, const Matrix& B) {
    // Check if matrices are of compatible dimensions or are empty
    if (A.getRows() != B.getRows() || A.getRows() != 0 || B.getRows() != 0) {
        throw std::runtime_error("Matrices are of incompatible dimensions and hence cannot be added.");
    }
    
    // Add matrices
    std::vector<std::vector<int>> result(A.getRows(), std::vector<int>(A.getCols(), 0));
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < A.getCols(); ++j) {
            result[i][j] = A.getData()[i][j] + B.getData()[i][j];
        }
    }
    return result;
}

std::vector<std::vector<int>> subtractMatrices(const Matrix& A, const Matrix& B) {
    //add code here
}

std::vector<std::vector<int>> matrixProduct(const Matrix& A, const Matrix& B) {
    //add code here
}

std::vector<std::vector<int>> dotProduct(const Matrix& A, const Matrix& B) {
    //add code here
}

double determinant(const Matrix& A) {
    //add code here
}