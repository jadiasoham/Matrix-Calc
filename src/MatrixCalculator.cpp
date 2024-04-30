#include <iostream>
#include "MatrixCalculator.h"

/**
 * @brief Adds two matrices element-wise.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of adding the matrices.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
Matrix MatrixCalculator::addMatrices(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
        throw std::runtime_error("Matrices are of incompatible dimensions and hence cannot be added.");
    }
    
    std::vector<std::vector<int>> result(A.getRows(), std::vector<int>(A.getCols(), 0));
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < A.getCols(); ++j) {
            result[i][j] = A.getData()[i][j] + B.getData()[i][j];
        }
    }
    return Matrix(result);
}

/**
 * @brief Subtracts one matrix from another element-wise.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of subtracting matrix B from matrix A.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
Matrix MatrixCalculator::subtractMatrices(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
        throw std::runtime_error("Matrices are of incompatible dimenesions and hence cannot be subtracted.");
    }

    std::vector<std::vector<int>> result(A.getRows(), std::vector<int>(A.getCols(), 0));
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < A.getCols(); ++j) {
            result[i][j] = A.getData()[i][j] - B.getData()[i][j];
        }
    }
    return Matrix(result);
}

/**
 * @brief Computes the element-wise product of two matrices.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of multiplying the matrices element-wise.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
Matrix MatrixCalculator::elementWiseProduct(const Matrix& A, const Matrix& B) {
    // Check compatibility: number of columns of A must be equal to number of rows of B
    if (A.getCols() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
        throw std::runtime_error("Matrices are not of compatible dimension and hence cannot be multiplied.");
    }

    // Initialize result matrix with appropriate dimensions
    std::vector<std::vector<int>> result(A.getRows(), std::vector<int>(B.getCols(), 0));

    // Compute matrix product
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < B.getCols(); ++j) {
            for (size_t k = 0; k < A.getCols(); ++k) {
                result[i][j] += A.getData()[i][k] * B.getData()[k][j];
            }
        }
    }
    return Matrix(result);
}

/**
 * @brief Computes the dot product of two matrices.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of the dot product operation.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
Matrix MatrixCalculator::dotProduct(const Matrix& A, const Matrix& B) {
    // Check compatibility:
    if (A.getCols() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
        throw std::runtime_error("Matrices are not of compatible dimension and hence cannot be multiplied.");
    }

    std::vector<std::vector<int>> result(A.getCols(), std::vector<int>(B.getCols(), 0));
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < B.getCols(); ++j) {
            result[i][j] = 0;
            for (size_t k = 0; k < B.getRows(); ++k) {
                result[i][j] += A.getData()[i][k] * B.getData()[k][j];
            }
        }
    }
    return Matrix(result);
}

int main() {
    return 0;
}
