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
std::vector<std::vector<int>> MatrixCalculator::addMatrices(const Matrix& A, const Matrix& B) {
    // Check if matrices are of compatible dimensions or are empty
    if (A.getRows() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
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

/**
 * @brief Subtracts one matrix from another element-wise.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of subtracting matrix B from matrix A.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
std::vector<std::vector<int>> MatrixCalculator::subtractMatrices(const Matrix& A, const Matrix& B) {
    if (A.getRows() != B.getRows() || A.getRows() == 0 || B.getRows() == 0) {
        throw std::runtime_error("Matrices are of incompatible dimenesions and hence cannot be subtracted.");
    }

    std::vector<std::vector<int>> result(A.getRows(), std::vector<int>(A.getCols(), 0));
    for (size_t i = 0; i < A.getRows(); ++i) {
        for (size_t j = 0; j < A.getCols(); ++j) {
            result[i][j] = A.getData()[i][j] - B.getData()[i][j];
        }
    }
    return result;
}

/**
 * @brief Computes the element-wise product of two matrices.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of multiplying the matrices element-wise.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
std::vector<std::vector<int>> MatrixCalculator::elementWiseProduct(const Matrix& A, const Matrix& B) {
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
    return result;
}

/**
 * @brief Computes the dot product of two matrices.
 * 
 * @param A First matrix.
 * @param B Second matrix.
 * @return The result of the dot product operation.
 * @throw std::runtime_error If matrices have incompatible dimensions or are empty.
 */
std::vector<std::vector<int>> MatrixCalculator::dotProduct(const Matrix& A, const Matrix& B) {
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
    return result;
}

/**
 * @brief Computes the determinant of a square matrix.
 * 
 * @param A The input square matrix.
 * @return The determinant of the matrix.
 * @throw std::runtime_error If matrix is not square or empty.
 */
double MatrixCalculator::determinant(const Matrix& A) {
    // Check if non empty square matrix
    if (A.getRows() != A.getCols() || A.getRows() == 0) {
        throw std::runtime_error("Deteminant can be calculate only for a non empty square matrix.");
    }

    int n = A.getRows();
    if (n == 1) {
        return A.getData()[0][0];
    }
    double det = 0.0;
    for (int i = 0; i < n; ++i) {
        std::vector<std::vector<int>> minor(n - 1, std::vector<int>(n - 1, 0));
        for (int j = 1; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (k < i) {
                    minor[j - 1][k] = A.getData()[j][k];
                } else if (k > i) {
                    minor[j - 1][k - 1] = A.getData()[j][k];
                }
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * A.getData()[0][i] * determinant(minor);
    }
    return det;
}

/**
 * @brief Displays the contents of a matrix.
 * 
 * @param A The matrix to be displayed.
 */
void MatrixCalculator::displayMatrix(const Matrix& A) {
    for (const auto& row : A.getData()) {
        for (int value : row) {
            std::cout<< value << "\t";
        }
        std::cout << std::endl;
    }
}