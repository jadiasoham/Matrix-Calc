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

/**
 * @brief Produces transpose of a matrix
 * 
 * @param matrix The matrix to be transposed
 * 
 * @return The transposed matrix
*/
Matrix MatrixCalculator::transposeMatrix(const Matrix& matrix) {
    size_t rows = matrix.getRows();
    size_t cols = matrix.getCols();

    // Initialize the transpose matrix with appropriate dimensions
    std::vector<std::vector<int>> transposedData(cols, std::vector<int>(rows));

    // Transpose the matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            transposedData[j][i] = matrix.getData()[i][j];
        }
    }
    return Matrix(transposedData);
}

/**
 * @brief Calculate the rank of a matrix
 * 
 * @param matrix The matrix whose rank is to be calculated
 * 
 * @return The rank of a matrix as an int
*/
int MatrixCalculator::rank(const Matrix& matrix) {
    // Transpose the matrix the bring it into the row echelon form
    Matrix transposed = transposeMatrix(matrix);

    size_t rows = transposed.getRows();
    size_t cols = transposed.getCols();

    size_t rank = 0;
    size_t col = 0;

    // Row reduction to bring it into row echelon form:
    for (size_t row = 0; row < rows; ++row) {
        // Find the first non zero element in the current column:
        while (col < cols && transposed.getData()[row][col] == 0) {
            ++col;
        }

        // If we find non zero elements, increment the rank and continue row reduction
        if (col < cols) {
            rank++;
            // Eliminate non zero elements below current pivot
            for (size_t i = row + 1; i < rows; ++i) {
                int factor = transposed.getData()[i][col] / transposed.getData()[row][col];
                for (size_t j = col; j < cols; ++j) {
                    transposed.getData()[i][j] -= factor * transposed.getData()[row][j];
                }
            }
        }
    }
    return rank;
}