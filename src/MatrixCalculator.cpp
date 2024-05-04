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

/**
 * @brief Finds the maximum off-diagonal element of a matrix.
 * 
 * @tparam T Type of elements in the matrix.
 * @param matrix The matrix.
 * @param p Reference to store the row index of the maximum off-diagonal element.
 * @param q Reference to store the column index of the maximum off-diagonal element.
 * @param maxOffDiagonal Reference to store the value of the maximum off-diagonal element.
 */
template<typename T>
void MatrixCalculator::findMaxOffDiagonal(const std::vector<std::vector<T>>& matrix, int& p, int& q, T& maxOffDiagonal) {
    int n matrix.size();
    maxOffDiagonal = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(matrix[i][j]) > maxOffDiagonal) {
                maxOffDiagonal = std::abs(matrix[i][j]);
                p  = i;
                q = j;
            }
        }
    }
}

/**
 * @brief Computes eigenvalues and eigenvectors of a square matrix using the Jacobi eigenvalue algorithm.
 * 
 * @tparam T Type of elements in the matrix.
 * @param matrix The input matrix.
 * @return A pair containing eigenvalues and corresponding eigenvectors.
 * @throw std::invalid_argument If the input matrix is not square.
 */
template<typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> MatrixCalculator::jacobiEigen(const std::vector<std::vector<T>>& matrix) {
    int n = matrix.size();

    //checks for square matrix:
    if (n != matrix[0].size()) {
        throw std::invalid_argument("Matrix must be square.");
    }

    //initialize the eigen vector matrix as an identity matrix:
    std::vector<std::vector<T>> eigenVectors(n, std::vector<T>(n, 0));
    for (int i = 0; i < n; ++i) {
        eigenVectors[i][i] = 1;
    }

    //initialize the tolerance and maximum iterations:
    const T tolerance = 1e-10;
    const T toleranceChange = 0.1 * tolerance // Maximum 10% of tolerance change to consider the algo converged 
    const int maxIter = 1000;

    //initialize the variables:
    int p, q;
    T maxOffDiagonal, prevMaxOffDiagonal, theta, t, c, s;

    //perform jacobi iterations:
    int iterations = 0;
    while(true) {
        //find max diagonal element:
        findMaxOffDiagonal(matrix, p, q, maxOffDiagonal);

        // Check for convergence:
        if (maxOffDiagonal < tolerance || std::abs(maxOffDiagonal - prevMaxOffDiagonal) < toleranceChange || iterations > maxIter) {
            break;
        }

        // Store the current max off-diagonal element for the next iteration:
        prevMaxOffDiagonal = maxOffDiagonal;

        //compute the rotation angle:
        theta = 0.5 * std::atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p]);

        //sine and cosine of the rotation angle:
        s = std::sin(theta);
        c = std::cos(theta);

        //compute elements of rotation matrix;
        T tau = s/(c + 1);
        t = s/(1 + c);
        c = 1/std::sqrt(1 + t * t);
        s = t * c;

        //perform similarity transformation:
        for (int i = 0; i < n; ++i) {
            if (i != p && i != q) {
                T a_ip = matrix[i][p];
                T a_iq = matrix[i][q];
                matrix[i][p] = c * a_ip - s * a_iq;
                matrix[p][i] = matrix[i][p];
                matrix[i][q] = c * a_iq - s * a_ip;
                matrix[q][i] = matrix[i][q];
            }
        }

        T a_pp = matrix[p][p];
        T a_qq = matrix[q][q];
        matrix[p][p] = (c * c) * a_pp - 2 * (c * s) * matrix[p][q] + (s * s) * a_qq;
        matrix[q][q] = (s * s) * a_pp + 2 * (c * s) * matrix[p][q] + (c * c) * a_qq;
        matrix[p][q] = matrix[q][p] = 0;

        //update the eigen vectors:
        for (int i = 0; i < n; ++i) {
            T v_ip = eigenVectors[i][p];
            T v_iq = eigenVectors[i][p];
            eigenVectors[i][p] = c * v_ip - s * v_iq;
            eigenVectors[i][q] = c * v_iq + s * v_ip;
        }
        iterations++;
    }
    //extract the eigenvalues:
    std::vector<T> eigenvalues(n);
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = matrix[i][i];
    }

    return std::make_pair(eigenvalues, eigenVectors);
}

/**
 * @brief Computes the eigenvalues and eigenvectors of a square matrix using the Jacobi eigenvalue algorithm.
 * 
 * This method calculates the eigenvalues and eigenvectors of a square matrix using the Jacobi eigenvalue algorithm,
 * which iteratively diagonalizes the input matrix by performing similarity transformations until convergence.
 * The algorithm returns a pair containing the eigenvalues matrix and the eigenvectors matrix.
 * 
 * @param A The square matrix for which eigenvalues and eigenvectors are to be computed.
 * @return A pair containing the eigenvalues matrix and the eigenvectors matrix.
 * @throw std::invalid_argument If the input matrix is not square.
 */
std::pair<Matrix, Matrix> MatrixCalculator::eigen(const Matrix& A) {
    auto result = jacobiEigen(A.getData());
    Matrix eigenvalues({result.first});
    Matrix eigenVectors(result.second);
    return std::make_pair(eigenvalues, eigenVectors);
}

int main() {
    return 0;
}