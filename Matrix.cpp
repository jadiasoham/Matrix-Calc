#include "Matrix.h"
#include <iostream>

/**
 * @brief Constructs a matrix from the given input data.
 * 
 * @param input The input data representing the matrix.
 */
Matrix::Matrix(const std::vector<std::vector<int>>& input) : data(input) {}
// Overload the constructor:
Matrix::Matrix(const std::vector<std::vector<double>>& input) {
    //Convert double data to int
    std::vector<std::vector<int>> intData(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        intData[i].resize(input.size());
        for (size_t j = 0; j < input[i].size(); ++j) {
            intData[i][j] = static_cast<int>(input[i][j]);
        }
    }
    data = intData;
}

/**
 * @brief Retrieves the data of the matrix.
 * 
 * @return The matrix data.
 */
std::vector<std::vector<int>> Matrix::getData() const {
    return data;
}

/**
 * @brief Retrieves the number of rows in the matrix.
 * 
 * @return The number of rows.
 */
int Matrix::getRows() const {
    return data.size();
}

/**
 * @brief Retrieves the number of columns in the matrix.
 * 
 * @return The number of columns.
 */
int Matrix::getCols() const {
    if (data.empty()) {
        return 0;
    }
    return data[0].size();
}

/**
 * @brief Computes the determinant of the matrix.
 * 
 * @return The determinant of the matrix.
 * @throw std::runtime_error If the matrix is not square or empty.
 */
double Matrix::determinant() {
    // Check if non-empty square matrix;
    if (data.size() != data[0].size() || data.empty()) {
        throw std::runtime_error("Determinant can only be calculated for a non empty square matrix");
    }

    int n = data.size();
    if (n == 1) {
        return data[0][0];
    }
    double det = 0.0;
    for (int i = 0; i < n; ++i) {
        std::vector<std::vector<int>> minor(n - 1, std::vector<int>(n - 1, 0));
        for (int j = 1; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (k < i) {
                    minor[j - 1][k] = data[j][k];
                } else if (k > i) {
                    minor[j - 1][k - 1] = data[j][k];
                }
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * data[0][i] * Matrix(minor).determinant();
    }
    return det;
}

/**
 * @brief Computes the transpose of the matrix.
 * 
 * @return The transposed matrix.
 */
Matrix Matrix::transpose() {
    size_t rows = data.size();
    size_t cols = data[0].size();

    // Initialize the transpose matrix with appropriate dimensions:
    std::vector<std::vector<int>> transposedData(cols, std::vector<int>(rows));

    // Transpose the matrix
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            transposedData[j][i] = data[i][j];
        }
    }
    return Matrix(transposedData);
}

/**
 * @brief Computes the rank of the matrix.
 * 
 * @return The rank of the matrix.
 */
int Matrix::rank() {
    std::vector<std::vector<int>> copyData = data; // Make a copy of the original matrix

    size_t rows = copyData.size();
    size_t cols = copyData[0].size();

    size_t rank = 0;
    size_t col = 0;

    // Row reduction to bring it into row echelon form:
    for (size_t row = 0; row < rows; ++row) {
        // Find the first non-zero element in the current column:
        while (col < cols && copyData[row][col] == 0) {
            ++col;
        }
        // If we find a non-zero element, increment the rank and continue the row reduction
        if (col < cols) {
            rank++;
            // Eliminate the non-zero elements below the current pivot
            for (size_t i = row + 1; i < rows; ++i) {
                int factor = copyData[i][col] / copyData[row][col];
                for (size_t j = col; j < cols; ++j) {
                    copyData[i][j] -= factor * copyData[row][j];
                }
            }
        }
    }
    return rank;
}


/**
 * @brief Displays the contents of the matrix.
 */
void Matrix::display() {
    for (const auto& row : data) {
        for (int value : row) {
            std::cout<< value <<"\t";
        }
        std::cout<<std::endl;
    }
}

/**
 * @brief Computes the inverse of a matrix using the Gauss-Jordan elimination.
 * 
 * @return The inverse of a matrix.
 * 
 * @throw std::runtime_error If matrix is not square or matrix is singular.
*/
Matrix Matrix::inverse() {
    //Check if matrix is square:
    if (data.size() != data[0].size()) {
        throw std::runtime_error("Inverse can only be computed for a square matrix.");
    }


    //Helper function to perform row operations and make elements below diagonal zero.
    auto makeElementsBelowZero = [this](size_t row, size_t col) {
        for (size_t i = 0; i < data.size(); ++i) {
            if (i != row) {
                double factor = data[i][col];
                for (size_t j = 0; j < data[i].size(); ++j) {
                    data[i][j] -= factor * data[row][j];
                }
            }
        }
    };
    // Helper function to perform row operations and make diagonal elements one.
    auto makeDiagonalOne = [this](size_t row, size_t col) {
        double pivot = data[row][col];
        for (size_t i = col; i < data[row].size(); ++i) {
            data[row][i] /= pivot;
        }
    };
    size_t n = data.size();
    // Create an augmented matrix [A|I]
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(2*n, 0));
    for (size_t i = 0; i < n; ++i) {
        augmentedMatrix[i][n + i] = 1; //Identity part
        for (size_t j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = data[i][j]; //Original matrix part
        }
    }
    // Apply Gauss-Jordan elimination:
    for (size_t i = 0; i < n; ++i) {
        // Find the pivot for the current row:
        size_t pivotRow = i;
        while (pivotRow < n && augmentedMatrix[pivotRow][i] == 0) {
            ++pivotRow;
        }
        if (pivotRow == n) {
            throw std::runtime_error("This is a singular matrix, hence inverse cannot be computed.");
        }
        // Swap the pivot row with the current row:
        if (pivotRow != i) {
            std::swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);
        }
        // Make diagonal element one:
        makeDiagonalOne(i, i);

        //Make elements below diagonal zero:
        makeElementsBelowZero(i, i);
    }
    std::vector<std::vector<double>> inverseMatrix(n, std::vector<double>(n, 0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inverseMatrix[i][j] = augmentedMatrix[i][j + n];
        }
    }
    return Matrix(inverseMatrix);
}