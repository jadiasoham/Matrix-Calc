#include "Matrix.h"
#include <iostream>

/**
 * @brief Constructs a matrix from the given input data.
 * 
 * @param input The input data representing the matrix.
 */
Matrix::Matrix(const std::vector<std::vector<int>>& input) : data(input) {}

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