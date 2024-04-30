# Matrix Calculator

This project consists of two main components: the Matrix class and the MatrixCalculator class. The Matrix class provides functionality for creating and manipulating matrices, while the MatrixCalculator class offers various matrix operations such as addition, subtraction, element-wise product, and dot product.

## Matrix Class

### Constructor

- `Matrix(const std::vector<std::vector<int>>& input)`: Constructs a matrix object from a given 2D vector of integers.

### Public Methods

- `std::vector<std::vector<int>> getData() const`: Returns the data stored in the matrix.
- `int getRows() const`: Returns the number of rows in the matrix.
- `int getCols() const`: Returns the number of columns in the matrix.
- `double determinant()`: Computes the determinant of the matrix.
- `Matrix transpose()`: Computes the transpose of the matrix.
- `int rank()`: Computes the rank of the matrix.
- `void display()`: Displays the contents of the matrix.
- `Matrix inverse()`: Computes the iverse of the matrix.

## MatrixCalculator Class

### Public Methods

- `Matrix addMatrices(const Matrix& A, const Matrix& B)`: Adds two matrices element-wise.
- `Matrix subtractMatrices(const Matrix& A, const Matrix& B)`: Subtracts one matrix from another element-wise.
- `Matrix elementWiseProduct(const Matrix& A, const Matrix& B)`: Computes the element-wise product of two matrices.
- `Matrix dotProduct(const Matrix& A, const Matrix& B)`: Computes the dot product of two matrices.

## Usage

### Matrix Class

#### Creating a Matrix

```cpp
#include "Matrix.h"

// Create a matrix with provided data
std::vector<std::vector<int>> data = {{1, 2}, {3, 4}};
Matrix A(data);
```

#### Accessing Matrix Properties

```cpp
#include "Matrix.h"

Matrix A(data);

// Get the number of rows and columns
int rows = A.getRows();
int cols = A.getCols();
```

#### Performing Matrix Operations

```cpp
#include "Matrix.h"

Matrix A(data);

// Compute determinant
double det = A.determinant();

// Transpose the matrix
Matrix transposed = A.transpose();

// Calculate the rank of the matrix
int matrixRank = A.rank();
```

#### Displaying Matrix Contents

```cpp
#include "Matrix.h"

Matrix A(data);

// Display matrix contents
A.display();
```

### MatrixCalculator Class

#### Creating a MatrixCalculator Object

```cpp
#include "MatrixCalculator.h"

MatrixCalculator calculator;
```

#### Performing Matrix Operations

```cpp
#include "MatrixCalculator.h"

// Add two matrices
Matrix sum = calculator.addMatrices(A, B);

// Subtract one matrix from another
Matrix difference = calculator.subtractMatrices(A, B);

// Compute element-wise product
Matrix elementWiseProduct = calculator.elementWiseProduct(A, B);

// Compute dot product
Matrix dotProduct = calculator.dotProduct(A, B);
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
