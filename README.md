# Matrix Calculator

This is a simple C++ program that provides basic matrix operations such as addition, subtraction, multiplication, and determinant calculation.

## Features

- Addition of two matrices
- Subtraction of two matrices
- Multiplication of two matrices
- Determinant calculation of a square matrix

## Usage

### Compilation

To compile the program, ensure you have a C++ compiler installed on your system. Then, navigate to the directory containing the source code and execute the following command:

```bash
g++ main.cpp -o matrix_calculator
```

### Running the Program

After compiling, run the executable file:

```bash
./matrix_calculator
```

The program will execute without any output as the `main` function is currently empty. To use the matrix calculator functions, you can include this header file in your own C++ program and instantiate the `MatrixCalculation` class.

### Example Usage

Here's an example of how to use the matrix calculator functions:

```cpp
#include "MatrixCalculation.h"
#include <iostream>

int main() {
    MatrixCalculation matrixCalc;

    // Example matrices
    std::vector<std::vector<int>> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::vector<std::vector<int>> B = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

    // Addition
    std::vector<std::vector<int>> sum = matrixCalc.addMatrices(A, B);
    std::cout << "Sum of matrices A and B:" << std::endl;
    printMatrix(sum);

    // Subtraction
    std::vector<std::vector<int>> difference = matrixCalc.subtractMatrices(A, B);
    std::cout << "Difference of matrices A and B:" << std::endl;
    printMatrix(difference);

    // Multiplication
    std::vector<std::vector<int>> product = matrixCalc.multiplyMatrices(A, B);
    std::cout << "Product of matrices A and B:" << std::endl;
    printMatrix(product);

    // Determinant
    int detA = matrixCalc.determinant(A);
    std::cout << "Determinant of matrix A: " << detA << std::endl;

    return 0;
}
```

## Contributing

Feel free to contribute to this project by forking the repository, making changes, and submitting a pull request. Suggestions, improvements, and bug fixes are welcome!

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
