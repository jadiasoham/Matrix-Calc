#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cmath>

class Matrix {
private:
    std::vector<std::vector<int>> data;

public:
    Matrix(const std::vector<std::vector<int>>& input);
    Matrix(const std::vector<std::vector<double>>& input);
    std::vector<std::vector<int>> getData() const;
    int getRows() const;
    int getCols() const;
    double determinant();
    Matrix transpose();
    int rank();
    double trace();
    void display();
    Matrix inverse();
};

#endif
