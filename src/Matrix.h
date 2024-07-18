#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cmath>

class Matrix {
private:
    std::vector<std::vector<double>> data;

public:
    Matrix(const std::vector<std::vector<int>>& input); //Constructor for int data
    Matrix(const std::vector<std::vector<double>>& input); //Constructor for double data
    std::vector<std::vector<double>> getData() const;
    int getRows() const;
    int getCols() const;
    double determinant();
    Matrix transpose();
    int rank();
    double trace();
    void display();
    Matrix inverse();
    Matrix replaceRows(int n, std::vector<double>& data);
    Matrix replaceCols(int n, std::vector<double>& data);
};

#endif
