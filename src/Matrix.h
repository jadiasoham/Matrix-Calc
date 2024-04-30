#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

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
    void display();
    Matrix inverse();
};

#endif
