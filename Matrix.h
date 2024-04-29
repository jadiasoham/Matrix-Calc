#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

using namespace std;

class Matrix {
private:
    vector<vector<int>> data;

public:
    Matrix(const vector<vector<int>>& input);
    vector<vector<int>> getData() const;
    int getRows() const;
    int getCols() const;
};

#endif