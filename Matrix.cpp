#include "Matrix.h"

using namespace std;

Matrix::Matrix(const vector<vector<int>>& input) : data(input) {}

vector<vector<int>> Matrix::getData() const {
    return data;
}

int Matrix::getRows() const {
    return data.size();
}

int Matrix::getCols() const {
    if (data.empty()) {
        return 0;
    }
    return data[0].size();
}