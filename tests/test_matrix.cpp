/*
This test file is for Matrix.cpp
It uses the google test framework
*/

#include <gtest/gtest.h>
#include "Matrix.h"

// Test case for the matrix class
TEST(MatrixTest, ConstructorAndGetters) {
  // Test constructor:
  std::vector<std::vector<int>> input = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  Matrix matrix(input);

  // Test getters:
  ASSERT_EQ(matrix.getRows(), 3);
  ASSERT_EQ(matrix.getCols(), 3);
  ASSERT_EQ(matrix.getData(), input);
}

// More test cases may be added if need arise.

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
