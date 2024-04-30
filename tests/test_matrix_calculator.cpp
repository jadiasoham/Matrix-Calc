#include <gtest/gtest.h>
#include "MatrixCalculator.h"

// Test case for MatrixCalculator class:
TEST(MatrixCalculatorTest, Addition) {
  std::vector<std::vector<int>> input1 = {{1, 2}, {3, 4}};
  std::vector<std::vector<int>> input2 = {{5, 6}, {7, 8}};

  Matrix A(input1);
  Matrix B(input2);

  MatrixCalculator calculator;
  Matrix result = calculator.addMatrices(A, B);

  std::vector<std::vector<int>> expected = {{6, 8}, {10, 12}};
  ASSERT_EQ(result.getData(), expected);
}

// Additional test cases may be added to add other functionalities.

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
