#include "matrix.h"

#include <iostream>

matrix::Matrix CreateVectorColumn(std::size_t number_of_rows);
matrix::Matrix CreateVectorRow(std::size_t number_of_columns);

int main() {
  using namespace matrix;
  std::cout << CreateVectorColumn(5);
  std::cout << CreateVectorRow(5);

  system("pause");
  return 0;
}

matrix::Matrix CreateVectorColumn(std::size_t number_of_rows) {
  return matrix::Matrix(number_of_rows, 1);
}

matrix::Matrix CreateVectorRow(std::size_t number_of_columns) {
  return matrix::Matrix(1, number_of_columns);
}
