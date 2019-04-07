#include <iostream>
#include <vector>
#include <random>
#include <cstddef>
#include <chrono>

using Matrix = std::vector<std::vector<double>>;

constexpr std::size_t kAlignment = 8;

Matrix CreateMatrix(std::size_t n);
Matrix CreateMatrixA(std::size_t n);
Matrix CreateMatrixB(std::size_t n);

Matrix operator*(const Matrix &lhs, const Matrix &rhs);
void IterateRows(
  const Matrix &lhs, const Matrix &rhs,
  Matrix &result,
  std::size_t row = 0
);
void IterateColumns(
  const Matrix &lhs, const Matrix &rhs,
  Matrix &result, std::size_t row,
  std::size_t col = 0
);
double CalculateElement(
  const Matrix &lhs, const Matrix &rhs,
  std::size_t row, std::size_t col,
  std::size_t k = 0
);

std::ostream& operator<<(std::ostream &stream, const Matrix &matrix);

int main() {
  std::size_t n{};
  std::cout << "Enter n: ";
  std::cin >> n;

  auto A = CreateMatrixA(n);
  auto B = CreateMatrixB(n);

  auto C = A * B;

  std::cout << A << B << C;

  return 0;
}

Matrix CreateMatrix(std::size_t n) {
  Matrix matrix(n);
  for (auto &row : matrix) {
    row.resize(n);
  }
  return matrix;
}

Matrix CreateMatrixA(std::size_t n) {
  auto A = CreateMatrix(n);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      A[i][i - j] = n - j;
    }
  }

  return A;
}

Matrix CreateMatrixB(std::size_t n) {
  auto B = CreateMatrix(n);

  // fill upper half
  for (std::size_t i = 0; i <= n/2; ++i) {
    for (std::size_t j = n - 1; j >= (n - 1 - i); --j) {
      B[i][j] = i * j + 1;
    }
  }

  // fill lower half
  for (std::size_t i = n/2 + 1; i < n; ++i) {
    for (std::size_t j = i; j < n; ++j) {
      B[i][j] = i * j + 1;
    }
  }

  return B;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs) {
  auto n = lhs.size();
  auto result = CreateMatrix(n);

  IterateRows(lhs, rhs, result);

  return result;
}

void IterateRows(
  const Matrix &lhs, const Matrix &rhs,
  Matrix &result,
  std::size_t row
) {
  if (row >= lhs.size()) {
    return;
  }
  IterateColumns(lhs, rhs, result, row);
  IterateRows(lhs, rhs, result, row + 1);
}

void IterateColumns(
  const Matrix &lhs, const Matrix &rhs,
  Matrix &result, std::size_t row,
  std::size_t col
) {
  if (col >= lhs.size()) {
    return;
  }

  result[row][col] = CalculateElement(lhs, rhs, row, col);
  IterateColumns(lhs, rhs, result, row, col + 1);
}

double CalculateElement(
  const Matrix &lhs, const Matrix &rhs,
  std::size_t row, std::size_t col,
  std::size_t k
) {
  if (k >= lhs.size()) {
    return 0.0;
  }
  return (lhs[row][k] * rhs[k][col]
      + CalculateElement(lhs, rhs, row, col, k + 1));
}

std::ostream& operator<<(std::ostream &stream, const Matrix &matrix) {
  stream << "Matrix:\n";
  for (const auto &row : matrix) {
    for (const auto &element : row) {
      std::cout.width(kAlignment);
      stream << element << ' ';
    }
    stream << std::endl;
  }
  return stream;
}
