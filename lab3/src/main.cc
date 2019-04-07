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
double CalculateElement(
  const Matrix &lhs, const Matrix &rhs,
  std::size_t row, std::size_t col,
  std::size_t current = 0
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

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      result[i][j] = CalculateElement(lhs, rhs, i, j);
    }
  }

  return result;
}

double CalculateElement(
  const Matrix &lhs, const Matrix &rhs,
  std::size_t row, std::size_t col,
  std::size_t current
) {
  return (current < lhs.size()
      ? lhs[row][current] * rhs[current][col]
          + CalculateElement(lhs, rhs, row, col, current + 1)
      : 0.0);
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
