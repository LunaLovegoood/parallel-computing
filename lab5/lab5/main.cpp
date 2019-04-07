#define CL_USE_DEPRECATED_OPENCL_2_0_APIS

#include "CL/cl.hpp"
#include "matrix.h"

#include <iostream>
#include <assert.h>

#define CALCULATION_TYPE float
#define DEFAULT_RANDOM_RANGE 0, 10

template <typename T>
Matrix<T> CreateVectorColumn(std::size_t size);
template <typename T>
Matrix<T> CreateRandomVectorColumn(std::size_t size);

template <typename T>
Matrix<T> Calculate_b(std::size_t size);
template <typename T>
Matrix<T> Calculate_C2(std::size_t size);

template <typename T>
Matrix<T> Calculate_y1(std::size_t size);
template <typename T>
Matrix<T> Calculate_y2(std::size_t size);
template <typename T>
Matrix<T> Calculate_Y3(std::size_t size);

int main() {
  assert(CheckOpenCLSupport(CL_DEVICE_TYPE_GPU));
  std::cout << "OpenCl is supported" << std::endl;
  
  std::size_t n{};
  std::cout << "Enter n: ";
  std::cin >> n;

  auto y1 = Calculate_y1<CALCULATION_TYPE>(n);
  auto y2 = Calculate_y2<CALCULATION_TYPE>(n);
  auto Y3 = Calculate_Y3<CALCULATION_TYPE>(n);

  auto result = (Y3 * Y3) * y2 + Y3 * (y1 + y2);

  std::cout << "Result:\n";
  std::cout << result << std::endl;

  system("pause");
  return 0;
}

template <typename T>
Matrix<T> CreateVectorColumn(std::size_t size) {
  return Matrix<T>::Empty(size, 1);
}

template <typename T>
Matrix<T> CreateRandomVectorColumn(std::size_t size) {
  return Matrix<T>::Random(size, 1, DEFAULT_RANDOM_RANGE);
}

template <typename T>
Matrix<T> Calculate_b(std::size_t size) {
  auto b = CreateVectorColumn<T>(size);

  for (std::size_t i = 1; i <= size; i++) {
    if (i % 2 == 0) {
      b[i - 1][1] = 1 / T(i*i + 2);
    } else {
      b[i - 1][1] = 1 / T(i);
    }
  }

  return b;
}

template <typename T>
Matrix<T> Calculate_C2(std::size_t size) {
  auto C2 = Matrix<T>::Random(size, size, DEFAULT_RANDOM_RANGE);

  for (std::size_t i = 1; i <= size; i++) {
    for (std::size_t j = 1; j <= size; j++) {
      C2[i - 1][j - 1] = 1 / T(i + 2*j);
    }
  }

  return C2;
}

template <typename T>
Matrix<T> Calculate_y1(std::size_t n) {
  auto A = Matrix<T>::Random(n, n, DEFAULT_RANDOM_RANGE);
  auto b = Calculate_b<T>(n);

  return A * b;
}

template <typename T>
Matrix<T> Calculate_y2(std::size_t n) {
  auto A1 = Matrix<T>::Random(n, n, DEFAULT_RANDOM_RANGE);
  auto b1 = CreateRandomVectorColumn<T>(n);
  auto c1 = CreateRandomVectorColumn<T>(n);

  return A1 * (b1 + c1);
}

template <typename T>
Matrix<T> Calculate_Y3(std::size_t n) {
  auto A2 = Matrix<T>::Random(n, n, DEFAULT_RANDOM_RANGE);
  auto B2 = Matrix<T>::Random(n, n, DEFAULT_RANDOM_RANGE);
  auto C2 = Calculate_C2<T>(n);

  return A2 * (B2 - C2);
}
