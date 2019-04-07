#include "logger.h"
#include "matrix.h"

#include <cctype>
#include <iostream>
#include <thread>

using namespace matrix;

Matrix CreateVectorColumn(std::size_t number_of_rows);
Matrix CreateVectorRow(std::size_t number_of_columns);

Matrix SolveEquation(const std::size_t n);

void GetInputData(
  Matrix *A,
  Matrix *A1,
  Matrix *b1,
  Matrix *c1,
  Matrix *B2,
  Matrix *C2
);
void InputMatrix(Matrix *matrix, const char *matrix_name);
Matrix Get_C2(std::size_t n);
Matrix Get_b(std::size_t n);

int main() {
  std::size_t n = 0;

  std::cout << "Please enter n:" << std::endl;
  std::cin >> n;

  Logger::getInstance().log(SolveEquation(n).to_string());

  system("pause");
  return 0;
}

Matrix SolveEquation(const std::size_t n) {
  /*****************************************************************************************
   * First level
   ******************************************************************************************/

  Matrix B2(n, n);
  Matrix A2(n, n);
  Matrix b1 = CreateVectorColumn(n);
  Matrix c1 = CreateVectorColumn(n);
  Matrix A1(n, n);
  Matrix A(n, n);
  Matrix C2{};
  Matrix b{};

  // Cij = (1 / (i + 2j))
  auto C2_thread = std::thread([&] { C2 = Get_C2(n); });
  // bi = (1 / (i*i + 1)) || bi = (1 / i)
  auto b_thread = std::thread([&] { b = Get_b(n); });

  GetInputData(&A, &A1, &b1, &c1, &B2, &C2);

  C2_thread.join();
  b_thread.join();

  /*****************************************************************************************
   * Second level
   ******************************************************************************************/

  Matrix B2_minus_C2{};
  Matrix b1_plus_c1{};

  auto B2_minus_C2_thread = std::thread([&] { B2_minus_C2 = B2 - C2; });
  auto b1_plus_c1_thread = std::thread([&] { b1_plus_c1 = b1 + c1; });

  B2_minus_C2_thread.join();
  b1_plus_c1_thread.join();

  /*****************************************************************************************
   * Third level
   ******************************************************************************************/

  Matrix y1{};
  Matrix y2{};
  Matrix Y3{};

  auto y1_thread = std::thread([&] { y1 = A * b; });
  auto y2_thread = std::thread([&] { y2 = A1 * b1_plus_c1; });
  auto Y3_thread = std::thread([&] { Y3 = A2 * B2_minus_C2; });

  y1_thread.join();
  y2_thread.join();
  Y3_thread.join();

  /*****************************************************************************************
   * Fourth level
   ******************************************************************************************/

  Matrix Y3_squared{};
  Matrix y1_plus_y2{};

  auto Y3_squared_thread = std::thread([&] { Y3_squared = Y3 * Y3; });
  auto y1_plus_y2_thread = std::thread([&] { y1_plus_y2 = y1 + y2; });

  Y3_squared_thread.join();
  y1_plus_y2_thread.join();

  /*****************************************************************************************
   * Fifth level
   ******************************************************************************************/

  Matrix Y3_squared_mupltiplied_by_y2{};
  Matrix Y3_multiplied_by_sum_of_y1_and_y2{};

  auto Y3_squared_mupltiplied_by_y2_thread =
    std::thread([&] { Y3_squared_mupltiplied_by_y2 = Y3_squared * y2; });
  auto Y3_multiplied_by_sum_of_y1_and_y2_thread =
      std::thread([&] { Y3_multiplied_by_sum_of_y1_and_y2 = Y3 * y1_plus_y2; });

  Y3_squared_mupltiplied_by_y2_thread.join();
  Y3_multiplied_by_sum_of_y1_and_y2_thread.join();

  /*****************************************************************************************
   * Sixth level
   ******************************************************************************************/

  return Y3_squared_mupltiplied_by_y2 + Y3_multiplied_by_sum_of_y1_and_y2;
}

void GetInputData(
  Matrix *A,
  Matrix *A1,
  Matrix *b1,
  Matrix *c1,
  Matrix *B2,            
  Matrix *C2
) {
  char is_random_data = '\0';
  std::cout << "Random data (Y/N): ";
  std::cin >> is_random_data;

  if (std::toupper(is_random_data) == 'Y') {
    A->fill_random(0, 10);
    A1->fill_random(0, 10);
    b1->fill_random(0, 10);
    c1->fill_random(0, 10);
    B2->fill_random(0, 10);
    C2->fill_random(0, 10);
  } else {
    InputMatrix(A, "A");
    InputMatrix(A1, "A1");
    InputMatrix(b1, "b1");
    InputMatrix(c1, "c1");
    InputMatrix(B2, "B2");
    InputMatrix(C2, "C2");
  }
}

void InputMatrix(Matrix *matrix, const char *matrix_name) {
  std::cout << "Please enter matrix " << matrix_name << std::endl;

  for (int i = 0; i < matrix->rows(); i++) {
    for (int j = 0; j < matrix->cols(); j++) {
      std::cout << "Element (" << i + 1 << ", " << j + 1 << ") = ";
      std::cin >> matrix->at(i, j);
    }
  }
}

Matrix Get_C2(std::size_t n) {
  Matrix C2(n, n);

  for (int i = 1; i <= C2.rows(); i++) {
    for (int j = 1; j <= C2.cols(); j++) {
      C2.at(i - 1, j - 1) = 1.0 / (i + 2 * j);
    }
  }

  return C2;
}

Matrix Get_b(std::size_t n) {
  Matrix b = CreateVectorColumn(n);

  for (int i = 1; i <= b.rows(); i++) {
    if ((i - 1) & 0x01) {  // odd
      b.at(i - 1, 0) = 1.0 / i;
    } else {  // even
      b.at(i - 1, 0) = 1.0 / (i * i + 2);
    }
  }

  return b;
}

Matrix CreateVectorColumn(std::size_t number_of_rows) {
  return Matrix(number_of_rows, 1);
}

Matrix CreateVectorRow(std::size_t number_of_columns) {
  return Matrix(1, number_of_columns);
}
