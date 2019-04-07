#ifndef LAB5_MATRIX_H_
#define LAB5_MATRIX_H_

#include "opencl_wrapper_functions.h"

#include <cstddef>
#include <string>
#include <memory>
#include <type_traits>
#include <utility>
#include <iostream>
#include <assert.h>
#include <random>
#include <chrono>

namespace {

// Default build options for program build
#if _WIN64 
  constexpr char kDefaultBuildOptions[] = "-cl-std=CL1.2 -D __x__64__";
#else
  constexpr char kDefaultBuildOptions[] = "-cl-std=CL1.2";
#endif

// Kernel paths
constexpr char kNegativeKernel[] = "kernels\\matrix\\negative.cl";
constexpr char kAdditionKernel[] = "kernels\\matrix\\addition.cl";
constexpr char kSubtractionKernel[] = "kernels\\matrix\\subtraction.cl";
constexpr char kMultiplicationKernel[] = "kernels\\matrix\\multiplication.cl";

// Matrix-with-Matrix kernel names without type suffix
constexpr char kAdditionKernelName[] = "matrix_add";
constexpr char kSubtractionKernelName[] = "matrix_subtract";
constexpr char kMultiplicationKernelName[] = "matrix_multiply";
constexpr char kNegativeKernelName[] = "negative";

// Matrix-with-scalar kernel names without type suffix
constexpr char kScalarAdditionKernelName[] = "matrix_scalar_add";
constexpr char kScalarMultiplicationKernelName[] = "matrix_scalar_multiply";

template <typename T>
struct dependent_false : std::false_type {};

// Used to get correct type suffix for the kernel name
template <typename T>
inline std::string GetTypeSuffix();

struct {
  std::unique_ptr<cl::Program> addition{nullptr};
  std::unique_ptr<cl::Program> subtraction{nullptr};
  std::unique_ptr<cl::Program> multiplication{nullptr};
  std::unique_ptr<cl::Program> negative{nullptr};
} kMatrixPrograms;

// Used to initialize program if it wasn't initialized before
#define LAZY_PROGRAM_INIT(program_name, kernel_path)                        \
  if (!kMatrixPrograms.program_name) {                                      \
    kMatrixPrograms.program_name.reset(                                     \
        new cl::Program(CreateProgram(kernel_path, kDefaultBuildOptions))); \
  }

#define OPERATION_INIT(rows, cols)                                  \
auto context = program.getInfo<CL_PROGRAM_CONTEXT>();               \
auto device = context.getInfo<CL_CONTEXT_DEVICES>().front();        \
                                                                    \
Matrix<T> result = Matrix<T>::Empty(rows, cols);                    \
                                                                    \
cl::Buffer lhsBuffer(context,                                       \
                     CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS |     \
                         CL_MEM_USE_HOST_PTR,                       \
                     lhs.size_, lhs.matrix_);                       \
cl::Buffer rhsBuffer(context,                                       \
                     CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS |     \
                         CL_MEM_USE_HOST_PTR,                       \
                     rhs.size_, rhs.matrix_);                       \
cl::Buffer resultBuffer(context,                                    \
                        CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY | \
                            CL_MEM_USE_HOST_PTR,                    \
                        result.size_, result.matrix_);              \
                                                                    \
cl::Kernel kernel(program, kernel_name.data());                     \
kernel.setArg(0, lhsBuffer);                                        \
kernel.setArg(1, rhsBuffer);                                        \
kernel.setArg(2, resultBuffer)

} // namespace

template <typename T>
class Matrix {
  static_assert(std::is_floating_point<T>::value,
    "Matrix element's type should be floating point type");

  class MatrixRow {
   public:
    MatrixRow(T *row_start) : row_start_(row_start) {}

    T& operator[](std::size_t j) { return *(row_start_ + j); }

   private:
    T *row_start_{nullptr};
  };

 public:
  Matrix() = default;
  Matrix(const Matrix &matrix);
  Matrix(Matrix &&matrix);
  ~Matrix();

  static Matrix Zeros(std::size_t order);
  
  static Matrix Empty(std::size_t rows, std::size_t cols);

  static Matrix From1D(const T *matrix, std::size_t rows, std::size_t cols);
  static Matrix From1D(const T *matrix, std::size_t order);

  static Matrix Random(std::size_t rows, std::size_t cols,
                       T lower_bound, T upper_bound);

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }

  // Accessor
  // Does NOT perform boundary checks
  MatrixRow operator[](std::size_t i) { return MatrixRow(matrix_ + i*cols_); }

  // Fills matrix with given value
  void Fill(T value);

  Matrix& operator=(Matrix matrix);

  // Swaps contents of two matrices
  template <typename T>
  friend void swap(Matrix<T> &first, Matrix<T> &second);
  
  Matrix<T> operator-() const;
  template <typename T>
  friend Matrix<T> operator+(const Matrix<T> &lhs, const Matrix<T> &rhs);
  template <typename T>
  friend Matrix<T> operator-(const Matrix<T> &lhs, const Matrix<T> &rhs);
  template <typename T>
  friend Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs);

  template <typename T>
  friend Matrix<T> operator+(const T scalar, const Matrix<T> &matrix);
  template <typename T>
  friend Matrix<T> operator+(const Matrix<T> &matrix, const T scalar);
  template <typename T>
  friend Matrix<T> operator-(const T scalar, const Matrix<T> &matrix);
  template <typename T>
  friend Matrix<T> operator-(const Matrix<T> &matrix, const T scalar);
  template <typename T>
  friend Matrix<T> operator*(const T scalar, const Matrix<T> &matrix);
  template <typename T>
  friend Matrix<T> operator*(const Matrix<T> &matrix, const T scalar);

 private:
  T *matrix_{ nullptr };
  std::size_t rows_{};
  std::size_t cols_{};
  std::size_t size_{}; // size in bytes

  Matrix(std::size_t rows, std::size_t cols, const T *matrix = nullptr);

  static Matrix ExecuteBinaryMatrixOperation(
      const Matrix<T> &lhs,
      const Matrix<T> &rhs,
      cl::Program program,
      const std::string &kernel_name
  );
  static Matrix ExecuteBinaryScalarMatrixOperation(
    const T scalar,
    const Matrix<T> matrix,
    cl::Program program,
    const std::string &kernel_name
  );
};

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &matrix) {
  if (matrix.matrix_) {
    matrix_ = new T[matrix.size_];
    memcpy(matrix_, matrix.matrix_, matrix.size_);
  }

  rows_ = matrix.rows_;
  cols_ = matrix.cols_;
  size_ = matrix.size_;
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> &&matrix) : Matrix() {
  swap(*this, matrix);
}

template <typename T>
Matrix<T>::Matrix(std::size_t rows, std::size_t cols, const T *matrix) {
  rows_ = rows;
  cols_ = cols;
  size_ = sizeof(T) * rows * cols;

  if (matrix) {
    matrix_ = new T[rows * cols];
    memcpy(matrix_, matrix, size_);
  } else {
    matrix_ = new T[rows * cols]{};
  }
}

template <typename T>
Matrix<T>::~Matrix() {
  delete[] matrix_;
}

template <typename T>
Matrix<T> Matrix<T>::Zeros(std::size_t order) {
  return Matrix<T>(order, order);
}

template <typename T>
Matrix<T> Matrix<T>::Empty(std::size_t rows, std::size_t cols) {
  return Matrix<T>(rows, cols);
}

template <typename T>
Matrix<T> Matrix<T>::From1D(const T *matrix, std::size_t rows,
                            std::size_t cols) {
  assert(matrix);
  return Matrix<T>(rows, cols, matrix);
}

template <typename T>
Matrix<T> Matrix<T>::From1D(const T *matrix, std::size_t order) {
  assert(matrix);
  return Matrix<T>(order, order, matrix);
}

template <typename T>
Matrix<T> Matrix<T>::Random(std::size_t rows, std::size_t cols,
                            T lower_bound, T upper_bound) {
  auto matrix = Matrix<T>::Empty(rows, cols);

  unsigned seed = static_cast<unsigned>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  std::mt19937 eng(seed);
  std::uniform_real_distribution<T> distribution(lower_bound, upper_bound);

  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      matrix.matrix_[i*cols + j] = distribution(eng);
    }
  }

  return matrix;
}

template <typename T>
void Matrix<T>::Fill(T value) {
  if (!matrix_) {
    return;
  }

  auto device = GetFirstGPUDevice();
  cl::Context context(device);

  cl::Buffer buffer(
      context,
      CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR,
      size_,
      matrix_
  );

  cl::CommandQueue cmd_queue(context, device);
  cmd_queue.enqueueFillBuffer(buffer, value, 0, size_);
  cmd_queue.enqueueReadBuffer(buffer, CL_TRUE, 0, size_, matrix_);
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> rhs) {
  swap(*this, rhs);
  return *this;
}

template <typename T>
void swap(Matrix<T> &first, Matrix<T> &second) {
  using std::swap;

  swap(first.matrix_, second.matrix_);
  swap(first.rows_, second.rows_);
  swap(first.cols_, second.cols_);
  swap(first.size_, second.size_);
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const {
  if (!matrix_) {
    return Matrix<T>();
  }

  LAZY_PROGRAM_INIT(negative, kNegativeKernel);

  auto program = *kMatrixPrograms.negative;
  auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
  auto device = context.getInfo<CL_CONTEXT_DEVICES>().front();

  Matrix<T> result = Matrix<T>::Empty(rows_, cols_);

  cl::Buffer inBuffer(
      context, CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_USE_HOST_PTR,
      size_, matrix_);
  cl::Buffer resultBuffer(
      context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR,
      result.size_, result.matrix_);

  std::string kernel_name = std::string(kNegativeKernelName) + GetTypeSuffix<T>();
  cl::Kernel kernel(program, kernel_name.data());
  kernel.setArg(0, inBuffer);      // in
  kernel.setArg(1, resultBuffer);  // result

  cl::CommandQueue cmd_queue(context, device);
  cmd_queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                 cl::NDRange(result.rows_ * result.cols_));
  cl::finish();

  return result;
}

template <typename T>
Matrix<T> operator+(const Matrix<T> &lhs, const Matrix<T> &rhs) {
  LAZY_PROGRAM_INIT(addition, kAdditionKernel);
  auto full_kernel_name =
      std::string(kAdditionKernelName) + GetTypeSuffix<T>();

  return Matrix<T>::ExecuteBinaryMatrixOperation(
      lhs, rhs, *kMatrixPrograms.addition, full_kernel_name
  );
}

template <typename T>
Matrix<T> operator-(const Matrix<T> &lhs, const Matrix<T> &rhs) {
  LAZY_PROGRAM_INIT(subtraction, kSubtractionKernel);
  auto full_kernel_name =
      std::string(kSubtractionKernelName) + GetTypeSuffix<T>();

  return Matrix<T>::ExecuteBinaryMatrixOperation(
      lhs, rhs, *kMatrixPrograms.subtraction, full_kernel_name
  );
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs) {
  if (!lhs.matrix_ || !rhs.matrix_ ||  // check if args are initialized
      lhs.cols_ != rhs.rows_           // and if dimensions are correct
  ) {
    return Matrix<T>();
  }
  
  LAZY_PROGRAM_INIT(multiplication, kMultiplicationKernel);
  auto kernel_name =
      std::string(kMultiplicationKernelName) + GetTypeSuffix<T>();

  auto program = *kMatrixPrograms.multiplication;
  
  OPERATION_INIT(lhs.rows_, rhs.cols_);
  kernel.setArg(3, lhs.cols_);
  kernel.setArg(4, rhs.cols_);

  cl::CommandQueue cmd_queue(context, device);
  cmd_queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                 cl::NDRange(result.rows_, result.cols_));
  cl::finish();

  return result;
}

template <typename T>
Matrix<T> operator+(const T scalar, const Matrix<T> &matrix) {
  LAZY_PROGRAM_INIT(addition, kAdditionKernel);
  auto full_kernel_name =
      std::string(kScalarAdditionKernelName) + GetTypeSuffix<T>();

  return Matrix<T>::ExecuteBinaryScalarMatrixOperation(
      scalar, matrix, *kMatrixPrograms.addition, full_kernel_name);
}

template <typename T>
Matrix<T> operator+(const Matrix<T> &matrix, const T scalar) {
  return (scalar + matrix);
}

template <typename T>
Matrix<T> operator-(const T scalar, const Matrix<T> &matrix) {
  return scalar + (-matrix);
}

template <typename T>
Matrix<T> operator-(const Matrix<T> &matrix, const T scalar) {
  return matrix + (-scalar);
}

template <typename T>
Matrix<T> operator*(const T scalar, const Matrix<T> &matrix) {
  LAZY_PROGRAM_INIT(addition, kMultiplicationKernel);
  auto full_kernel_name =
      std::string(kScalarMultiplicationKernelName) + GetTypeSuffix<T>();

  return Matrix<T>::ExecuteBinaryScalarMatrixOperation(
      scalar, matrix, *kMatrixPrograms.multiplication, full_kernel_name);
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &matrix, const T scalar) {
  return scalar * matrix;
}

template <typename T>
Matrix<T> Matrix<T>::ExecuteBinaryMatrixOperation(
    const Matrix<T> &lhs,
    const Matrix<T> &rhs,
    cl::Program program,
    const std::string &kernel_name
) {
  if (!lhs.matrix_ || !rhs.matrix_ ||  // check if args are initialized
      lhs.rows_ != rhs.rows_ ||        // and if dimensions
      lhs.cols_ != rhs.cols_           // are correct
  ) {
    return Matrix<T>();
  }

  OPERATION_INIT(lhs.rows_, lhs.cols_);

  cl::CommandQueue cmd_queue(context, device);
  cmd_queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                 cl::NDRange(result.rows_ * result.cols_));
  cl::finish();

  return result;
}

template <typename T>
Matrix<T> Matrix<T>::ExecuteBinaryScalarMatrixOperation(
    const T scalar,
    const Matrix<T> matrix,
    cl::Program program,
    const std::string &kernel_name
) {
  if (!matrix.matrix_) {
    return Matrix<T>();
  }

  auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
  auto device = context.getInfo<CL_CONTEXT_DEVICES>().front();

  Matrix<T> result = Matrix<T>::Empty(matrix.rows_, matrix.cols_);

  cl::Buffer inBuffer(
      context, CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_USE_HOST_PTR,
      matrix.size_, matrix.matrix_);
  cl::Buffer resultBuffer(
      context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY | CL_MEM_USE_HOST_PTR,
      result.size_, result.matrix_);

  cl::Kernel kernel(program, kernel_name.data());
  kernel.setArg(0, inBuffer);
  kernel.setArg(1, resultBuffer);
  kernel.setArg(2, scalar);

  cl::CommandQueue cmd_queue(context, device);
  cmd_queue.enqueueNDRangeKernel(kernel, cl::NullRange,
                                 cl::NDRange(result.rows_ * result.cols_));
  cl::finish();

  return result;
}

template <typename T>
std::ostream& operator<<(std::ostream &stream, Matrix<T> matrix) {
  for (std::size_t i = 0; i < matrix.rows(); i++) {
    for (std::size_t j = 0; j < matrix.cols(); j++) {
      std::cout << matrix[i][j] << ' ';
    }
    std::cout << std::endl;
  }
  return stream;
}

namespace {

template <typename T>
inline std::string GetTypeSuffix() {
  if constexpr (std::is_same<T, float>::value) {
    return std::string("_float");

  } else if constexpr (std::is_same<T, double>::value) {
    return std::string("_double");

  } else {
    static_assert(dependent_false<T>::value,
        "Type is not available in kernel");
  }
}

} // namespace

#undef LAZY_PROGRAM_INIT
#undef CL_INIT

#endif  // LAB5_MATRIX_H_
