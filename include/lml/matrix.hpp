#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace lml {

template <typename T>

class Matrix {
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<T> data_;

public:
  Matrix(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols) {}

  T &operator()(std::size_t row, std::size_t col) {
    return data_[row * cols_ + col];
  }
  const T &operator()(std::size_t row, std::size_t col) const {
    return data_[row * cols_ + col];
  }

  std::size_t cols() const { return cols_; }

  std::size_t rows() const { return rows_; }

  Matrix<T> operator+(const Matrix<T> &other) const {
    Matrix<T> result(rows_, cols_);

    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
      result.data_[i] = data_[i] + other.data_[i];
    }

    return result;
  }

  Matrix<T> operator-(const Matrix<T> &other) const {
    Matrix<T> result(rows_, cols_);

    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
      result.data_[i] = data_[i] - other.data_[i];
    }

    return result;
  }

  Matrix<T> operator+=(const Matrix<T> &other) {

    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
      data_[i] += other.data_[i];
    }

    return *this;
  }

  Matrix<T> operator-=(const Matrix<T> &other) {

    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
      data_[i] -= other.data_[i];
    }

    return *this;
  }

  Matrix<T> operator*=(const T scalar) {
    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
      data_[i] *= scalar;
    }

    return *this;
  }

  Matrix<T> transpose() const {
    Matrix<T> result(cols_, rows_);

    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        result(j, i) = (*this)(i, j);
      }
    }

    return result;
  }

  T detrecursive() const {

    if (rows_ == 1) {
      return (*this)(0, 0);
    } else if (rows_ == 2) {
      return (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);
    }

    T result = T{0};
    for (std::size_t k = 0; k < cols_; k++) {

      Matrix<T> A(rows_ - 1, cols_ - 1);

      for (std::size_t i = 0; i < rows_; i++) {
        for (std::size_t j = 1; j < cols_; j++) {
          if (i < k) {
            A(i, j - 1) = (*this)(i, j);
          } else if (i > k) {
            A(i - 1, j - 1) = (*this)(i, j);
          }
        }
      }
      T sgn = (k % 2 == 0) ? T{1} : T{-1};
      result += sgn * (*this)(0, k) * A.detGaus();
    }

    return result;
  }

  Matrix<T> rowswap(const std::size_t r1, const std::size_t r2) {

    for (std::size_t i = 0; i < cols_; i++) {
      T temp = (*this)(r1, i);
      (*this)(r1, i) = (*this)(r2, i);
      (*this)(r2, i) = temp;
    }
    return *this;
  }

  Matrix<T> rowmult(const T a, const std::size_t r) {
    for (std::size_t i = 0; i < cols_; i++) {
      (*this)(r, i) *= a;
    }
    return *this;
  }

  Matrix<T> rowadd(const std::size_t r1, const std::size_t r2, const T c) {
    for (std::size_t j = 0; j < cols_; j++) {
      (*this)(r2, j) += c * (*this)(r1, j);
    }
    return *this;
  }

  Matrix<T> rowreduce(const std::size_t r1, const std::size_t r2,
                      const std::size_t pivotcol) {
    if ((*this)(r2, pivotcol) == T{0}) {
      return *this;
    }
    T c = (*this)(r2, pivotcol) / (*this)(r1, pivotcol);
    (*this)(r2, pivotcol) = T{0};
    for (std::size_t j = pivotcol + 1; j < rows_; j++) {
      (*this)(r2, j) -= c * (*this)(r1, j);
    }
    return *this;
  }

  T detGaus() const {
    Matrix<T> A = (*this);
    T result = T{0};
    for (std::size_t k = 0; k < rows_; k++) {
      std::size_t pivot = k;

      for (std::size_t i = k + 1; i < rows_; ++i) {
        if (std::abs(A(i, k)) > std::abs(A(pivot, k))) {
          pivot = i;
        }
      }

      if (A(pivot, k) == T{0}) {
        return T{0};
      }

      if (pivot != k) {
        A.rowswap(k, pivot);
        result = -result;
      }

      for (std::size_t i = k + 1; i < rows_; ++i) {
        A.rowreduce(k, i, k);
      }

      result = A(k, k) * result;
    }
    return result;
  }

  T det() const {
    if (cols_ != rows_) {
      throw std::invalid_argument("Matrix must be square.");
    }

    return (*this).detGaus();
  }

}; // end class matrix

template <typename T>

Matrix<T> operator*(T scalar, const Matrix<T> &matrix) {
  Matrix<T> result = matrix;
  result *= scalar;
  return result;
};

template <typename T>

Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B) {
  if (A.cols() != B.rows()) {
    throw std::invalid_argument("Matrix dimensions don't match");
  }
  Matrix<T> result(B.cols(), A.rows());

  for (std::size_t i = 0; i < B.cols(); ++i) {
    for (std::size_t j = 0; j < A.rows(); ++j) {
      for (std::size_t k = 0; k < A.cols(); ++k) {
        result(i, j) += A(i, k) * B(k, j);
      }
    }
  }
  return result;
}

template <typename T>
Matrix<T> zeros(const std::size_t rows, const std::size_t cols) {
  return Matrix<T>(rows, cols);
}

template <typename T> Matrix<T> identity(const std::size_t n) {
  Matrix<T> result(n, n);

  for (std::size_t i = 0; i < n; ++i) {
    result(i, i) = T{1};
  }
  return result;
}

template <typename T>
std::vector<T> solveAxb(const Matrix<T> &M, const std::vector<T> &rhs) {
  Matrix<T> A = M;
  std::vector<T> b = rhs;

  if (A.rows() != b.size()) {
    throw std::invalid_argument("Matrix dimensions don't match");
  }

  std::size_t k = 0;
  std::size_t r = 0;
  std::vector<std::array<std::size_t, 2>> pivots;

  while (k < A.cols() &&
         r < A.rows()) { // bring A into unnormalized row echelon form

    std::size_t pivotrow = r;

    for (std::size_t i = r; i < A.rows(); i++) {

      if (std::abs(A(i, k)) > std::abs(A(pivotrow, k))) {
        pivotrow = i;
      }
    }

    if (A(pivotrow, k) == T{0}) {
      ++k;
      continue;
    }

    A.rowswap(pivotrow, r);
    std::swap(b[pivotrow], b[r]);
    pivots.push_back({r, k});

    for (std::size_t i = r + 1; i < A.rows(); i++) {
      T c = A(i, k) / A(r, k);
      A(i, k) = T{0};

      for (std::size_t j = k + 1; j < A.cols(); j++) {
        A(i, j) -= c * A(r, j);
      }
      b[i] -= c * b[r];
    }
    ++r;
    ++k;

  } // end while
  if (pivots.empty()) {
    for (std::size_t i = 0; i < b.size(); i++) {
      if (b[i] != T{0}) {
        return {};
      }
    }
    return std::vector<T>(A.cols(), T{0});
  }
  std::array<std::size_t, 2> lastpivot = pivots.back();
  std::size_t lastpivotrow = lastpivot[0];

  std::vector<T> result(
      A.cols(),
      T{0}); // result is initialized to 0, meaning free variables are set to 0

  for (std::size_t i = A.rows() - 1; i > lastpivotrow; i--) {
    if (b[i] != 0) {
      return std::vector<T>();
    }
  }
  for (std::size_t i = pivots.size(); i-- > 0;) { // perform back substitution
    auto [pivotrow, pivotcol] = pivots[i];
    result[pivotcol] = b[pivotrow] / A(pivotrow, pivotcol);
    for (std::size_t i1 = pivotrow; i1-- > 0;) {
      b[i1] -= result[pivotcol] * A(i1, pivotcol);
    }
  }
  return result;
}
} // namespace lml
