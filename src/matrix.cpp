#include <cstddef>
#include <vector>

#include "lml/matrix.hpp"

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
};
