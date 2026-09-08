#include <cstddef>
#include <vector>

template <typename T> class Matrix {
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<T> data_;

public:
  Matrix(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols) {}
};
