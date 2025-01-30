#ifndef MATRIX_H
#define MATRIX_SQUARE_MATRIX_IMPLEMENTED

#include <stdexcept>
#include <iostream>
#include <initializer_list>

class MatrixIsDegenerateError : public std::runtime_error {
 public:
  MatrixIsDegenerateError() : std::runtime_error("MatrixIsDegenerateError") {
  }
};

class MatrixOutOfRange : public std::out_of_range {
 public:
  MatrixOutOfRange() : std::out_of_range("MatrixOutOfRange") {
  }
};

template <typename T, size_t N, size_t M>
class Matrix {
 private:
 public:
  T matrix_[N][M];

  void ShowMatrix() const {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        std::cout << matrix_[i][j] << " ";
      }
      std::cout << '\n';
    }
  }

  size_t RowsNumber() const {
    return N;
  }

  size_t ColumnsNumber() const {
    return M;
  }

  T& operator()(const size_t n, const size_t m) {
    return matrix_[n][m];
  }

  const T& operator()(const size_t n, const size_t m) const {
    return matrix_[n][m];
  }

  T& At(const size_t n, const size_t m) {
    if (n >= N ||m >= M) {
      throw MatrixOutOfRange{};
    }

    return matrix_[n][m];
  }

  const T& At(const size_t n, const size_t m) const {
    if (n >= N ||m >= M) {
      throw MatrixOutOfRange{};
    }

    return matrix_[n][m];
  }

  Matrix& operator+=(const Matrix& more) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matrix_[i][j] += more(i, j);
      }
    }
    return *this;
  }

  Matrix& operator-=(const Matrix& more) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matrix_[i][j] -= more(i, j);
      }
    }
    return *this;
  }

  Matrix& operator*=(const Matrix<T, M, M>& more) {
    Matrix third;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        third(i, j) = 0;
        for (size_t index = 0; index < M; ++index) {
          third(i, j) += matrix_[i][index] * more(index, j);
        }
      }
    }
    *this = third;
    return *this;
  }

  friend Matrix operator+(const Matrix& first, const Matrix& second) {
    Matrix third = first;
    third += second;
    return third;
  }

  friend Matrix operator-(const Matrix& first, const Matrix& second) {
    Matrix third = first;
    third -= second;
    return third;
  }

  bool operator==(const Matrix& more) const {
    bool flag = true;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        if (matrix_[i][j] != more(i, j)) {
          flag = false;
        }
      }
    }
    return flag;
  }

  bool operator!=(const Matrix& more) const {
    return !(*this == more);
  }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& tmp_matrix) {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        if (j > 0) {
          os << ' ';
        }
        os << tmp_matrix(i, j);
      }
      os << '\n';
    }
    return os;
  }

  friend std::istream& operator>>(std::istream& is, Matrix& tmp_matrix) {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        is >> tmp_matrix(i, j);
      }
    }
    return is;
  }
};

template <typename T, size_t N, size_t M, typename U>
Matrix<T, N, M>& operator*=(Matrix<T, N, M>& matrix, const U& val) {
  auto value = static_cast<T>(val);
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++) {
      matrix.matrix_[i][j] *= value;
    }
  }
  return matrix;
}

template <typename T, size_t N, size_t M, typename U>
Matrix<T, N, M> operator*(const Matrix<T, N, M>& matrix, const U& val) {
  Matrix<T, N, M> third = matrix;
  third *= val;
  return third;
}

template <typename T, size_t N, size_t M, typename U>
Matrix<T, N, M> operator*(const U& val, const Matrix<T, N, M>& matrix) {
  return matrix * val;
}

template <typename T, size_t N, size_t M, typename U>
Matrix<T, N, M>& operator/=(Matrix<T, N, M>& tmp_matrix, const U& val) {
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < M; j++) {
      tmp_matrix(i, j) /= val;
    }
  }
  return tmp_matrix;
}

template <typename T, size_t N, size_t M, typename U>
Matrix<T, N, M> operator/(const Matrix<T, N, M>& tmp_matrix, const U& val) {
  Matrix<T, N, M> third = tmp_matrix;
  third /= val;
  return third;
}

template <typename T, size_t N, size_t M, size_t K>
Matrix<T, N, K> operator*(const Matrix<T, N, M>& first, const Matrix<T, M, K>& second) {
  Matrix<T, N, K> third;
  for (size_t n = 0; n < N; n++) {
    for (size_t k = 0; k < K; k++) {
      third.matrix_[n][k] = 0;
      for (size_t m = 0; m < M; m++) {
        third.matrix_[n][k] += first.matrix_[n][m] * second.matrix_[m][k];
      }
    }
  }
  return third;
}

template <typename T, size_t N, size_t M>
Matrix<T, M, N> GetTransposed(const Matrix<T, N, M>& to_trans_matrix) {
  Matrix<T, M, N> third;

  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < N; j++) {
      third(i, j) = to_trans_matrix(j, i);
    }
  }
  return third;
}

#endif
