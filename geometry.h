#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <variant>

namespace raster {
/***************************************************向量*******************************************/
// #if (__cplusplus >= 201703)
// template <typename T, size_t N> struct Vec : public std::array<T, N> {};
// #else
template <typename T, size_t N> struct Vec {
  std::array<T, N> data{};

  T &operator[](size_t i) noexcept { return data[i]; }
  constexpr const T &operator[](size_t i) const noexcept { return data[i]; }
};
//#endif

template <typename T> struct Vec<T, 2> {
  union {
    std::array<T, 2> data{};
    struct {
      T x;
      T y;
    };
  };
  T &operator[](size_t i) noexcept { return data[i]; }
  constexpr const T &operator[](size_t i) const noexcept { return data[i]; }
};

template <typename T> struct Vec<T, 3> {
  union {
    std::array<T, 3> data{};
    struct {
      T x;
      T y;
      T z;
    };
  };
  T &operator[](size_t i) noexcept { return data[i]; }
  constexpr const T &operator[](size_t i) const noexcept { return data[i]; }
};

template <typename T> struct Vec<T, 4> {
  union {
    std::array<T, 4> data{};
    struct {
      T x;
      T y;
      T z;
      T w;
    };
  };
  T &operator[](size_t i) noexcept { return data[i]; }
  constexpr const T &operator[](size_t i) const noexcept { return data[i]; }
};

// Vec运算
template <typename T, size_t N>
[[nodiscard]] constexpr inline bool operator==(const Vec<T, N> left,
                                               const Vec<T, N> right) {
  return left.data == right.data;
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline bool operator!=(const Vec<T, N> left,
                                               const Vec<T, N> right) {
  return !(left == right);
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator+(const Vec<T, N> left,
                                         const Vec<T, N> right) {
  Vec<T, N> v;
  std::transform(left.data.cbegin(), left.data.cend(), right.data.cbegin(),
                 v.data.begin(), std::plus<>{});
  return v;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator-(const Vec<T, N> left,
                                         const Vec<T, N> right) {
  Vec<T, N> v;
  std::transform(left.data.cbegin(), left.data.cend(), right.data.cbegin(),
                 v.data.begin(), std::minus<>{});
  return v;
}

//不是点乘，不是叉乘
template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator*(const Vec<T, N> left,
                                         const Vec<T, N> right) {
  Vec<T, N> v;
  std::transform(left.data.cbegin(), left.data.cend(), right.data.cbegin(),
                 v.data.begin(), std::multiplies<>{});
  return v;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator/(const Vec<T, N> left,
                                         const Vec<T, N> right) {
  Vec<T, N> v;
  std::transform(left.data.cbegin(), left.data.cend(), right.data.cbegin(),
                 v.data.begin(), std::divides<>{});
  return v;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator*(const Vec<T, N> src_v, T x) {
  Vec<T, N> v;
  std::transform(src_v.data.cbegin(), src_v.data.cend(), v.data.begin(),
                 [x](T element) { return element * x; });
  return v;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator*(T x, const Vec<T, N> src_v) {
  Vec<T, N> v;
  return v * x;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator/(const Vec<T, N> src_v, T x) {
  Vec<T, N> v;
  std::transform(src_v.data.cbegin(), src_v.data.cend(), v.data.begin(),
                 [x](T element) { return element / x; });
  return v;
}

template <typename T, size_t N>
[[nodiscard]] inline Vec<T, N> operator/(T x, const Vec<T, N> src_v) {
  Vec<T, N> v;
  std::transform(src_v.data.cbegin(), src_v.data.cend(), v.data.begin(),
                 [x](T element) { return x / element; });
  return v;
}

// Vec函数
//|v|^2
template <typename T, size_t N> inline T vecLenSquare(const Vec<T, N> &v) {
  T sum = std::accumulate(v.data.cbegin(), v.data.cend(), 0,
                          [](T a, T b) { return a + (b * b); });
  return sum;
}

//|a|
template <typename T, size_t N> inline T vecLen(const Vec<T, N> &v) {
  return std::sqrt(vecLenSquare(v));
}

template <typename T, size_t N>
inline Vec<T, N> vecNormalize(const Vec<T, N> &v) {
  return v / vecLen(v);
}

template <typename T, size_t N>
inline T vecDot(const Vec<T, N> &v1, const Vec<T, N> &v2) {
  return std::inner_product(v1.data.cbegin(), v1.data.cend(), v2.data.cbegin(),
                            0);
}

template <typename T, size_t N>
inline auto vecCross(const Vec<T, N> &v1, const Vec<T, N> &v2) {
  if constexpr (N == 2) {
    return T{v1.x * v2.y - v1.y * v2.x};
  } else if constexpr (N == 3) {
    return Vec<T, 3>{v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                     v1.x * v2.y - v1.y * v2.x};
  } else if constexpr (N == 4) {
    return Vec<T, 4>{v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                     v1.x * v2.y - v1.y * v2.x, v1.w};
  } else {
    static_assert(!sizeof(std::decay_t<T>),
                  "only support arithmetic type except char type");
  }
}

/***************************************************矩阵*******************************************/

// using Matrix = std::array<std::array<T, ROW>, COL>;
template <typename T, size_t ROW, size_t COL>
struct Matrix : public std::array<std::array<T, ROW>, COL> {
  // 取一行
  Vec<T, COL> row(size_t r) const { return Vec<T, COL>{this->at(r)}; }
  // 取一列
  Vec<T, ROW> col(size_t c) const {
    Vec<T, ROW> v;
    for (size_t r = 0; r < ROW; r++) {
      v.data[r] = *this[r][c];
    }

    return v;
  }

  void setRow(size_t r, const Vec<T, COL> &v) { *this[r] = v.data; }

  void setCol(size_t c, const Vec<T, ROW> &v) {
    for (size_t r = 0; r < ROW; r++) {
      *this[r][c] = v[r];
    }
  }

  // 取得删除某行和某列的子矩阵：子式
  inline Matrix<T, ROW - 1, COL - 1> getMinor(size_t row, size_t col) const {
    static_assert((ROW > 1) && (COL > 1), "Matrix too little to getMinor");

    Matrix<T, ROW - 1, COL - 1> mat;
    for (size_t r = 0; r < ROW - 1; r++) {
      for (size_t c = 0; c < COL - 1; c++) {
        mat[r][c] = *this[r < row ? r : r + 1][c < col ? c : c + 1];
      }
    }
    return mat;
  }

  // 取得转置矩阵
  inline Matrix<T, COL, ROW> transpose() const {
    Matrix<T, COL, ROW> mat;
    for (size_t r = 0; r < ROW; r++) {
      for (size_t c = 0; c < COL; c++)
        mat[c][r] = *this[r][c];
    }
    return mat;
  }
};

// Matrix运算
template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator+(const Matrix<T, ROW, COL> &left,
                                     const Matrix<T, ROW, COL> &right) {
  Matrix<T, ROW, COL> sum;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      sum[r][c] = left[r][c] + right[r][c];
  }
  return sum;
}

template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator-(const Matrix<T, ROW, COL> &left,
                                     const Matrix<T, ROW, COL> &right) {
  Matrix<T, ROW, COL> sum;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      sum[r][c] = left[r][c] - right[r][c];
  }
  return sum;
}

template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator*(const Matrix<T, ROW, COL> &left,
                                     const Matrix<T, ROW, COL> &right) {
  Matrix<T, ROW, COL> sum;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      sum[r][c] = vecDot(left.row(r), right.col(c));
  }
  return sum;
}

} // namespace raster
