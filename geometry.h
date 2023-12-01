#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <valarray> //not use
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
      v.data[r] = (*this)[r][c];
    }

    return v;
  }

  void setRow(size_t r, const Vec<T, COL> &v) { (*this)[r] = v.data; }

  void setCol(size_t c, const Vec<T, ROW> &v) {
    for (size_t r = 0; r < ROW; r++) {
      (*this)[r][c] = v[r];
    }
  }

  // 取得删除某行和某列的子矩阵：子式
  inline Matrix<T, ROW - 1, COL - 1> getMinor(size_t row, size_t col) const {
    static_assert((ROW > 1) && (COL > 1), "Matrix too little to getMinor");

    Matrix<T, ROW - 1, COL - 1> mat;
    for (size_t r = 0; r < ROW - 1; r++) {
      for (size_t c = 0; c < COL - 1; c++) {
        mat[r][c] = (*this)[r < row ? r : r + 1][c < col ? c : c + 1];
      }
    }
    return mat;
  }

  // 取得转置矩阵
  inline Matrix<T, COL, ROW> transpose() const {
    Matrix<T, COL, ROW> mat;
    for (size_t r = 0; r < ROW; r++) {
      for (size_t c = 0; c < COL; c++)
        mat[c][r] = (*this)[r][c];
    }
    return mat;
  }
};

// Matrix运算
template <typename T, size_t ROW, size_t COL>
[[nodiscard]] constexpr inline bool
operator==(const Matrix<T, ROW, COL> &left, const Matrix<T, ROW, COL> &right) {
  return std::equal(left.begin(), left.end(), right.begin());
}

template <typename T, size_t ROW, size_t COL>
[[nodiscard]] constexpr inline bool
operator!=(const Matrix<T, ROW, COL> &left, const Matrix<T, ROW, COL> &right) {
  return !(left == right);
}

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
  Matrix<T, ROW, COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      mat[r][c] = left[r][c] - right[r][c];
  }
  return mat;
}

// mat*mat
template <typename T, size_t ROW, size_t COL, size_t NEW_COL>
inline Matrix<T, ROW, COL> operator*(const Matrix<T, ROW, COL> &left,
                                     const Matrix<T, COL, NEW_COL> &right) {
  Matrix<T, ROW, NEW_COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < NEW_COL; c++)
      mat[r][c] = vecDot(left.row(r), right.col(c));
  }
  return mat;
}

// mat*x
template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator*(const Matrix<T, ROW, COL> &left, T x) {
  Matrix<T, ROW, COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      mat[r][c] = left[r][c] * x;
  }
  return mat;
}

// x*mat
template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator*(T x, const Matrix<T, ROW, COL> &right) {
  return right * x;
}

// mat*v
template <typename T, size_t ROW, size_t COL>
inline Vec<T, ROW> operator*(const Matrix<T, ROW, COL> &left,
                             const Vec<T, COL> &right) {
  Vec<T, ROW> v;
  for (size_t i = 0; i < ROW; i++) {
    v[i] = vecDot(right, left.row(i));
  }
  return v;
}

// v*mat
template <typename T, size_t ROW, size_t COL>
inline Vec<T, COL> operator*(const Vec<T, ROW> &left,
                             const Matrix<T, ROW, COL> &right) {
  Vec<T, COL> v;
  for (size_t i = 0; i < COL; i++) {
    v[i] = vecDot(left, right.col(i));
  }
  return v;
}

template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator/(const Matrix<T, ROW, COL> &left, T x) {
  Matrix<T, ROW, COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      mat[r][c] = left[r][c] / x;
  }
  return mat;
}

template <typename T, size_t ROW, size_t COL>
inline Matrix<T, ROW, COL> operator/(T x, const Matrix<T, ROW, COL> &right) {
  Matrix<T, ROW, COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      mat[r][c] = x / right[r][c];
  }
  return mat;
}

template <typename T, size_t ROW, size_t COL>
constexpr inline T matCofactor(const Matrix<T, ROW, COL> &m, size_t row,
                               size_t col);

// Matrix函数
//行列式求值
template <typename T, size_t ROW, size_t COL>
constexpr inline T matDet(const Matrix<T, ROW, COL> &m) {
  static_assert(ROW == COL, "The rows and columns of a matrix are different.");

  if constexpr (ROW == 1) {
    return m[0][0];
  } else if constexpr (ROW == 2) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
  } else {
    T sum = 0;
    for (size_t i = 0; i < ROW; i++)
      sum += m[0][i] * matCofactor(m, 0, i);
    return sum;
  }
}

// 余子式
template <typename T, size_t ROW, size_t COL>
constexpr inline T matCofactor(const Matrix<T, ROW, COL> &m, size_t row,
                               size_t col) {
  static_assert(ROW == COL, "The rows and columns of a matrix are different.");
  if constexpr (ROW == 1) {
    return 0;
  } else {
    return matDet(m.getMinor(row, col)) * (((row + col) % 2) ? -1 : 1);
  }
}

// 伴随矩阵
template <typename T, size_t ROW, size_t COL>
constexpr inline T matAdjoint(const Matrix<T, ROW, COL> &m) {
  Matrix<T, ROW, COL> mat;
  for (size_t r = 0; r < ROW; r++) {
    for (size_t c = 0; c < COL; c++)
      mat[r][c] = matrix_cofactor(m, c, r);
  }
  return mat;
}

// 求逆矩阵
template <typename T, size_t ROW, size_t COL>
constexpr inline T matInvert(const Matrix<T, ROW, COL> &m) {
  Matrix<T, ROW, COL> mat = matAdjoint(m);
  T det = vecDot(m.Row(0), mat.Col(0));
  return mat / det;
}

using Vec2f = Vec<float, 2>;
using Vec2i = Vec<int32_t, 2>;
using Vec3f = Vec<float, 3>;
using Vec3i = Vec<int32_t, 3>;
using Vec4f = Vec<float, 4>;
using Vec4i = Vec<int32_t, 4>;

using Mat3x3f = Matrix<float, 3, 3>;
using Mat4x4f = Matrix<float, 4, 4>;

// 3D数学运算
constexpr inline Mat4x4f matIdentity() {
  return Mat4x4f{std::array<float, 4>{1, 0, 0, 0},
                 {0, 1, 0, 0},
                 {0, 0, 1, 0},
                 {0, 0, 0, 1}};
}

constexpr inline Mat4x4f matTranslate(float x, float y, float z) {
  return Mat4x4f{std::array<float, 4>{1, 0, 0, 0},
                 {0, 1, 0, 0},
                 {0, 0, 1, 0},
                 {x, y, z, 1}};
}

constexpr inline Mat4x4f matScale(float x, float y, float z) {
  return Mat4x4f{std::array<float, 4>{x, 0, 0, 0},
                 {0, y, 0, 0},
                 {0, 0, z, 0},
                 {0, 0, 0, 1}};
}

constexpr inline Mat4x4f matRotate(float x, float y, float z, float theta) {
  float qsin = static_cast<float>(std::sin(theta * 0.5f));
  float qcos = static_cast<float>(std::cos(theta * 0.5f));
  float w = qcos;
  Vec3f vec = vecNormalize(Vec3f{x, y, z});
  x = vec.x * qsin;
  y = vec.y * qsin;
  z = vec.z * qsin;

  return Mat4x4f{
      std::array<float, 4>{1 - 2 * y * y - 2 * z * z, 2 * x * y + 2 * w * z,
                           2 * x * z - 2 * w * y, 0},
      std::array<float, 4>{2 * x * y - 2 * w * z, 1 - 2 * x * x - 2 * z * z,
                           2 * y * z + 2 * w * x, 0},
      std::array<float, 4>{2 * x * z + 2 * w * y, 2 * y * z - 2 * w * x,
                           1 - 2 * x * x - 2 * y * y, 0},
      std::array<float, 4>{0, 0, 0, 1}};
}

inline Mat4x4f matLookat(const Vec3f &eye, const Vec3f &at, const Vec3f &up) {
  Vec3f zaxis = vecNormalize(at - eye);
  Vec3f xaxis = vecNormalize(vecCross(up, zaxis));
  Vec3f yaxis = vecCross(zaxis, xaxis);

  return Mat4x4f{
      std::array<float, 4>{xaxis.x, yaxis.x, zaxis.x, 0},
      {xaxis.y, yaxis.y, zaxis.y, 0},
      {xaxis.z, yaxis.z, zaxis.z, 0},
      {-vecDot(eye, xaxis), -vecDot(eye, yaxis), -vecDot(eye, zaxis), 1}};
}

inline Mat4x4f matPerspecTive(float fovy, float aspect, float zn, float zf) {
  float fax = 1.0f / static_cast<float>(tan(fovy * 0.5f));
  Mat4x4f mat{};
  mat[0][0] = fax / aspect;
  mat[1][1] = fax;
  mat[2][2] = zf / (zf - zn);
  mat[3][2] = -zn * zf / (zf - zn);
  mat[2][3] = 1;
  return mat;
}

} // namespace raster
