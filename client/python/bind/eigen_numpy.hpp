#pragma once
/**
 * Helpers: eOn row-major Eigen matrices <-> numpy ndarrays (nanobind).
 */
#include "Eigen.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

using NpF64 =
    nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>;
using NpI32 =
    nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu>;
using NpI64 =
    nb::ndarray<nb::numpy, int64_t, nb::c_contig, nb::device::cpu>;

inline AtomMatrix atom_matrix_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 2 || arr.shape(1) != 3) {
    throw std::invalid_argument(
        "expected float64 array of shape (n, 3), got ndim=" +
        std::to_string(arr.ndim()));
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  AtomMatrix m(n, 3);
  const double *p = arr.data();
  std::copy(p, p + n * 3, m.data());
  return m;
}

inline Matrix3d matrix3_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 2 || arr.shape(0) != 3 || arr.shape(1) != 3) {
    throw std::invalid_argument("expected float64 array of shape (3, 3)");
  }
  Matrix3d m;
  std::copy(arr.data(), arr.data() + 9, m.data());
  return m;
}

inline VectorXd vector_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 1) {
    throw std::invalid_argument("expected 1-d float64 array");
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  VectorXd v(n);
  std::copy(arr.data(), arr.data() + n, v.data());
  return v;
}

inline VectorXi vectori_from_numpy_i64(const NpI64 &arr) {
  if (arr.ndim() != 1) {
    throw std::invalid_argument("expected 1-d int64 array");
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  VectorXi v(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    v(i) = static_cast<int>(arr.data()[i]);
  }
  return v;
}

// Own a heap buffer so the ndarray remains valid after return.
template <typename EigenT>
inline nb::ndarray<nb::numpy, double, nb::c_contig>
matrix_to_numpy(const EigenT &m) {
  const size_t rows = static_cast<size_t>(m.rows());
  const size_t cols = static_cast<size_t>(m.cols());
  double *buf = new double[rows * cols];
  std::copy(m.data(), m.data() + rows * cols, buf);
  // capsule deletes the buffer
  nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<double *>(p); });
  return nb::ndarray<nb::numpy, double, nb::c_contig>(
      buf, {rows, cols}, owner);
}

template <typename Derived>
inline nb::ndarray<nb::numpy, double, nb::c_contig>
vector_to_numpy(const Eigen::MatrixBase<Derived> &v) {
  const size_t n = static_cast<size_t>(v.size());
  double *buf = new double[n];
  for (size_t i = 0; i < n; ++i) {
    buf[i] = static_cast<double>(v(static_cast<Eigen::Index>(i)));
  }
  nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<double *>(p); });
  return nb::ndarray<nb::numpy, double, nb::c_contig>(buf, {n}, owner);
}

inline nb::ndarray<nb::numpy, int64_t, nb::c_contig>
vectori_to_numpy(const VectorXi &v) {
  const size_t n = static_cast<size_t>(v.size());
  int64_t *buf = new int64_t[n];
  for (size_t i = 0; i < n; ++i) {
    buf[i] = static_cast<int64_t>(v(static_cast<Eigen::Index>(i)));
  }
  nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<int64_t *>(p); });
  return nb::ndarray<nb::numpy, int64_t, nb::c_contig>(buf, {n}, owner);
}

} // namespace eonc::pybind
