/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once
#include "client/Eigen.h"
#include <range/v3/all.hpp>
#include <stdexcept>
#include <vector>

namespace eonc::mat {

template <typename Mask> auto filtered_indices(Mask const &mask) {
  using ranges::begin;
  using ranges::end;
  using ranges::views::filter;
  using ranges::views::iota;

  int const N = mask.size();
  auto c = iota(0, N) | filter([&mask](auto const &i) { return mask[i]; });
  return std::vector(begin(c), end(c));
}

struct AtomMask {
  Vector<size_t> atomIndices;
  Vector<bool> freeMask;
  std::vector<int> freeIndices;

  AtomMask(const Vector<size_t> &indices, const Vector<bool> &mask)
      : atomIndices(indices),
        freeMask(mask) {
    if (indices.size() != mask.size()) {
      throw std::invalid_argument("Indices and mask must have the same size.");
    }
    freeIndices = filtered_indices(freeMask);
  }
};

class FreeFixer {
public:
  FreeFixer(const AtomMask &mask)
      : free_indices_(mask.freeIndices) {}

  // Method to apply the mask to a specific matrix (e.g., positions)
  template <typename Derived>
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                Derived::ColsAtCompileTime>
  operator()(const Eigen::MatrixBase<Derived> &matrix) const {
    return this->getRows(matrix, free_indices_);
  }

private:
  std::vector<int> free_indices_;

  // Utility function to get rows from a matrix based on indices
  template <typename Derived>
  auto getRows(const Eigen::MatrixBase<Derived> &matrix,
               const std::vector<int> &indices) const {
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                  Derived::ColsAtCompileTime>
        filtered(indices.size(), matrix.cols());
    for (size_t i = 0; i < indices.size(); ++i) {
      filtered.row(i) = matrix.row(indices[i]);
    }
    return filtered;
  }
};

} // namespace eonc::mat
