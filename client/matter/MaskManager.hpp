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
#include "client/HelperFunctions.h"
#include <stdexcept>

namespace eonc::mat {

struct AtomMask {
  const size_t ncol;
  Vector<size_t> atomIndices;
  Vector<bool> freeMask;
  Vector<size_t> freeIndices;

  AtomMask(const Vector<size_t> &indices, const Vector<bool> &mask,
           const size_t cols = 3)
      : ncol(cols),
        atomIndices(indices),
        freeMask(mask),
        freeIndices(
            atomIndices(eonc::helper_functions::filtered_indices(freeMask))) {
    if (indices.size() != mask.size()) {
      throw std::invalid_argument("Indices and mask must have the same size.");
    }
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
    for (auto index : free_indices_) {
      if (index >= matrix.rows()) {
        throw std::out_of_range("Index out of bounds");
      }
    }

    return getRows(matrix, free_indices_);
  }

private:
  Vector<size_t> free_indices_;

  // Utility function to get rows from a matrix based on indices
  template <typename Derived>
  auto getRows(const Eigen::MatrixBase<Derived> &matrix,
               const Vector<size_t> &indices) const {
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                  Derived::ColsAtCompileTime>
        filtered(indices.size(), matrix.cols());
    for (auto i = 0; i < indices.size(); ++i) {
      filtered.row(i) = matrix.row(indices[i]);
    }
    return filtered;
  }
};

} // namespace eonc::mat
