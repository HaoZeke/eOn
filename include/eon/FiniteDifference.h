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

#include "Eigen.h"

#include <cctype>
#include <cmath>
#include <string>
#include <string_view>

namespace eonc {

/// Real finite-difference scheme for the assembled Hessian and for
/// Lanczos/Davidson Hessian-vector products. Complex-step is not a scheme:
/// cutoffs and external potentials are not holomorphic.
enum class FdScheme { OneSided, Central, Fourth };

inline FdScheme parseFdScheme(std::string_view scheme) {
  std::string lower;
  lower.reserve(scheme.size());
  for (unsigned char c : scheme) {
    lower.push_back(static_cast<char>(std::tolower(c)));
  }
  if (lower == "central") {
    return FdScheme::Central;
  }
  if (lower == "fourth" || lower == "fourth_order" || lower == "central4") {
    return FdScheme::Fourth;
  }
  return FdScheme::OneSided;
}

/// Step for an assembled Hessian.
/// ``Main.finite_difference`` is the optimizer step and defaults to 0.01 Å.
/// A one-sided difference at that length biased a soft eigenvalue by 40
/// percent. The assembled Hessian uses the requested step, and never a step
/// above 1e-3 Å.
inline double hessianFiniteDifference(double requested) {
  if (!(requested > 0.0) || !std::isfinite(requested)) {
    return requested;
  }
  constexpr double kMaxHessianStep = 1e-3;
  return requested > kMaxHessianStep ? kMaxHessianStep : requested;
}

/// Derivative of a sampled force map along one real step of length ``dr``.
/// Fourth-order central:
/// ``[-f(+2h) + 8 f(+h) - 8 f(-h) + f(-2h)] / (12 h)``.
template <class M>
M fdForceDerivative(FdScheme scheme, double dr, const M &f0, const M &fPlus,
                    const M &fMinus, const M &fPlus2, const M &fMinus2) {
  M slope = fPlus;
  switch (scheme) {
  case FdScheme::Fourth:
    slope = (-fPlus2 + 8.0 * fPlus - 8.0 * fMinus + fMinus2) / (12.0 * dr);
    break;
  case FdScheme::Central:
    slope = (fPlus - fMinus) / (2.0 * dr);
    break;
  case FdScheme::OneSided:
    slope = (fPlus - f0) / dr;
    break;
  }
  return slope;
}

/// Energy Hessian-vector product ``-dF``. ``eval(scale)`` is the force at
/// ``x + scale * dr * v``. One-sided and central read only the samples they
/// need; fourth reads ``scale`` in ``{1, -1, 2, -2}``.
template <class Eval>
VectorXd fdHessianVector(FdScheme scheme, double dr, const VectorXd &force0,
                         Eval &&eval) {
  const VectorXd zero = VectorXd::Zero(force0.size());
  switch (scheme) {
  case FdScheme::Fourth: {
    const VectorXd fPlus = eval(1.0);
    const VectorXd fMinus = eval(-1.0);
    const VectorXd fPlus2 = eval(2.0);
    const VectorXd fMinus2 = eval(-2.0);
    return -fdForceDerivative(scheme, dr, force0, fPlus, fMinus, fPlus2,
                              fMinus2);
  }
  case FdScheme::Central: {
    const VectorXd fPlus = eval(1.0);
    const VectorXd fMinus = eval(-1.0);
    return -fdForceDerivative(scheme, dr, force0, fPlus, fMinus, zero, zero);
  }
  case FdScheme::OneSided:
    break;
  }
  const VectorXd fPlus = eval(1.0);
  return -fdForceDerivative(FdScheme::OneSided, dr, force0, fPlus, zero, zero,
                            zero);
}

} // namespace eonc
