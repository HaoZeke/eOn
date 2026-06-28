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
// Leng et al. JCP 138, 094110 (2013) Locally Optimal Rotation (Algorithm I).
// FD Hessian-vector products match Lanczos/ImprovedDimer energy-Hessian convention:
//   H v ≈ -(F(x+δv) - F(x)) / δ  with F = forces (not gradients).

#include "LORRotation.h"
#include "HelperFunctions.h"
#include "SafeMath.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

using namespace eonc::helpers;

LORRotation::LORRotation(std::shared_ptr<Matter> matter,
                         const Parameters &params,
                         std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  auto x1Pot = (pot->needsPerImageInstance() && params.main_options.parallel)
                   ? eonc::helpers::makePotential(params)
                   : pot;
  x0 = std::make_shared<Matter>(pot, params);
  x1 = std::make_shared<Matter>(x1Pot, params);
  *x0 = *matter;
  *x1 = *matter;
  totalForceCalls = 0;
  statsRotations = 0;
  eigenvector.resize(3 * matter->numberOfAtoms());
  eigenvector.setZero();
}

VectorXd LORRotation::hessianVector(const VectorXd & /*g0 unused*/,
                                    const VectorXd &x0_r, const VectorXd &v,
                                    const VectorXd &freeMask, double delta) {
  // Force-based FD (same sign convention as Lanczos::compute).
  VectorXd dir = v.array() * freeMask.array();
  const double nrm = dir.norm();
  if (nrm < 1e-14) {
    return VectorXd::Zero(v.size());
  }
  dir /= nrm;

  VectorXd F0 = x0->getForcesV();
  x1->setPositionsV(x0_r + delta * dir);
  VectorXd F1 = x1->getForcesV();
  totalForceCalls += 1;

  // H_energy · dir  ≈  -(F1 - F0) / delta
  VectorXd Hv = -(F1 - F0) / delta;
  Hv = Hv.array() * freeMask.array();
  // Scale so Hv is for unit v (dir is unit); caller passes general v — return H·v
  // with v not necessarily unit: H·v = ||v|| * H·(v/||v||) but we used unit dir
  // from masked v, so Hv is H·unit. Restore magnitude of projected v:
  return Hv * nrm;
}

void LORRotation::compute(std::shared_ptr<Matter> matter,
                          AtomMatrix initialDirectionAtomMatrix) {
  const long nAtoms = matter->numberOfAtoms();
  const int dim = static_cast<int>(3 * nAtoms);
  const VectorXd freeMask = matter->getFreeV();

  VectorXd N = VectorXd::Map(initialDirectionAtomMatrix.data(), dim);
  N = N.array() * freeMask.array();
  if (N.norm() < 1e-10) {
    N.setRandom();
    N = N.array() * freeMask.array();
  }
  N.normalize();

  *x0 = *matter;
  *x1 = *matter;
  const VectorXd x0_r = x0->getPositionsV();
  const double delta = params.main_options.finiteDifference;

  // Paper Baker runs use ~5–20 rotations; never use geometry max_iterations.
  const int rotmax = std::clamp(
      static_cast<int>(params.dimer_options.rotations_max > 0
                           ? params.dimer_options.rotations_max
                           : 20),
      1, 50);

  // Paper stops rotation when ||F_perp|| is small (they use ~0.1 eV/Å in VASP
  // tests). Align with torque_min / converged_angle without exiting on noise.
  const double residualTol = std::max(
      1e-3, std::min(0.1, params.dimer_options.torque_min));

  curvatureHistory.clear();
  convergedOnResidual = false;
  statsRotations = 0;
  totalForceCalls = 0;

  // Center force for bookkeeping (one call; also used if we ever need F0-only paths)
  (void)x0->getForcesV();
  totalForceCalls += 1;

  auto applyMask = [&](VectorXd &v) {
    v = v.array() * freeMask.array();
  };

  auto unitize = [&](VectorXd &v) -> double {
    applyMask(v);
    const double n = v.norm();
    if (n > 1e-14) {
      v /= n;
    }
    return n;
  };

  // H N — one FD along N (unit)
  VectorXd HN = hessianVector(VectorXd(), x0_r, N, freeMask, delta);
  applyMask(HN);
  double CN = N.dot(HN);
  curvatureHistory.push_back(CN);

  VectorXd F = HN - CN * N; // F_⊥ = (I - N Nᵀ) H N
  applyMask(F);
  double Fnorm = F.norm();

  QUILL_LOG_INFO(log,
                 "[LOR] iter=0 ||F_perp||={:.6e} C_N={:.6f} (Algorithm I start)",
                 Fnorm, CN);

  if (Fnorm < residualTol) {
    convergedOnResidual = true;
    eigenvalue = CN;
    eigenvector = N;
    QUILL_LOG_INFO(log, "[LOR] converged on residual at start");
    return;
  }

  VectorXd Theta = F / Fnorm;
  VectorXd HTheta = hessianVector(VectorXd(), x0_r, Theta, freeMask, delta);
  applyMask(HTheta);

  // 2×2 Rayleigh–Ritz in span{N, Θ} (orthonormal)
  Eigen::Matrix2d A2;
  A2(0, 0) = N.dot(HN);
  A2(0, 1) = N.dot(HTheta);
  A2(1, 0) = A2(0, 1);
  A2(1, 1) = Theta.dot(HTheta);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2(A2);
  Eigen::Vector2d coeffs = es2.eigenvectors().col(0); // smallest eigenpair
  double a = coeffs(0);
  double b = coeffs(1);

  // N^(2) = a N + b Θ  (normalize; scale HN by same factor — critical for translation)
  VectorXd Nlin = a * N + b * Theta;
  const double nN = unitize(Nlin);
  VectorXd HNlin = a * HN + b * HTheta;
  if (nN > 1e-14) {
    HNlin /= nN;
  }
  N = Nlin;
  HN = HNlin;

  // P^(2) = Θ, HP = HΘ (unit Θ already)
  VectorXd P = Theta;
  VectorXd HP = HTheta;

  CN = N.dot(HN);
  curvatureHistory.push_back(CN);
  F = HN - CN * N;
  applyMask(F);
  Fnorm = F.norm();
  statsRotations = 1;

  QUILL_LOG_INFO(log,
                 "[LOR] iter=1 (2x2) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} b={:.4f}",
                 Fnorm, CN, a, b);

  for (int k = 2; k <= rotmax; ++k) {
    if (Fnorm < residualTol) {
      convergedOnResidual = true;
      QUILL_LOG_INFO(log, "[LOR] converged residual iter={} ||F||={:.6e}", k,
                     Fnorm);
      break;
    }

    Theta = F / Fnorm;
    // One new force evaluation: H · Θ
    HTheta = hessianVector(VectorXd(), x0_r, Theta, freeMask, delta);
    applyMask(HTheta);

    // Keep historical P as unit vector (paper: may be non-orthogonal to N, Θ)
    VectorXd P3 = P;
    const double pNrm = unitize(P3);
    VectorXd HP3 = HP;
    if (pNrm > 1e-14 && std::abs(pNrm - 1.0) > 1e-12) {
      // P already unit from prior step; HP consistent with unit P
    }
    (void)pNrm;

    Eigen::Matrix3d A3 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d B3 = Eigen::Matrix3d::Identity();
    const VectorXd basis[3] = {N, Theta, P3};
    const VectorXd Hbasis[3] = {HN, HTheta, HP3};
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
        A3(i, j) = basis[i].dot(Hbasis[j]);
        A3(j, i) = A3(i, j);
        if (i != j) {
          B3(i, j) = basis[i].dot(basis[j]);
          B3(j, i) = B3(i, j);
        }
      }
    }

    // Symmetrize A for numerical noise (FD H not exactly symmetric in practice)
    A3 = 0.5 * (A3 + A3.transpose());

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> ges(A3, B3);
    if (ges.info() != Eigen::Success) {
      QUILL_LOG_WARNING(log,
                        "[LOR] 3x3 generalized eigen failed at iter={}; stop",
                        k);
      break;
    }
    const Eigen::Vector3d c3 = ges.eigenvectors().col(0);
    a = c3(0);
    b = c3(1);
    const double c = c3(2);

    // N^{k+1} = a N + b Θ + c P ; P^{k+1} = b Θ + c P  (Algorithm I)
    VectorXd Nnew = a * N + b * Theta + c * P3;
    VectorXd Pnew = b * Theta + c * P3;
    VectorXd HNnew = a * HN + b * HTheta + c * HP3;
    VectorXd HPnew = b * HTheta + c * HP3;

    const double nNew = unitize(Nnew);
    if (nNew > 1e-14) {
      HNnew /= nNew;
    }
    const double pNew = unitize(Pnew);
    if (pNew > 1e-14) {
      HPnew /= pNew;
    }

    N = Nnew;
    P = Pnew;
    HN = HNnew;
    HP = HPnew;

    const double CNnew = N.dot(HN);
    if (!curvatureHistory.empty() && CNnew > curvatureHistory.back() + 1e-3) {
      QUILL_LOG_INFO(log,
                     "[LOR] curvature increased iter={} C_prev={:.6f} C_new={:.6f} "
                     "(FD/anharmonic; continue)",
                     k, curvatureHistory.back(), CNnew);
    }
    curvatureHistory.push_back(CNnew);
    CN = CNnew;

    F = HN - CN * N;
    applyMask(F);
    Fnorm = F.norm();
    statsRotations = k;

    QUILL_LOG_INFO(log,
                   "[LOR] iter={} (3x3) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} "
                   "b={:.4f} c={:.4f}",
                   k, Fnorm, CN, a, b, c);

    if (Fnorm < residualTol) {
      convergedOnResidual = true;
      break;
    }
  }

  // Final FD refresh of H·N so eigenvalue matches a true FD curvature (reduces
  // accumulated translation error on anharmonic PES before climb).
  HN = hessianVector(VectorXd(), x0_r, N, freeMask, delta);
  applyMask(HN);
  CN = N.dot(HN);
  curvatureHistory.push_back(CN);

  eigenvalue = CN;
  eigenvector = N;
  QUILL_LOG_INFO(log,
                 "[LOR] done rotations={} force_calls={} C_N={:.6f} "
                 "converged_residual={}",
                 statsRotations, totalForceCalls, eigenvalue,
                 convergedOnResidual);
}
