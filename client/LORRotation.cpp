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
#include <limits>
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
  const double residualTol = std::max(1e-3, params.dimer_options.torque_min);
  auto relativeResidual = [](double fnorm, double cn) {
    return fnorm / (std::abs(cn) + 1.0);
  };

  curvatureHistory.clear();
  convergedOnResidual = false;
  double bestCN = std::numeric_limits<double>::infinity();
  VectorXd bestN = N;
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
  if (CN < bestCN) {
    bestCN = CN;
    bestN = N;
  }

  VectorXd F = HN - CN * N; // F_⊥ = (I - N Nᵀ) H N
  applyMask(F);
  double Fnorm = F.norm();

  QUILL_LOG_INFO(log,
                 "[LOR] iter=0 ||F_perp||={:.6e} C_N={:.6f} (Algorithm I start)",
                 Fnorm, CN);

  if (relativeResidual(Fnorm, CN) < residualTol) {
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
  if (CN < bestCN) {
    bestCN = CN;
    bestN = N;
  }
  F = HN - CN * N;
  applyMask(F);
  Fnorm = F.norm();
  statsRotations = 1;

  QUILL_LOG_INFO(log,
                 "[LOR] iter=1 (2x2) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} b={:.4f}",
                 Fnorm, CN, a, b);

  for (int k = 2; k <= rotmax; ++k) {
    if (relativeResidual(Fnorm, CN) < residualTol) {
      convergedOnResidual = true;
      QUILL_LOG_INFO(log, "[LOR] converged residual iter={} ||F||={:.6e}", k,
                     Fnorm);
      break;
    }

    Theta = F / Fnorm;
    // One new force evaluation: H · Θ
    HTheta = hessianVector(VectorXd(), x0_r, Theta, freeMask, delta);
    applyMask(HTheta);

    // Orthonormalize P against N, Θ for a well-conditioned 3×3 (paper allows
    // non-ortho B; FD noise + GEP then blows up coefficients).
    VectorXd P3 = P - N.dot(P) * N - Theta.dot(P) * Theta;
    applyMask(P3);
    VectorXd HP3 = HP;
    const double pNrm = P3.norm();
    if (pNrm < 1e-8) {
      // Degenerate P — fall back to 2×2 in span{N, Θ}
      Eigen::Matrix2d A2b;
      A2b(0, 0) = N.dot(HN);
      A2b(0, 1) = N.dot(HTheta);
      A2b(1, 0) = A2b(0, 1);
      A2b(1, 1) = Theta.dot(HTheta);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2b(A2b);
      Eigen::Vector2d c2 = es2b.eigenvectors().col(0);
      a = c2(0);
      b = c2(1);
      VectorXd Nnew = a * N + b * Theta;
      VectorXd HNnew = a * HN + b * HTheta;
      const double nNew = unitize(Nnew);
      if (nNew > 1e-14) {
        HNnew /= nNew;
      }
      P = Theta;
      HP = HTheta;
      N = Nnew;
      HN = HNnew;
      CN = N.dot(HN);
      curvatureHistory.push_back(CN);
      if (CN < bestCN) {
        bestCN = CN;
        bestN = N;
      }
      F = HN - CN * N;
      applyMask(F);
      Fnorm = F.norm();
      statsRotations = k;
      continue;
    }
    P3 /= pNrm;
    // HP3: project H P into P3 direction approximately
    HP3 = HP - N.dot(HP) * N - Theta.dot(HP) * Theta;
    applyMask(HP3);

    Eigen::Matrix3d A3 = Eigen::Matrix3d::Zero();
    const VectorXd basis[3] = {N, Theta, P3};
    const VectorXd Hbasis[3] = {HN, HTheta, HP3};
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
        A3(i, j) = basis[i].dot(Hbasis[j]);
        A3(j, i) = A3(i, j);
      }
    }
    A3 = 0.5 * (A3 + A3.transpose());
    // Orthonormal basis → standard eigenproblem (stable vs generalized B)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ges(A3);
    if (ges.info() != Eigen::Success) {
      QUILL_LOG_WARNING(log,
                        "[LOR] 3x3 eigen failed at iter={}; stop",
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
    // Paper: under quadratic PES + force translation, C is non-increasing.
    // Reject updates that increase C (FD/GEP noise) — restores history monotonicity
    // and keeps the mode in the softest-mode basin (mode agreement with Lanczos/CG).
    if (!curvatureHistory.empty() && CNnew > curvatureHistory.back() + 1e-4) {
      QUILL_LOG_INFO(log,
                     "[LOR] curvature stall iter={} C_prev={:.6f} C_new={:.6f} "
                     "(reject update; keep prior mode)",
                     k, curvatureHistory.back(), CNnew);
      // Revert to pre-update state: restore from best tracked softest mode
      N = bestN;
      CN = bestCN;
      HN = hessianVector(VectorXd(), x0_r, N, freeMask, delta);
      applyMask(HN);
      CN = N.dot(HN);
      F = HN - CN * N;
      applyMask(F);
      Fnorm = F.norm();
      statsRotations = k;
      break;
    }
    curvatureHistory.push_back(CNnew);
    CN = CNnew;
    if (CN < bestCN) {
      bestCN = CN;
      bestN = N;
    }

    F = HN - CN * N;
    applyMask(F);
    Fnorm = F.norm();
    statsRotations = k;

    QUILL_LOG_INFO(log,
                   "[LOR] iter={} (3x3) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} "
                   "b={:.4f} c={:.4f}",
                   k, Fnorm, CN, a, b, c);

    if (relativeResidual(Fnorm, CN) < residualTol) {
      convergedOnResidual = true;
      break;
    }
  }

  // Optional FD refresh only if it does not worsen the softest curvature.
  {
    VectorXd HNfd = hessianVector(VectorXd(), x0_r, bestN, freeMask, delta);
    applyMask(HNfd);
    const double Cfd = bestN.dot(HNfd);
    if (Cfd <= bestCN + 1e-4) {
      bestCN = Cfd;
      HN = HNfd;
      N = bestN;
      CN = Cfd;
      if (curvatureHistory.empty() ||
          Cfd <= curvatureHistory.back() + 1e-4) {
        curvatureHistory.push_back(Cfd);
      }
    } else {
      N = bestN;
      CN = bestCN;
    }
  }
  eigenvalue = CN;
  eigenvector = N;
  QUILL_LOG_INFO(log,
                 "[LOR] done rotations={} force_calls={} C_N={:.6f} "
                 "converged_residual={} best_tracked",
                 statsRotations, totalForceCalls, eigenvalue,
                 convergedOnResidual);
}
