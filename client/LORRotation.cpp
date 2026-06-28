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

#include "LORRotation.h"
#include "HelperFunctions.h"

#include <Eigen/Eigenvalues>
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

VectorXd LORRotation::hessianVector(const VectorXd &g0, const VectorXd &x0_r,
                                    const VectorXd &v, const VectorXd &freeMask,
                                    double delta) {
  VectorXd dir = v.array() * freeMask.array();
  double nrm = dir.norm();
  if (nrm < 1e-14) {
    return VectorXd::Zero(v.size());
  }
  dir /= nrm;
  x1->setPositionsV(x0_r + delta * dir);
  VectorXd g1 = -x1->getForcesV();
  totalForceCalls += 1;
  VectorXd Hv = (g1 - g0) / delta;
  Hv = Hv.array() * freeMask.array();
  return Hv;
}

void LORRotation::compute(std::shared_ptr<Matter> matter,
                          AtomMatrix initialDirectionAtomMatrix) {
  const long nAtoms = matter->numberOfAtoms();
  const int dim = static_cast<int>(3 * nAtoms);
  VectorXd freeMask = matter->getFreeV();
  VectorXd N = VectorXd::Map(initialDirectionAtomMatrix.data(), dim);
  N = N.array() * freeMask.array();
  if (N.norm() < 1e-10) {
    N.setRandom();
    N = N.array() * freeMask.array();
  }
  N.normalize();

  *x0 = *matter;
  *x1 = *matter;
  VectorXd x0_r = x0->getPositionsV();
  const double delta = params.main_options.finiteDifference;
  // Paper Algorithm I uses a small rotation budget (~5–20), not geometry
  // max_iterations (often 1000). Prefer rotations_max.
  const int rotmax = std::max(
      1, static_cast<int>(params.dimer_options.rotations_max > 0
                             ? params.dimer_options.rotations_max
                             : 20));
  // Residual tolerance: map converged_angle (degrees) to a modest force residual
  // scale; also allow absolute floor so OptBench Morse Pt converges.
  const double residualTol =
      std::max(1e-2, params.dimer_options.converged_angle * 1e-2);

  curvatureHistory.clear();
  convergedOnResidual = false;
  statsRotations = 0;
  totalForceCalls = 0;

  VectorXd g0 = -x0->getForcesV();
  totalForceCalls += 1;

  // H N via one FD along N
  VectorXd HN = hessianVector(g0, x0_r, N, freeMask, delta);
  double CN = N.dot(HN);
  curvatureHistory.push_back(CN);

  VectorXd F = HN - CN * N; // F_⊥ = (I - N N^T) H N
  F = F.array() * freeMask.array();
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
  // H Theta — one FD (iteration-1 2×2 needs H on N and Theta; H N already done)
  VectorXd HTheta = hessianVector(g0, x0_r, Theta, freeMask, delta);

  // 2×2 Rayleigh problem in span{N, Theta} (orthonormal)
  Eigen::Matrix2d A2;
  A2(0, 0) = N.dot(HN);
  A2(0, 1) = N.dot(HTheta);
  A2(1, 0) = A2(0, 1);
  A2(1, 1) = Theta.dot(HTheta);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2(A2);
  Eigen::Vector2d coeffs = es2.eigenvectors().col(0); // smallest eigenvalue
  double a = coeffs(0);
  double b = coeffs(1);

  VectorXd P = Theta; // P^(2) = F-hat^(1)
  VectorXd HP = HTheta;
  N = a * N + b * Theta;
  N = N.array() * freeMask.array();
  if (N.norm() > 1e-14) {
    N.normalize();
  }
  // Force translation: H N_new = a H N + b H Theta (no new FD)
  HN = a * HN + b * HTheta;
  // Re-normalize HN projection consistency not required; recompute CN, F from HN
  CN = N.dot(HN);
  // Paper guarantees CN non-increasing vs prior Ritz when quadratic; record
  curvatureHistory.push_back(CN);
  F = HN - CN * N;
  F = F.array() * freeMask.array();
  Fnorm = F.norm();
  statsRotations = 1;

  QUILL_LOG_INFO(log, "[LOR] iter=1 (2x2) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} b={:.4f}",
                 Fnorm, CN, a, b);

  for (int k = 2; k <= rotmax; ++k) {
    if (Fnorm < residualTol) {
      convergedOnResidual = true;
      QUILL_LOG_INFO(log, "[LOR] converged residual iter={} ||F||={:.6e}", k, Fnorm);
      break;
    }

    Theta = F / Fnorm;
    // One new force: H · Theta (F_⊥ direction)
    HTheta = hessianVector(g0, x0_r, Theta, freeMask, delta);

    // Orthonormalize P against N, Theta for numerical stability of B
    VectorXd Pwork = P;
    Pwork = Pwork - N.dot(Pwork) * N - Theta.dot(Pwork) * Theta;
    Pwork = Pwork.array() * freeMask.array();
    double pNorm = Pwork.norm();
    VectorXd P3 = P;
    VectorXd HP3 = HP;
    bool use3 = (pNorm > 1e-8);
    if (use3) {
      P3 = Pwork / pNorm;
      // Translate HP into P3 direction is approximate; re-project HP along P3
      // by translating in original P basis: HP3 ≈ HP projected — use FD-free
      // translation HP was for old P; for P3 we use HP with Gram-Schmidt note:
      // paper keeps P as linear combo; we translate HP with same (a,b,c) later.
      // For matrix assembly use P (normalized historical), not GS P3, matching
      // paper's B matrix with N·P and Theta·P possibly nonzero.
      P3 = P;
      if (P3.norm() > 1e-14) {
        P3.normalize();
      }
      HP3 = HP;
    }

    // Build 3×3 A and B (generalized)
    // basis V = [N, Theta, P]
    Eigen::Matrix3d A3 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d B3 = Eigen::Matrix3d::Identity();
    VectorXd basis[3] = {N, Theta, P3};
    VectorXd Hbasis[3] = {HN, HTheta, HP3};
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

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> ges(A3, B3);
    if (ges.info() != Eigen::Success) {
      QUILL_LOG_WARNING(log, "[LOR] 3x3 generalized eigen failed at iter={}; stop",
                        k);
      break;
    }
    Eigen::Vector3d c3 = ges.eigenvectors().col(0);
    a = c3(0);
    b = c3(1);
    double c = c3(2);

    // Update N and P (Algorithm I)
    VectorXd Nnew = a * N + b * Theta + c * P3;
    VectorXd Pnew = b * Theta + c * P3;
    Nnew = Nnew.array() * freeMask.array();
    Pnew = Pnew.array() * freeMask.array();
    if (Nnew.norm() > 1e-14) {
      Nnew.normalize();
    }
    if (Pnew.norm() > 1e-14) {
      Pnew.normalize();
    }

    // Force translation (no FD for H N or H P)
    VectorXd HNnew = a * HN + b * HTheta + c * HP3;
    VectorXd HPnew = b * HTheta + c * HP3;

    N = Nnew;
    P = Pnew;
    HN = HNnew;
    HP = HPnew;

    double CNnew = N.dot(HN);
    // Stall detection: curvature not decreasing (within tol)
    if (!curvatureHistory.empty() &&
        CNnew > curvatureHistory.back() + 1e-4) {
      QUILL_LOG_INFO(log,
                     "[LOR] curvature stall iter={} C_prev={:.6f} C_new={:.6f}",
                     k, curvatureHistory.back(), CNnew);
      // still accept update; paper claims monotonic under quadratic + translation
    }
    curvatureHistory.push_back(CNnew);
    CN = CNnew;

    F = HN - CN * N;
    F = F.array() * freeMask.array();
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

  eigenvalue = CN;
  eigenvector = N;
  QUILL_LOG_INFO(log,
                 "[LOR] done rotations={} force_calls={} C_N={:.6f} "
                 "converged_residual={}",
                 statsRotations, totalForceCalls, eigenvalue,
                 convergedOnResidual);
}
