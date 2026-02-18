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
#ifdef WITH_FEATOMIC

#include "SoapNudgedElasticBand.h"
#include "fpe_handler.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>

// --- Constructors ---

SoapNudgedElasticBand::SoapNudgedElasticBand(
    std::shared_ptr<Matter> initialPassed, std::shared_ptr<Matter> finalPassed,
    const Parameters &parametersPassed, std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(initialPassed, finalPassed, parametersPassed, potPassed),
      soap_engine_(parametersPassed.soap_neb_options.cutoff_radius,
                   parametersPassed.soap_neb_options.smoothing_width,
                   parametersPassed.soap_neb_options.density_width,
                   parametersPassed.soap_neb_options.max_angular,
                   parametersPassed.soap_neb_options.max_radial) {
  soap_log_ = spdlog::get("combi");
  SPDLOG_LOGGER_INFO(soap_log_, "[SoapNEB] Initialized SOAP-space NEB with "
                                 "cutoff={}, l_max={}, n_max={}, "
                                 "soap_space_optimizer={}",
                     parametersPassed.soap_neb_options.cutoff_radius,
                     parametersPassed.soap_neb_options.max_angular,
                     parametersPassed.soap_neb_options.max_radial,
                     parametersPassed.soap_neb_options.soap_space_optimizer);
}

SoapNudgedElasticBand::SoapNudgedElasticBand(
    std::vector<Matter> initPath, const Parameters &parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(std::move(initPath), parametersPassed, potPassed),
      soap_engine_(parametersPassed.soap_neb_options.cutoff_radius,
                   parametersPassed.soap_neb_options.smoothing_width,
                   parametersPassed.soap_neb_options.density_width,
                   parametersPassed.soap_neb_options.max_angular,
                   parametersPassed.soap_neb_options.max_radial) {
  soap_log_ = spdlog::get("combi");
}

// --- Factory method ---

std::shared_ptr<ObjectiveFunction>
SoapNudgedElasticBand::createObjectiveFunction() {
  if (params.soap_neb_options.soap_space_optimizer) {
    SPDLOG_LOGGER_INFO(soap_log_,
                       "[SoapNEB] Using SOAP-space objective function");
    return std::make_shared<SoapSpaceNEBObjectiveFunction>(this, params);
  }
  // Default: Cartesian optimizer with SOAP-projected forces
  return NudgedElasticBand::createObjectiveFunction();
}

// --- Helper: Matter -> torch tensors ---

void SoapNudgedElasticBand::matterToTensors(const Matter &m,
                                             torch::Tensor &positions,
                                             torch::Tensor &types,
                                             torch::Tensor &cell,
                                             torch::Tensor &pbc) const {
  int nAtoms = m.numberOfAtoms();

  // Positions: [N, 3] f64
  const AtomMatrix &pos = m.getPositions();
  positions = torch::from_blob(const_cast<double *>(pos.data()), {nAtoms, 3},
                               torch::kFloat64)
                  .clone();

  // Atomic types: [N] i32
  VectorXi atomicNrs = m.getAtomicNrs();
  std::vector<int32_t> types_vec(nAtoms);
  for (int i = 0; i < nAtoms; i++) {
    types_vec[i] = static_cast<int32_t>(atomicNrs[i]);
  }
  types = torch::tensor(types_vec, torch::kInt32);

  // Cell: [3, 3] f64
  Matrix3d cellMat = m.getCell();
  cell = torch::from_blob(const_cast<double *>(cellMat.data()), {3, 3},
                           torch::kFloat64)
             .clone();

  // PBC: check if cell norms are non-negligible
  auto cell_norms = torch::norm(cell, 2, /*dim=*/1);
  pbc = cell_norms.abs() > 1e-9;
}

// --- Recompute SOAP data for all images ---

void SoapNudgedElasticBand::recomputeSoapData() {
  soap_descriptors_.resize(numImages + 2);
  soap_jacobians_.resize(numImages + 2);

  for (long i = 0; i <= numImages + 1; i++) {
    torch::Tensor positions, types, cell, pbc;
    matterToTensors(*path[i], positions, types, cell, pbc);

    if (i >= 1 && i <= numImages) {
      // Intermediate images: computeJacobian also produces the descriptor
      soap_jacobians_[i] =
          soap_engine_.computeJacobian(positions, types, cell, pbc);
      soap_descriptors_[i] = soap_engine_.lastDescriptor();
    } else {
      // Endpoints: only need descriptor
      soap_descriptors_[i] = soap_engine_.compute(positions, types, cell, pbc);
    }
  }

  // Cache descriptor dimension
  if (descriptor_dim_ == 0 && soap_descriptors_[0].defined()) {
    descriptor_dim_ = static_cast<int>(soap_descriptors_[0].size(0));
  }
}

// --- Recompute descriptors only (no Jacobians) for setPositions consistency ---

void SoapNudgedElasticBand::recomputeDescriptorsOnly() {
  soap_descriptors_.resize(numImages + 2);

  for (long i = 0; i <= numImages + 1; i++) {
    torch::Tensor positions, types, cell, pbc;
    matterToTensors(*path[i], positions, types, cell, pbc);
    soap_descriptors_[i] = soap_engine_.compute(positions, types, cell, pbc);
  }
}

// --- SOAP-space NEB force update ---

void SoapNudgedElasticBand::updateForces() {
  // featomic / torch / metatomic trigger harmless FE_INVALID
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  // 1. Update Cartesian forces for all intermediate images
  for (long i = 1; i <= numImages; i++) {
    path[i]->getForces();
  }

  // 2. Recompute SOAP descriptors and Jacobians
  recomputeSoapData();

  // 3. Find highest-energy non-endpoint image
  auto first = path.begin() + 1;
  auto last = path.begin() + numImages + 1;
  auto it = std::max_element(
      first, last,
      [](const std::shared_ptr<Matter> &a, const std::shared_ptr<Matter> &b) {
        return a->getPotentialEnergy() < b->getPotentialEnergy();
      });
  maxEnergyImage = std::distance(path.begin(), it);

  // Reset climbing image if CI is disabled
  if (!params.neb_options.climbing_image.enabled) {
    climbingImage = 0;
  }

  // Spring constant
  double k = params.neb_options.spring.constant;

  // Resize SOAP NEB forces storage if in SOAP-space optimizer mode
  if (params.soap_neb_options.soap_space_optimizer) {
    soap_neb_forces_.resize(numImages + 2);
  }

  // 4. Projection loop: compute NEB forces in SOAP space
  for (long i = 1; i <= numImages; i++) {
    // SOAP descriptors for this image and neighbors
    auto &S_i = soap_descriptors_[i];
    auto &S_prev = soap_descriptors_[i - 1];
    auto &S_next = soap_descriptors_[i + 1];
    auto &J_i = soap_jacobians_[i]; // [D, 3N]

    // Energies
    double energy = path[i]->getPotentialEnergy();
    double energyPrev = path[i - 1]->getPotentialEnergy();
    double energyNext = path[i + 1]->getPotentialEnergy();

    // SOAP-space distances
    auto diff_next = S_next - S_i;
    auto diff_prev = S_i - S_prev;
    double distNext = diff_next.norm().item<double>();
    double distPrev = diff_prev.norm().item<double>();

    // SOAP-space tangent (improved tangent method, energy-weighted)
    torch::Tensor tau_S;
    if (energyNext > energy && energy > energyPrev) {
      // Ascending: tangent toward next image
      tau_S = diff_next;
    } else if (energyPrev > energy && energy > energyNext) {
      // Descending: tangent toward previous image (from prev to current)
      tau_S = diff_prev;
    } else {
      // At an extremum: energy-weighted bisection
      double energyDiffPrev = energyPrev - energy;
      double energyDiffNext = energyNext - energy;
      double minDE = std::min(std::abs(energyDiffPrev), std::abs(energyDiffNext));
      double maxDE = std::max(std::abs(energyDiffPrev), std::abs(energyDiffNext));

      if (energyDiffPrev > energyDiffNext) {
        tau_S = diff_next * minDE + diff_prev * maxDE;
      } else {
        tau_S = diff_next * maxDE + diff_prev * minDE;
      }
    }

    // Normalize tangent
    double tauNorm = tau_S.norm().item<double>();
    if (tauNorm > 1e-10) {
      tau_S = tau_S / tauNorm;
    } else {
      tau_S = diff_next;
      tauNorm = tau_S.norm().item<double>();
      if (tauNorm > 1e-10) {
        tau_S = tau_S / tauNorm;
      }
    }

    // Project Cartesian force into SOAP space: F_S = J @ f_cart
    AtomMatrix forceCart = path[i]->getForces(); // [N, 3]
    auto f_cart_flat =
        torch::from_blob(const_cast<double *>(forceCart.data()),
                         {3 * atoms}, torch::kFloat64)
            .clone();
    auto F_S = torch::mv(J_i, f_cart_flat); // D-vector

    // Force components in SOAP space
    double F_par_mag = torch::dot(F_S, tau_S).item<double>();
    auto F_par_S = F_par_mag * tau_S;
    auto F_perp_S = F_S - F_par_S;

    // Spring force in SOAP space
    auto F_spring_S = k * (distNext - distPrev) * tau_S;

    // NEB force in SOAP space
    torch::Tensor F_NEB_S;
    if (params.neb_options.climbing_image.enabled &&
        static_cast<long>(i) == static_cast<long>(maxEnergyImage)) {
      climbingImage = maxEnergyImage;
      // Climbing image: invert tangential component
      F_NEB_S = F_S - 2.0 * F_par_mag * tau_S;
    } else {
      // Standard NEB: perpendicular potential + parallel spring
      F_NEB_S = F_perp_S + F_spring_S;
    }

    // Store SOAP-space NEB forces for the SOAP-space optimizer
    if (params.soap_neb_options.soap_space_optimizer) {
      soap_neb_forces_[i] = F_NEB_S.clone();
    }

    // Map back to Cartesian via Jacobian transpose: F_cart_NEB = J^T @ F_NEB_S
    auto F_NEB_cart = torch::mv(J_i.t(), F_NEB_S); // [3N]

    // Copy to projectedForce
    auto F_ptr = F_NEB_cart.contiguous().data_ptr<double>();
    for (int a = 0; a < atoms; a++) {
      (*projectedForce[i])(a, 0) = F_ptr[3 * a + 0];
      (*projectedForce[i])(a, 1) = F_ptr[3 * a + 1];
      (*projectedForce[i])(a, 2) = F_ptr[3 * a + 2];
    }

    // Store the SOAP tangent in the Cartesian tangent slot for diagnostics:
    // map tau_S back to Cartesian via J^T @ tau_S
    auto tau_cart = torch::mv(J_i.t(), tau_S);
    auto tau_cart_norm = tau_cart.norm().item<double>();
    if (tau_cart_norm > 1e-10) {
      tau_cart = tau_cart / tau_cart_norm;
    }
    auto tau_ptr = tau_cart.contiguous().data_ptr<double>();
    for (int a = 0; a < atoms; a++) {
      (*tangent[i])(a, 0) = tau_ptr[3 * a + 0];
      (*tangent[i])(a, 1) = tau_ptr[3 * a + 1];
      (*tangent[i])(a, 2) = tau_ptr[3 * a + 2];
    }

    // Zero net translational force (for free-atom systems)
    if (path[i]->numberOfFreeAtoms() == path[i]->numberOfAtoms()) {
      for (int j = 0; j <= 2; j++) {
        double translationMag = projectedForce[i]->col(j).sum();
        int natoms_col = projectedForce[i]->col(j).size();
        projectedForce[i]->col(j).array() -=
            translationMag / static_cast<double>(natoms_col);
      }
    }
  }

  // Flag that forces are fresh
  movedAfterForceCall = false;

  fpeh.restore_fpe();
}

// ==========================================================================
// SoapSpaceNEBObjectiveFunction implementation
// ==========================================================================

VectorXd SoapSpaceNEBObjectiveFunction::getGradient(bool /*fdstep*/) {
  if (soap_neb_->movedAfterForceCall) {
    soap_neb_->updateForces();
  }

  int D = soap_neb_->descriptor_dim_;
  VectorXd grad(D * soap_neb_->numImages);

  for (long i = 1; i <= soap_neb_->numImages; i++) {
    auto &F_S = soap_neb_->soap_neb_forces_[i];
    auto F_ptr = F_S.contiguous().data_ptr<double>();
    // Gradient = -force (optimizer minimizes)
    for (int d = 0; d < D; d++) {
      grad[(i - 1) * D + d] = -F_ptr[d];
    }
  }

  return grad;
}

double SoapSpaceNEBObjectiveFunction::getEnergy() {
  double energy = 0;
  for (long i = 1; i <= soap_neb_->numImages; i++) {
    energy += soap_neb_->path[i]->getPotentialEnergy();
  }
  return energy;
}

void SoapSpaceNEBObjectiveFunction::setPositions(VectorXd x) {
  // featomic / torch trigger harmless FE_INVALID
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  soap_neb_->movedAfterForceCall = true;

  int D = soap_neb_->descriptor_dim_;
  int nCart = 3 * soap_neb_->atoms;

  for (long i = 1; i <= soap_neb_->numImages; i++) {
    // Extract target SOAP descriptor for this image
    auto S_new_eigen = x.segment((i - 1) * D, D);
    auto S_new = torch::from_blob(const_cast<double *>(S_new_eigen.data()),
                                   {D}, torch::kFloat64)
                     .clone();

    // Current descriptor and Jacobian
    auto &S_cur = soap_neb_->soap_descriptors_[i];
    auto &J_i = soap_neb_->soap_jacobians_[i]; // [D, 3N]

    // Delta in SOAP space
    auto dS = S_new - S_cur; // [D]

    // Compute ΔR = J^+ · ΔS via regularized normal equations.
    // J is [D, 3N] with D >> 3N, so J^T J is small [3N, 3N].
    // Tikhonov regularization prevents amplification of near-null-space
    // directions (translations, rotations, permutational equivalences).
    // ΔR = (J^T J + λI)^{-1} J^T ΔS
    auto JtJ = torch::mm(J_i.t(), J_i);             // [3N, 3N]
    // λ controls regularization strength: larger λ → more conservative
    // Cartesian steps, effectively a trust region in Cartesian space.
    auto reg = 1e-2 * torch::eye(nCart, torch::kFloat64);
    auto Jt_dS = torch::mv(J_i.t(), dS);            // [3N]
    auto dR = torch::mv(torch::inverse(JtJ + reg), Jt_dS); // [3N]

    // Cap the Cartesian displacement to max_move per atom.
    // The LBFGS step capping operates on SOAP components which don't
    // directly constrain Cartesian motion; this safeguard prevents
    // atoms from moving into unphysical configurations.
    double maxMove = params.optimizer_options.max_move;
    auto dR_cont = dR.contiguous();
    auto dR_ptr = dR_cont.data_ptr<double>();
    double maxAtomDisp = 0.0;
    for (int a = 0; a < soap_neb_->atoms; a++) {
      double disp = std::sqrt(dR_ptr[3 * a + 0] * dR_ptr[3 * a + 0] +
                              dR_ptr[3 * a + 1] * dR_ptr[3 * a + 1] +
                              dR_ptr[3 * a + 2] * dR_ptr[3 * a + 2]);
      if (disp > maxAtomDisp)
        maxAtomDisp = disp;
    }
    double scale = 1.0;
    if (maxAtomDisp > maxMove) {
      scale = maxMove / maxAtomDisp;
    }

    // Update Cartesian positions
    const AtomMatrix &curPos = soap_neb_->path[i]->getPositions();
    AtomMatrix newPos(soap_neb_->atoms, 3);
    for (int a = 0; a < soap_neb_->atoms; a++) {
      newPos(a, 0) = curPos(a, 0) + scale * dR_ptr[3 * a + 0];
      newPos(a, 1) = curPos(a, 1) + scale * dR_ptr[3 * a + 1];
      newPos(a, 2) = curPos(a, 2) + scale * dR_ptr[3 * a + 2];
    }
    soap_neb_->path[i]->setPositions(newPos);
  }

  // Recompute SOAP descriptors (without Jacobians) for getPositions consistency
  soap_neb_->recomputeDescriptorsOnly();

  fpeh.restore_fpe();
}

VectorXd SoapSpaceNEBObjectiveFunction::getPositions() {
  int D = soap_neb_->descriptor_dim_;
  VectorXd pos(D * soap_neb_->numImages);

  for (long i = 1; i <= soap_neb_->numImages; i++) {
    auto &S_i = soap_neb_->soap_descriptors_[i];
    auto S_ptr = S_i.contiguous().data_ptr<double>();
    for (int d = 0; d < D; d++) {
      pos[(i - 1) * D + d] = S_ptr[d];
    }
  }

  return pos;
}

int SoapSpaceNEBObjectiveFunction::degreesOfFreedom() {
  return soap_neb_->descriptor_dim_ * soap_neb_->numImages;
}

bool SoapSpaceNEBObjectiveFunction::isConverged() {
  return getConvergence() < params.neb_options.force_tolerance;
}

double SoapSpaceNEBObjectiveFunction::getConvergence() {
  return soap_neb_->convergenceForce();
}

VectorXd SoapSpaceNEBObjectiveFunction::difference(VectorXd a, VectorXd b) {
  // No PBC in SOAP descriptor space
  return a - b;
}

#endif // WITH_FEATOMIC
