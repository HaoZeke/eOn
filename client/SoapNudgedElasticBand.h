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

#ifdef WITH_FEATOMIC

#include "NudgedElasticBand.h"
#include "SoapDescriptorEngine.h"

#include <torch/torch.h>
#include <vector>

/// NEB operating entirely in SOAP descriptor space.
///
/// All path-geometric operations (distances, tangents, spring forces, force
/// decomposition) happen in R^D, where D is the SOAP power spectrum dimension.
/// Cartesian coordinates are maintained at each image only for energy/force
/// evaluation via the metatomic potential.
///
/// At convergence with climbing image, the Cartesian forces are verified to be
/// zero (modulo rigid-body/permutation null space of the Jacobian).
class SoapNudgedElasticBand : public NudgedElasticBand {
public:
  SoapNudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                        std::shared_ptr<Matter> finalPassed,
                        const Parameters &parametersPassed,
                        std::shared_ptr<Potential> potPassed);

  SoapNudgedElasticBand(std::vector<Matter> initPath,
                        const Parameters &parametersPassed,
                        std::shared_ptr<Potential> potPassed);

  ~SoapNudgedElasticBand() = default;

  /// Override: NEB force computation in SOAP space.
  void updateForces(void) override;

private:
  SoapDescriptorEngine soap_engine_;

  /// Per-image SOAP descriptors, S_i (D-vector each)
  std::vector<torch::Tensor> soap_descriptors_;

  /// Per-image Jacobians, J_i ([D, 3N] each)
  std::vector<torch::Tensor> soap_jacobians_;

  /// Recompute all SOAP descriptors and Jacobians from current positions.
  void recomputeSoapData();

  /// Build torch tensors from a Matter object for featomic computation.
  void matterToTensors(const Matter &m, torch::Tensor &positions,
                       torch::Tensor &types, torch::Tensor &cell,
                       torch::Tensor &pbc) const;

  std::shared_ptr<spdlog::logger> soap_log_;
};

/// Objective function wrapper for SOAP-NEB.
///
/// The optimizer still operates on Cartesian positions via
/// getPositions/setPositions, but getGradient() returns the SOAP-projected
/// NEB forces. The difference() method uses Cartesian PBC differences.
class SoapNEBObjectiveFunction : public NEBObjectiveFunction {
public:
  SoapNEBObjectiveFunction(SoapNudgedElasticBand *nebPassed,
                           const Parameters &parametersPassed)
      : NEBObjectiveFunction(nebPassed, parametersPassed),
        soap_neb_{nebPassed} {}

  ~SoapNEBObjectiveFunction() = default;

private:
  SoapNudgedElasticBand *soap_neb_;
};

#endif // WITH_FEATOMIC
