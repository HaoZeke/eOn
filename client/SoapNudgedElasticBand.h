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

class SoapSpaceNEBObjectiveFunction; // forward declaration

/// NEB operating entirely in SOAP descriptor space.
///
/// All path-geometric operations (distances, tangents, spring forces, force
/// decomposition) happen in R^D, where D is the SOAP power spectrum dimension.
/// Cartesian coordinates are maintained at each image only for energy/force
/// evaluation via the metatomic potential.
///
/// Two modes of operation controlled by soap_space_optimizer:
///  - false (default): Cartesian optimizer, SOAP-projected NEB forces
///  - true: Optimizer works in SOAP descriptor space via pseudoinverse
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

protected:
  /// Override: create SOAP-space or Cartesian objective function.
  std::shared_ptr<ObjectiveFunction> createObjectiveFunction() override;

private:
  friend class SoapSpaceNEBObjectiveFunction;

  SoapDescriptorEngine soap_engine_;

  /// Per-image SOAP descriptors, S_i (D-vector each)
  std::vector<torch::Tensor> soap_descriptors_;

  /// Per-image Jacobians, J_i ([D, 3N] each)
  std::vector<torch::Tensor> soap_jacobians_;

  /// Per-image SOAP-space NEB forces, F_NEB_S_i (D-vector each).
  /// Only populated when soap_space_optimizer is true.
  std::vector<torch::Tensor> soap_neb_forces_;

  /// Cached descriptor dimension (set on first updateForces call).
  int descriptor_dim_{0};

  /// Recompute all SOAP descriptors and Jacobians from current positions.
  void recomputeSoapData();

  /// Recompute SOAP descriptors only (no Jacobians) for specified images.
  void recomputeDescriptorsOnly();

  /// Build torch tensors from a Matter object for featomic computation.
  void matterToTensors(const Matter &m, torch::Tensor &positions,
                       torch::Tensor &types, torch::Tensor &cell,
                       torch::Tensor &pbc) const;

  std::shared_ptr<spdlog::logger> soap_log_;
};

/// Objective function that operates entirely in SOAP descriptor space.
///
/// The optimizer sees concatenated SOAP descriptors as "positions" and
/// SOAP-space NEB forces as "gradients". Cartesian updates happen via
/// pseudoinverse mapping: ΔR = J^+ · ΔS.
class SoapSpaceNEBObjectiveFunction : public ObjectiveFunction {
public:
  SoapSpaceNEBObjectiveFunction(SoapNudgedElasticBand *nebPassed,
                                const Parameters &parametersPassed)
      : ObjectiveFunction(nullptr, parametersPassed),
        soap_neb_{nebPassed} {}

  ~SoapSpaceNEBObjectiveFunction() = default;

  VectorXd getGradient(bool fdstep = false) override;
  double getEnergy() override;
  void setPositions(VectorXd x) override;
  VectorXd getPositions() override;
  int degreesOfFreedom() override;
  bool isConverged() override;
  double getConvergence() override;
  VectorXd difference(VectorXd a, VectorXd b) override;

private:
  SoapNudgedElasticBand *soap_neb_;
};

#endif // WITH_FEATOMIC
