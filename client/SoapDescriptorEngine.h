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

#include <torch/torch.h>

#include <featomic/torch.hpp>
#include <metatomic/torch.hpp>
#include <metatensor/torch.hpp>

#include <string>
#include <vector>

/// Computes averaged SOAP power spectrum descriptors and the Jacobian
/// dS_avg/dR for the SOAP-NEB force transformations.
///
/// The descriptor is the average over all atoms of the per-atom SOAP power
/// spectrum vector, yielding a single D-dimensional representation of the
/// entire configuration. The Jacobian [D, 3N] maps Cartesian perturbations
/// to descriptor perturbations.
class SoapDescriptorEngine {
public:
  SoapDescriptorEngine(double cutoff, double smoothing_width,
                       double density_width, int max_angular, int max_radial);

  /// Compute the averaged SOAP descriptor for a configuration.
  /// Returns: D-dimensional 1D tensor.
  torch::Tensor compute(torch::Tensor positions, // [N,3] f64
                         torch::Tensor types,     // [N] i32
                         torch::Tensor cell,      // [3,3] f64
                         torch::Tensor pbc);      // [3] bool

  /// Compute the full Jacobian J = dS_avg/dR, shape [D, 3N].
  /// Extracts analytic gradients from featomic's gradient blocks.
  /// Also stores the descriptor, accessible via lastDescriptor().
  torch::Tensor computeJacobian(torch::Tensor positions, // [N,3] f64
                                 torch::Tensor types,     // [N] i32
                                 torch::Tensor cell,      // [3,3] f64
                                 torch::Tensor pbc);      // [3] bool

  /// Dimension of the averaged SOAP descriptor.
  /// Only valid after at least one compute() call.
  int descriptorDim() const { return descriptor_dim_; }

  /// Access the descriptor computed during the last computeJacobian() call.
  const torch::Tensor &lastDescriptor() const { return last_descriptor_; }

private:
  torch::intrusive_ptr<featomic_torch::CalculatorHolder> calculator_;
  double cutoff_;
  double smoothing_width_;
  double density_width_;
  int max_angular_;
  int max_radial_;
  int descriptor_dim_{0};

  /// Descriptor from the last computeJacobian() call.
  torch::Tensor last_descriptor_;
};

#endif // WITH_FEATOMIC
