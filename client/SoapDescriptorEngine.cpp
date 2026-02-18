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

#include "SoapDescriptorEngine.h"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

SoapDescriptorEngine::SoapDescriptorEngine(double cutoff,
                                           double smoothing_width,
                                           double density_width,
                                           int max_angular, int max_radial)
    : cutoff_(cutoff), smoothing_width_(smoothing_width),
      density_width_(density_width), max_angular_(max_angular),
      max_radial_(max_radial) {

  // Build the JSON parameter string for featomic's soap_power_spectrum
  std::string params_json = fmt::format(
      R"({{
    "cutoff": {{
        "radius": {},
        "smoothing": {{
            "type": "ShiftedCosine",
            "width": {}
        }}
    }},
    "density": {{
        "type": "Gaussian",
        "width": {}
    }},
    "basis": {{
        "type": "TensorProduct",
        "max_angular": {},
        "radial": {{
            "type": "Gto",
            "max_radial": {}
        }}
    }}
}})",
      cutoff_, smoothing_width_, density_width_, max_angular_, max_radial_);

  calculator_ = torch::make_intrusive<featomic_torch::CalculatorHolder>(
      "soap_power_spectrum", params_json);

  auto log = spdlog::get("combi");
  if (log) {
    log->info("[SoapDescriptorEngine] Initialized SOAP power spectrum: "
              "cutoff={}, smoothing_width={}, density_width={}, "
              "l_max={}, n_max={}",
              cutoff_, smoothing_width_, density_width_, max_angular_,
              max_radial_);
  }
}

torch::Tensor SoapDescriptorEngine::compute(torch::Tensor positions,
                                             torch::Tensor types,
                                             torch::Tensor cell,
                                             torch::Tensor pbc) {
  // Create a fresh leaf tensor for autograd.  We always detach+clone so that
  // register_autograd sees a leaf with requires_grad=true, and we can later
  // read .grad() from this exact tensor in computeJacobian().
  auto pos =
      positions.to(torch::kFloat64).to(torch::kCPU).detach().clone();
  pos.requires_grad_(true);
  last_positions_leaf_ = pos;

  auto typ = types.to(torch::kInt32).to(torch::kCPU);
  auto cel = cell.to(torch::kFloat64).to(torch::kCPU);
  auto pb = pbc.to(torch::kBool).to(torch::kCPU);

  // Create metatomic System (same pattern as MetatomicPotential)
  auto system =
      torch::make_intrusive<metatomic_torch::SystemHolder>(typ, pos, cel, pb);

  // Set up calculator options: request position gradients
  auto options =
      torch::make_intrusive<featomic_torch::CalculatorOptionsHolder>();
  options->gradients = {"positions"};

  // Compute SOAP descriptors
  auto tensor_map =
      calculator_->compute({system}, options);

  // Register autograd to connect positions->descriptor graph
  tensor_map = featomic_torch::register_autograd(
      {system}, tensor_map, {"positions"});

  // Extract and average over atoms
  // The TensorMap has blocks keyed by (species_1, species_2) pairs.
  // We concatenate all blocks' values and average over atoms.
  auto keys = tensor_map->keys();
  int64_t n_blocks = keys->count();

  // Collect per-atom descriptors from all blocks
  std::vector<torch::Tensor> block_means;
  int64_t total_props = 0;

  for (int64_t b = 0; b < n_blocks; b++) {
    auto block =
        metatensor_torch::TensorMapHolder::block_by_id(tensor_map, b);
    auto values = block->values(); // [n_atoms_in_block, n_components, n_props]

    // Reshape to [n_atoms, -1] by flattening component and property dimensions
    auto n_samples = values.size(0);
    auto flat = values.reshape({n_samples, -1}); // [n_atoms, D_block]

    // Sum over atoms in this block (some blocks may have different atoms)
    auto block_sum = flat.sum(0); // [D_block]
    block_means.push_back(block_sum);
    total_props += flat.size(1);
  }

  // Concatenate all block contributions and divide by number of atoms
  auto descriptor = torch::cat(block_means, 0); // [D_total]
  int64_t n_atoms = positions.size(0);
  descriptor = descriptor / static_cast<double>(n_atoms);

  descriptor_dim_ = static_cast<int>(descriptor.size(0));

  return descriptor;
}

torch::Tensor SoapDescriptorEngine::computeJacobian(torch::Tensor positions,
                                                      torch::Tensor types,
                                                      torch::Tensor cell,
                                                      torch::Tensor pbc) {
  int64_t n_atoms = positions.size(0);
  int64_t n_cart = 3 * n_atoms;

  // compute() creates a fresh leaf in last_positions_leaf_ and builds
  // the autograd graph through register_autograd.
  auto descriptor = compute(positions, types, cell, pbc);
  int64_t D = descriptor.size(0);

  // Collect gradients from the leaf tensor that compute() actually used.
  auto &pos = last_positions_leaf_;

  // Build Jacobian row by row: J[d, :] = d(S_d) / d(R_flat)
  auto jacobian = torch::zeros({D, n_cart}, torch::kFloat64);

  for (int64_t d = 0; d < D; d++) {
    // Zero gradients from previous iteration
    if (pos.grad().defined()) {
      pos.grad().zero_();
    }

    // Backprop through the d-th descriptor component
    descriptor[d].backward({}, /*retain_graph=*/true);

    if (pos.grad().defined()) {
      // grad is [N, 3], flatten to [3N]
      jacobian[d] = pos.grad().reshape({-1}).clone();
    }
  }

  return jacobian;
}

#endif // WITH_FEATOMIC
