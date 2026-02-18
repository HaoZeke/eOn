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
#include "fpe_handler.h"

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

  // featomic triggers harmless FE_INVALID during radial basis precomputation
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();
  calculator_ = torch::make_intrusive<featomic_torch::CalculatorHolder>(
      "soap_power_spectrum", params_json);
  fpeh.restore_fpe();

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
  // featomic / torch operations produce harmless FE_INVALID
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  auto pos = positions.to(torch::kFloat64).to(torch::kCPU).detach().clone();
  auto typ = types.to(torch::kInt32).to(torch::kCPU);
  auto cel = cell.to(torch::kFloat64).to(torch::kCPU);
  auto pb = pbc.to(torch::kBool).to(torch::kCPU);

  auto system =
      torch::make_intrusive<metatomic_torch::SystemHolder>(typ, pos, cel, pb);

  // Descriptor-only: no gradients needed
  auto options =
      torch::make_intrusive<featomic_torch::CalculatorOptionsHolder>();

  auto tensor_map = calculator_->compute({system}, options);

  // Extract and average over atoms
  // The TensorMap has blocks keyed by (species_1, species_2) pairs.
  // We concatenate all blocks' values and average over atoms.
  auto keys = tensor_map->keys();
  int64_t n_blocks = keys->count();

  std::vector<torch::Tensor> block_means;

  for (int64_t b = 0; b < n_blocks; b++) {
    auto block =
        metatensor_torch::TensorMapHolder::block_by_id(tensor_map, b);
    auto values = block->values(); // [n_atoms_in_block, n_components, n_props]

    auto n_samples = values.size(0);
    auto flat = values.reshape({n_samples, -1}); // [n_atoms, D_block]

    auto block_sum = flat.sum(0); // [D_block]
    block_means.push_back(block_sum);
  }

  auto descriptor = torch::cat(block_means, 0); // [D_total]
  int64_t n_atoms = positions.size(0);
  descriptor = descriptor / static_cast<double>(n_atoms);

  descriptor_dim_ = static_cast<int>(descriptor.size(0));

  fpeh.restore_fpe();
  return descriptor;
}

torch::Tensor SoapDescriptorEngine::computeJacobian(torch::Tensor positions,
                                                      torch::Tensor types,
                                                      torch::Tensor cell,
                                                      torch::Tensor pbc) {
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  int64_t n_atoms = positions.size(0);
  int64_t n_cart = 3 * n_atoms;

  auto pos = positions.to(torch::kFloat64).to(torch::kCPU).detach().clone();
  auto typ = types.to(torch::kInt32).to(torch::kCPU);
  auto cel = cell.to(torch::kFloat64).to(torch::kCPU);
  auto pb = pbc.to(torch::kBool).to(torch::kCPU);

  auto system =
      torch::make_intrusive<metatomic_torch::SystemHolder>(typ, pos, cel, pb);

  auto options =
      torch::make_intrusive<featomic_torch::CalculatorOptionsHolder>();
  options->gradients = {"positions"};

  auto tensor_map = calculator_->compute({system}, options);

  // Extract the Jacobian directly from featomic's analytic gradient blocks
  // instead of doing D individual backward() calls through autograd.
  auto keys = tensor_map->keys();
  int64_t n_blocks = keys->count();

  std::vector<torch::Tensor> block_descs;
  std::vector<torch::Tensor> block_jacs;

  for (int64_t b = 0; b < n_blocks; b++) {
    auto block =
        metatensor_torch::TensorMapHolder::block_by_id(tensor_map, b);
    auto values = block->values();
    auto n_samples = values.size(0);
    auto flat = values.reshape({n_samples, -1});
    int64_t D_block = flat.size(1);

    block_descs.push_back(flat.sum(0)); // [D_block]

    // Position gradient block: values [n_grad, 3, ...components..., n_props]
    // For power spectrum (no extra components): [n_grad, 3, n_props]
    // Samples columns: ["sample", "structure", "atom"] or ["sample", "atom"]
    auto grad_block =
        metatensor_torch::TensorBlockHolder::gradient(block, "positions");
    auto grad_values = grad_block->values();
    auto grad_samples = grad_block->samples();
    int64_t n_grad = grad_values.size(0);

    auto grad_flat = grad_values.reshape({n_grad, 3, -1}); // [n_grad, 3, D_block]

    // Find the "atom" column in gradient samples (always last)
    int atom_col = static_cast<int>(grad_samples->names().size()) - 1;
    auto samples_tensor = grad_samples->values(); // [n_grad, n_cols] int32

    // Accumulate: J_block[d, 3*atom + xyz] = sum over grad samples with that atom
    auto J_block = torch::zeros({D_block, n_cart}, torch::kFloat64);

    for (int64_t g = 0; g < n_grad; g++) {
      int64_t atom_idx = samples_tensor[g][atom_col].item<int32_t>();
      for (int xyz = 0; xyz < 3; xyz++) {
        J_block.index({torch::indexing::Slice(), 3 * atom_idx + xyz}) +=
            grad_flat.index({g, xyz, torch::indexing::Slice()});
      }
    }

    block_jacs.push_back(J_block);
  }

  // Concatenate blocks and normalize by atom count
  auto jacobian = torch::cat(block_jacs, 0) / static_cast<double>(n_atoms);
  last_descriptor_ =
      torch::cat(block_descs, 0) / static_cast<double>(n_atoms);
  descriptor_dim_ = static_cast<int>(last_descriptor_.size(0));

  fpeh.restore_fpe();
  return jacobian;
}

#endif // WITH_FEATOMIC
