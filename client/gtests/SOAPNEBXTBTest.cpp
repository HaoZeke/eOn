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

/// Integration test for SOAP-NEB with XTB potential on a 9-atom molecule.
///
/// Verifies that the SOAP descriptor-space NEB produces a converged band
/// with at least one extremum, matching the CI-NEB XTB test but operating
/// in SOAP space rather than Cartesian space.

#include "SoapNudgedElasticBand.h"
#include "catch2/catch_amalgamated.hpp"
#include <spdlog/sinks/null_sink.h>
#include <spdlog/spdlog.h>

namespace tests {

struct SOAPNEBLoggerSetup {
  SOAPNEBLoggerSetup() {
    if (!spdlog::get("combi")) {
      auto null_sink = std::make_shared<spdlog::sinks::null_sink_mt>();
      auto logger = std::make_shared<spdlog::logger>("combi", null_sink);
      spdlog::register_logger(logger);
    }
  }
};
static SOAPNEBLoggerSetup _soapneb_logger_setup;

TEST_CASE("SOAP-NEB XTB integration", "[neb][xtb][soap]") {
  Parameters params;
  params.potential_options.potential = PotType::XTB;
  params.xtb_options.paramset = "GFN2xTB";
  params.xtb_options.acc = 1.0;
  params.xtb_options.elec_temperature = 300.0;
  params.xtb_options.maxiter = 250;

  params.neb_options.image_count = 10;
  params.neb_options.spring.weighting.enabled = true;
  params.neb_options.spring.weighting.k_min = 0.972;
  params.neb_options.spring.weighting.k_max = 9.72;
  params.neb_options.spring.weighting.trigger = 0.5;
  params.neb_options.initialization.method = NEBInit::LINEAR;
  params.neb_options.endpoints.minimize = false;
  params.neb_options.climbing_image.enabled = true;
  params.neb_options.climbing_image.converged_only = true;
  params.neb_options.climbing_image.trigger_force = 0.5;
  params.neb_options.climbing_image.trigger_factor = 0.8;
  params.neb_options.force_tolerance = 0.0514221;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.max_iterations = 100;
  params.optimizer_options.max_move = 0.1;

  // Enable SOAP-space NEB
  params.soap_neb_options.enabled = true;
  params.soap_neb_options.cutoff_radius = 6.0;
  params.soap_neb_options.smoothing_width = 0.5;
  params.soap_neb_options.density_width = 0.3;
  params.soap_neb_options.max_angular = 4;
  params.soap_neb_options.max_radial = 6;

  auto pot = helper_functions::makePotential(params.potential_options.potential,
                                             params);
  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  std::string reactFile("reactant.con");
  std::string prodFile("product.con");
  initial->con2matter(reactFile);
  final_state->con2matter(prodFile);

  auto neb = std::make_unique<SoapNudgedElasticBand>(initial, final_state,
                                                      params, pot);
  auto status = neb->compute();

  REQUIRE(static_cast<int>(status) ==
          static_cast<int>(NudgedElasticBand::NEBStatus::GOOD));

  neb->findExtrema();
  REQUIRE(neb->numExtrema >= 1);
}

TEST_CASE("SOAP-space optimizer XTB integration", "[neb][xtb][soap][soapopt]") {
  Parameters params;
  params.potential_options.potential = PotType::XTB;
  params.xtb_options.paramset = "GFN2xTB";
  params.xtb_options.acc = 1.0;
  params.xtb_options.elec_temperature = 300.0;
  params.xtb_options.maxiter = 250;

  params.neb_options.image_count = 10;
  params.neb_options.spring.weighting.enabled = false;
  params.neb_options.spring.constant = 1.0;
  params.neb_options.initialization.method = NEBInit::LINEAR;
  params.neb_options.endpoints.minimize = false;
  params.neb_options.climbing_image.enabled = true;
  params.neb_options.climbing_image.converged_only = true;
  params.neb_options.climbing_image.trigger_force = 0.5;
  params.neb_options.climbing_image.trigger_factor = 0.8;
  params.neb_options.force_tolerance = 0.0514221;
  params.neb_options.max_iterations = 500;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.max_iterations = 500;
  params.optimizer_options.max_move = 0.05;
  params.optimizer_options.lbfgs.auto_scale = false;
  params.optimizer_options.lbfgs.inverse_curvature = 0.01;

  // Enable SOAP-space NEB with SOAP-space optimizer
  params.soap_neb_options.enabled = true;
  params.soap_neb_options.soap_space_optimizer = true;
  params.soap_neb_options.cutoff_radius = 6.0;
  params.soap_neb_options.smoothing_width = 0.5;
  params.soap_neb_options.density_width = 0.3;
  params.soap_neb_options.max_angular = 4;
  params.soap_neb_options.max_radial = 6;

  auto pot = helper_functions::makePotential(params.potential_options.potential,
                                             params);
  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  std::string reactFile("reactant.con");
  std::string prodFile("product.con");
  initial->con2matter(reactFile);
  final_state->con2matter(prodFile);

  auto neb = std::make_unique<SoapNudgedElasticBand>(initial, final_state,
                                                      params, pot);

  auto status = neb->compute();

  REQUIRE(static_cast<int>(status) ==
          static_cast<int>(NudgedElasticBand::NEBStatus::GOOD));

  neb->findExtrema();
  REQUIRE(neb->numExtrema >= 1);
}

} /* namespace tests */
