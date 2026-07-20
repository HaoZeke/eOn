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

#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <cmath>
#include <filesystem>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {

bool file_exists(const char *path) { return std::filesystem::exists(path); }

/// Packaging builds may ship with_fortran pots that fail to dlopen when
/// EON_POTENTIALS_PATH is wrong; treat unloaded pot (energy ~ 0) as skip.
void require_loaded_pot(double energy, const char *name) {
  if (!std::isfinite(energy) || std::abs(energy) < 1e-12) {
    SKIP("potential "
         << name
         << " returned no energy (fortran .so missing or EON_POTENTIALS_PATH)");
  }
}

void require_pos_con() {
  if (!file_exists("pos.con")) {
    SKIP("pos.con missing (expected meson workdir systems/si_diamond)");
  }
}

} // namespace

TEST_CASE("SW potential returns finite energy and forces on Si diamond",
          "[pot][sw][si]") {
  require_pos_con();
  Parameters params;
  params.potential_options.potential = PotType::SW_SI;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_sw_si.dat):
  double energy = matter->getPotentialEnergy();
  require_loaded_pot(energy, "si");
  REQUIRE(energy == Catch::Approx(-16.204955).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(0.005232).epsilon(1e-2));

  // Newton's third law
  AtomMatrix forces = matter->getForces();
  Eigen::Vector3d totalF = forces.colwise().sum();
  REQUIRE(totalF.norm() < 1e-6);
}

TEST_CASE("Tersoff energy matches SVN on Si diamond", "[pot][tersoff][si]") {
  require_pos_con();
  Parameters params;
  params.potential_options.potential = PotType::TERSOFF_SI;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_tersoff_si.dat):
  double energy = matter->getPotentialEnergy();
  require_loaded_pot(energy, "si");
  REQUIRE(energy == Catch::Approx(-17.440266).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(1.082719).epsilon(1e-3));
}

TEST_CASE("EDIP energy matches SVN on Si diamond", "[pot][edip][si]") {
  require_pos_con();
  Parameters params;
  params.potential_options.potential = PotType::EDIP;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_edip.dat):
  double energy = matter->getPotentialEnergy();
  require_loaded_pot(energy, "si");
  REQUIRE(energy == Catch::Approx(-18.838135).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(0.730971).epsilon(1e-3));
}

TEST_CASE("Lenosky energy matches SVN on Si diamond", "[pot][lenosky][si]") {
  require_pos_con();
  Parameters params;
  params.potential_options.potential = PotType::LENOSKY_SI;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  // SVN reference (data/reference/point_lenosky_si.dat):
  double energy = matter->getPotentialEnergy();
  require_loaded_pot(energy, "si");
  REQUIRE(energy == Catch::Approx(-17.284558).epsilon(1e-4));

  double maxForce = matter->getForces().rowwise().norm().maxCoeff();
  REQUIRE(maxForce == Catch::Approx(0.456933).epsilon(1e-3));
}

TEST_CASE("SW and Tersoff give different energies on same Si system",
          "[pot][si]") {
  require_pos_con();
  Parameters params;

  params.potential_options.potential = PotType::SW_SI;
  auto pot_sw = eonc::helpers::makePotential(params);
  auto m1 = std::make_shared<Matter>(pot_sw, params);
  m1->con2matter(std::string("pos.con"));
  double e_sw = m1->getPotentialEnergy();
  require_loaded_pot(e_sw, "sw_si");

  params.potential_options.potential = PotType::TERSOFF_SI;
  auto pot_tersoff = eonc::helpers::makePotential(params);
  auto m2 = std::make_shared<Matter>(pot_tersoff, params);
  m2->con2matter(std::string("pos.con"));
  double e_tersoff = m2->getPotentialEnergy();
  require_loaded_pot(e_tersoff, "tersoff_si");

  REQUIRE(e_sw != e_tersoff);
}

TEST_CASE("SW minimization energy matches SVN", "[pot][sw][si][minimization]") {
  require_pos_con();
  Parameters params;
  params.potential_options.potential = PotType::SW_SI;
  params.optimizer_options.method = OptType::LBFGS;
  params.optimizer_options.converged_force = 0.001;
  params.optimizer_options.max_iterations = 200;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("pos.con"));

  matter->relax(false, false, false, "sw_test", "sw_test");

  // SVN reference: minimized energy = -16.204961
  {
    double e = matter->getPotentialEnergy();
    require_loaded_pot(e, "sw_si");
    REQUIRE(e == Catch::Approx(-16.204961).epsilon(1e-4));
  }

  TEST_CASE("SW CG minimization matches SVN",
            "[pot][sw][si][minimization][cg]") {
    require_pos_con();
    Parameters params;
    params.potential_options.potential = PotType::SW_SI;
    params.optimizer_options.method = OptType::CG;
    params.optimizer_options.converged_force = 0.001;
    params.optimizer_options.max_iterations = 200;
    auto pot = eonc::helpers::makePotential(params);
    auto matter = std::make_shared<Matter>(pot, params);
    matter->con2matter(std::string("pos.con"));

    matter->relax(false, false, false, "sw_cg_test", "sw_cg_test");

    // SVN reference (data/reference/minimization_sw_cg.dat):
    // energy = -16.204961, 9 force calls
    {
      double e = matter->getPotentialEnergy();
      require_loaded_pot(e, "sw_si");
      REQUIRE(e == Catch::Approx(-16.204961).epsilon(1e-4));
    }
  }

} /* namespace tests */
