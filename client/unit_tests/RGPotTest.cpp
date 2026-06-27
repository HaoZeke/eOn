#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <cmath>
#include <cstdlib>
#include <memory>

using namespace Catch::Matchers;

namespace tests {

namespace {

bool env_nonempty(const char *name) {
  const char *v = std::getenv(name);
  return v != nullptr && v[0] != '\0';
}

} // namespace

TEST_CASE("RGPot nwchemc energy+forces via shipped force()",
          "[PotTest][RGPot][nwchemc]") {
  REQUIRE(env_nonempty("NWCHEMC_LIBRARY") || env_nonempty("RGPOT_NWCHEMC_ENGINE") ||
          env_nonempty("RGPOT_NWCHEM_ENGINE"));

  Parameters params{};
  params.potential_options.potential = PotType::RGPot;
  params.rgpot_options.backend = "nwchemc";
  params.rgpot_options.basis = "sto-3g";
  params.rgpot_options.theory = "scf";
  params.rgpot_options.scf_type = "rhf";
  params.rgpot_options.charge = 0;
  params.rgpot_options.multiplicity = 1;

  auto pot = eonc::helpers::makePotential(params.potential_options.potential,
                                          params);
  REQUIRE(pot != nullptr);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(matter->con2matter("pos.con"));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  bool any_force = false;
  for (long i = 0; i < matter->numberOfAtoms(); ++i) {
    for (int j = 0; j < 3; ++j) {
      REQUIRE(std::isfinite(forces(i, j)));
      if (std::abs(forces(i, j)) > 1e-12) {
        any_force = true;
      }
    }
  }
  REQUIRE(any_force);
}

TEST_CASE("RGPot cpmdc BLYP energy+forces via shipped force()",
          "[PotTest][RGPot][cpmdc]") {
  REQUIRE(env_nonempty("CPMDC_LIBRARY") || env_nonempty("RGPOT_CPMDC_ENGINE") ||
          env_nonempty("RGPOT_CPMD_ENGINE"));

  Parameters params{};
  params.potential_options.potential = PotType::RGPot;
  params.rgpot_options.backend = "cpmdc";
  params.rgpot_options.functional = "BLYP";
  params.rgpot_options.cutoff_ry = 70.0;
  params.rgpot_options.charge = 0;
  params.rgpot_options.multiplicity = 1;

  auto pot = eonc::helpers::makePotential(params.potential_options.potential,
                                          params);
  REQUIRE(pot != nullptr);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(matter->con2matter("pos.con"));
  REQUIRE(matter->numberOfAtoms() >= 3);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  bool any_force = false;
  for (long i = 0; i < matter->numberOfAtoms(); ++i) {
    for (int j = 0; j < 3; ++j) {
      REQUIRE(std::isfinite(forces(i, j)));
      if (std::abs(forces(i, j)) > 1e-12) {
        any_force = true;
      }
    }
  }
  REQUIRE(any_force);
}

} // namespace tests
