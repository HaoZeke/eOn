#include "../MatrixHelpers.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <cmath>
#include <cstdlib>
#include <memory>

#ifdef WITH_RGPOT
#include "potentials/Rgpot/RgpotPot.h"
#endif

using namespace Catch::Matchers;

namespace tests {

namespace {

bool env_nonempty(const char *name) {
  const char *v = std::getenv(name);
  return v != nullptr && v[0] != '\0';
}

} // namespace

TEST_CASE("RgpotPot in-process nwchemc force (no potserv)",
          "[PotTest][RGPOT][nwchemc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // Engine is a runtime dlopen dep; packaging CI may build WITH_RGPOT without
  // a local libnwchemc — skip rather than fail the default suite.
  if (!(env_nonempty("NWCHEMC_LIBRARY") ||
        env_nonempty("RGPOT_NWCHEMC_ENGINE") ||
        env_nonempty("RGPOT_NWCHEM_ENGINE"))) {
    SKIP("no NWChem embed library in environment "
         "(set NWCHEMC_LIBRARY / RGPOT_NWCHEMC_ENGINE)");
  }

  Parameters params{};
  params.potential_options.potential = PotType::RGPOT;
  params.rgpot_options.backend = "nwchemc";
  params.rgpot_options.basis = "sto-3g";
  params.rgpot_options.theory = "scf";
  params.rgpot_options.scf_type = "rhf";
  params.rgpot_options.charge = 0;
  params.rgpot_options.multiplicity = 1;

  auto pot =
      eonc::helpers::makePotential(params.potential_options.potential, params);
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter("pos.con")));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  // Two successive forces must both succeed (warm multi-call path)
  double energy2 = 0.0;
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy2, nullptr,
             matter->getCell().data());
  REQUIRE(std::isfinite(energy2));
  REQUIRE(std::abs(energy - energy2) < 1e-4);
#endif
}

TEST_CASE("RgpotPot in-process cpmdc force (no potserv)",
          "[PotTest][RGPOT][cpmdc]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // Engine is a runtime dlopen dep; packaging CI may build WITH_RGPOT without
  // a local libcpmdc — skip rather than fail the default suite.
  if (!(env_nonempty("CPMDC_LIBRARY") || env_nonempty("RGPOT_CPMDC_ENGINE") ||
        env_nonempty("RGPOT_CPMD_ENGINE"))) {
    SKIP("no CPMD embed library in environment "
         "(set CPMDC_LIBRARY / RGPOT_CPMDC_ENGINE)");
  }

  Parameters params{};
  params.potential_options.potential = PotType::RGPOT;
  params.rgpot_options.backend = "cpmdc";
  params.rgpot_options.functional = "BLYP";
  params.rgpot_options.cutoff_ry = 70.0;
  params.rgpot_options.charge = 0;
  params.rgpot_options.multiplicity = 1;

  auto pot =
      eonc::helpers::makePotential(params.potential_options.potential, params);
  REQUIRE(pot != nullptr);
  REQUIRE(pot->getType() == PotType::RGPOT);

  auto matter = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter("pos.con")));
  REQUIRE(matter->numberOfAtoms() == 9);

  double energy = 0.0;
  AtomMatrix forces = MatrixXd::Zero(matter->numberOfAtoms(), 3);
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy, nullptr,
             matter->getCell().data());

  REQUIRE(std::isfinite(energy));
  REQUIRE(std::abs(energy) > 1e-6);

  // Two successive forces must both succeed (warm multi-call path)
  double energy2 = 0.0;
  pot->force(matter->numberOfAtoms(), matter->getPositions().data(),
             matter->getAtomicNrs().data(), forces.data(), &energy2, nullptr,
             matter->getCell().data());
  REQUIRE(std::isfinite(energy2));
  REQUIRE(std::abs(energy - energy2) < 1e-4);
#endif
}

TEST_CASE("RgpotPot multi-rank driver instance is reused by makePotential",
          "[PotTest][RGPOT][mpi-driver]") {
#ifndef WITH_RGPOT
  SKIP("built without WITH_RGPOT");
#else
  // Structural lifetime fix: rank 0 must not construct a second pot while
  // workers sit in mpiWorkerLoop. Factory returns the held shared instance.
  // No live NWChem required for the holder null path; engine path optional.
  REQUIRE(RgpotPot::mpiDriverInstance() == nullptr);
  RgpotPot::clearMpiDriverInstance();
  REQUIRE(RgpotPot::mpiDriverInstance() == nullptr);

  // mpiWorldSize is 1 without EONMPI / without MPI_Init — single-rank path.
  REQUIRE(RgpotPot::mpiWorldSize() >= 1);
  REQUIRE(RgpotPot::mpiWorldRank() >= 0);

  if (!(env_nonempty("NWCHEMC_LIBRARY") ||
        env_nonempty("RGPOT_NWCHEMC_ENGINE") ||
        env_nonempty("RGPOT_NWCHEM_ENGINE"))) {
    SKIP("no NWChem embed library; holder null path only");
  }

  Parameters params{};
  params.potential_options.potential = PotType::RGPOT;
  params.rgpot_options.backend = "nwchemc";
  params.rgpot_options.basis = "sto-3g";
  params.rgpot_options.theory = "scf";
  params.rgpot_options.scf_type = "rhf";

  auto driver = std::make_shared<RgpotPot>(params);
  if (!driver->engineAvailable()) {
    SKIP("engine not available at runtime");
  }
  RgpotPot::setMpiDriverInstance(driver);
  REQUIRE(RgpotPot::mpiDriverInstance().get() == driver.get());

  auto from_factory = eonc::helpers::makePotential(
      params.potential_options.potential, params);
  REQUIRE(from_factory.get() == driver.get());

  RgpotPot::clearMpiDriverInstance();
  REQUIRE(RgpotPot::mpiDriverInstance() == nullptr);
  auto fresh = eonc::helpers::makePotential(
      params.potential_options.potential, params);
  REQUIRE(fresh.get() != driver.get());
  REQUIRE(fresh != nullptr);
#endif
}

} // namespace tests
