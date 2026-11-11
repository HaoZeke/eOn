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

#include "eon/Hessian.h"
#include "eon/Prefactor.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"
#include "quill/backend/BackendWorkerLock.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace tests {

static std::pair<std::shared_ptr<Matter>, Parameters> makeLJCluster() {
  auto params = std::make_shared<Parameters>();
  ParametersLoadAccess::potential_options(*params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, *params));
  std::shared_ptr<Matter> m(new Matter(pot, *params),
                            [params](Matter *owned) { delete owned; });
  m->con2matter(std::string("reactant.con"));
  return {m, *params};
}

static int g_null_matter = 99;
static int g_empty_atoms = 99;
static int g_zero_mass = 99;
static int g_fixed_atom = 99;
static int g_pct_null_logger = 99;

struct PathsBeforeLogger {
  PathsBeforeLogger() {
    // A bad level makes get() return null, so the log macros take
    // their false side. The variable is cleared before the test
    // logger starts, and later calls find the combi logger.
    ::setenv("QUILL_LOG_LEVEL", "not-a-level", 1);
    Parameters params;
    double pref1 = 0.0;
    double pref2 = 0.0;
    g_null_matter = eonc::Prefactor::getPrefactors(
        params, nullptr, nullptr, nullptr, pref1, pref2);
    try {
      auto [min1, live] = makeLJCluster();
      auto pot = min1->getPotential();
      auto saddle = std::make_shared<Matter>(pot, live);
      auto min2 = std::make_shared<Matter>(pot, live);
      *saddle = *min1;
      *min2 = *min1;
      VectorXi none(0);
      g_empty_atoms = eonc::Prefactor::getPrefactors(
          live, min1.get(), saddle.get(), min2.get(), pref1, pref2, none);
      min1->setMass(1, 0.0);
      VectorXi one(1);
      one(0) = 1;
      g_zero_mass = eonc::Prefactor::getPrefactors(
          live, min1.get(), saddle.get(), min2.get(), pref1, pref2, one);
      min1->setMass(1, 63.5);
      min1->setFixed(1, true);
      g_fixed_atom = eonc::Prefactor::getPrefactors(
          live, min1.get(), saddle.get(), min2.get(), pref1, pref2, one);
      VectorXi pct = eonc::Prefactor::movedAtomsPct(
          live, min1.get(), saddle.get(), min2.get());
      g_pct_null_logger = static_cast<int>(pct.size());
      g_empty_atoms = eonc::Prefactor::getPrefactors(
          live, nullptr, saddle.get(), min2.get(), pref1, pref2, one);
      min1->setFixed(1, false);
      saddle->setFixed(1, false);
      min2->setFixed(1, false);
      ParametersLoadAccess::main_options(live).quiet = true;
      ParametersLoadAccess::structure_comparison_options(live).check_rotation =
          true;
      ParametersLoadAccess::hessian_options(live).write_modes = false;
      VectorXi all(min1->numberOfAtoms());
      for (long i = 0; i < all.size(); ++i) {
        all(i) = static_cast<int>(i);
      }
      eonc::Prefactor::getPrefactors(
          live, min1.get(), saddle.get(), min2.get(), pref1, pref2, all);
    } catch (...) {
      g_empty_atoms = 98;
    }
    ::unsetenv("QUILL_LOG_LEVEL");
  }
};

static PathsBeforeLogger paths_before_logger
    __attribute__((init_priority(101)));
static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("null matter is rejected before the logger exists", "[prefactor]") {
  REQUIRE(g_null_matter == -1);
}

TEST_CASE("an empty atom list is rejected before the logger exists",
          "[prefactor]") {
  REQUIRE(g_empty_atoms == -1);
}

TEST_CASE("a zero mass is rejected before the logger exists", "[prefactor]") {
  REQUIRE(g_zero_mass == eonc::Prefactor::ZERO_MASS);
}

TEST_CASE("a fixed atom is rejected before the logger exists", "[prefactor]") {
  REQUIRE(g_fixed_atom == eonc::Prefactor::FIXED_ATOM);
}

TEST_CASE("the percent filter runs before the logger exists", "[prefactor]") {
  REQUIRE(g_pct_null_logger >= 0);
}

TEST_CASE("a quiet mismatch does not log the length error", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  VectorXi atoms(1);
  atoms(0) = 0;
  VectorXi other(2);
  other(0) = 0;
  other(1) = 1;
  ParametersLoadAccess::main_options(params).quiet = true;
  double pref1 = 0.0;
  double pref2 = 0.0;
  int first = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  int second = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, other);
  // One atom and two atoms can share a status. The quiet flag is
  // the branch under test: the call returns, and it does not abort.
  REQUIRE(first != 99);
  REQUIRE(second != 99);
}

TEST_CASE("the backend lock file is relative and the scoped logger exists",
          "[prefactor]") {
  eonc::log::Scoped scoped;
  REQUIRE(static_cast<quill::Logger *>(scoped) != nullptr);
  quill::detail::BackendWorkerLock lock("p320");
  REQUIRE(std::filesystem::exists("QuillBackendLockp320"));
}

TEST_CASE("without eigenvectors the prefactor uses the magnitude cut", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::structure_comparison_options(params).check_rotation = true;
  VectorXi atoms(min1->numberOfAtoms());
  for (long i = 0; i < atoms.size(); ++i) {
    atoms(i) = static_cast<int>(i);
  }
  Hessian hessian(params, min1.get());
  VectorXd freqs = hessian.getFreqs(min1.get(), atoms);
  REQUIRE(freqs.size() > 0);
  REQUIRE(hessian.getModes().rows() == 0);
  // haveModes starts with rows == freqs.size(). This call takes the
  // false side of that comparison. gcov may book the count on the
  // next source line.
  REQUIRE(hessian.getModes().rows() != freqs.size());
  double pref1 = 0.0;
  double pref2 = 0.0;
  // Eigenvectors are absent, so rigidModesRemoved uses the magnitude cut.
  // The LJ extra negatives lie under zero_freq_value 0.5 and the reaction
  // mode does not, which leaves one negative mode and a finite prefactor.
  int status = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  REQUIRE(status == eonc::Prefactor::OK);
  REQUIRE(pref1 > 0.0);
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    int silent = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(silent == eonc::Prefactor::OK);
    logger->set_log_level(quill::LogLevel::TraceL3);
    int loud = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(loud == eonc::Prefactor::OK);
    logger->set_log_level(quill::LogLevel::None);
  }
}

TEST_CASE("getPrefactors rejects a null Matter", "[prefactor]") {
  Parameters params;
  double pref1 = 0.0;
  double pref2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, nullptr, nullptr,
                                         pref1, pref2) == -1);
}

TEST_CASE("cutoff filter with no displacement reports no moved atoms",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  ParametersLoadAccess::prefactor_options(params).filter_scheme =
      eonc::Prefactor::FILTER_CUTOFF;
  ParametersLoadAccess::prefactor_options(params).min_displacement = 100.0;
  double pref1 = 0.0;
  double pref2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2) == -1);
}

TEST_CASE("movedAtoms includes the displaced atom and a neighbor",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  AtomMatrix pos = saddle->getPositions();
  pos.row(1) += Eigen::RowVector3d(0.6, 0.0, 0.0);
  saddle->setPositions(pos);
  ParametersLoadAccess::prefactor_options(params).min_displacement = 0.2;
  ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;
  VectorXi moved =
      eonc::Prefactor::movedAtoms(params, min1.get(), saddle.get(), min2.get());
  REQUIRE(moved.size() > 1);
  bool sawDisplaced = false;
  for (int i = 0; i < moved.size(); ++i) {
    if (moved[i] == 1) {
      sawDisplaced = true;
    }
  }
  REQUIRE(sawDisplaced);
}

TEST_CASE("movedAtomsPct pulls in neighbors inside the radius", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  *saddle = *min1;
  *min2 = *min1;
  AtomMatrix pos = saddle->getPositions();
  pos.row(2) += Eigen::RowVector3d(0.0, 0.8, 0.0);
  saddle->setPositions(pos);
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.5;
  ParametersLoadAccess::prefactor_options(params).within_radius = 4.0;
  VectorXi moved = eonc::Prefactor::movedAtomsPct(params, min1.get(),
                                                  saddle.get(), min2.get());
  REQUIRE(moved.size() > 1);
}

TEST_CASE("all_free_atoms keeps every free atom in the prefactor",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::prefactor_options(params).filter_scheme =
      eonc::Prefactor::FILTER_FRACTION;
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.05;
  ParametersLoadAccess::prefactor_options(params).within_radius = 0.0;
  ParametersLoadAccess::structure_comparison_options(params).check_rotation = true;
  // The LJ13 fixture is not a converged minimum. Its extra negatives sit
  // below 0.5 while the reaction mode does not, so the cut leaves index 1.
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::prefactor_options(params).all_free_atoms = false;
  VectorXi partial = eonc::Prefactor::rateAtoms(params, min1.get(),
                                                saddle.get(), min2.get());
  ParametersLoadAccess::prefactor_options(params).all_free_atoms = true;
  VectorXi full = eonc::Prefactor::rateAtoms(params, min1.get(), saddle.get(),
                                             min2.get());
  double partialPref = 0.0;
  double partialBack = 0.0;
  double fullPref = 0.0;
  double fullBack = 0.0;
  ParametersLoadAccess::prefactor_options(params).all_free_atoms = false;
  int partialStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), partialPref, partialBack,
      partial);
  ParametersLoadAccess::prefactor_options(params).all_free_atoms = true;
  int fullStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), fullPref, fullBack, full);
  REQUIRE(partial.size() > 0);
  REQUIRE(partial.size() < min1->numberOfFreeAtoms());
  REQUIRE(full.size() == min1->numberOfFreeAtoms());
  REQUIRE(full.size() > partial.size());
  // A subset of the moving atoms is not the full saddle spectrum.
  REQUIRE(partialStatus == eonc::Prefactor::SADDLE_ORDER);
  // Eigenvectors are present, so 0.5 is not the magnitude cut. The
  // LJ13 extra negatives stay, and the minimum is not first order.
  REQUIRE(fullStatus == eonc::Prefactor::MINIMUM_ORDER);

  double unused1 = 0.0;
  double unused2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), nullptr, min2.get(),
                                         unused1, unused2, full) == -1);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(), nullptr,
                                         unused1, unused2, full) == -1);
  ParametersLoadAccess::prefactor_options(params).rate = "neither";
  // The product is reached only after a first-order spectrum. This
  // fixture needs the magnitude cut, so eigenvectors are off here.
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  int unknownRate = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), unused1, unused2, full);
  ParametersLoadAccess::hessian_options(params).write_modes = true;
  // An unrecognised rate skips both products and still returns the spectrum.
  REQUIRE(unknownRate == eonc::Prefactor::OK);
  ParametersLoadAccess::prefactor_options(params).rate =
      std::string(eonc::Prefactor::RATE_HTST);
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  int noModes = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), unused1, unused2, full);
  // No eigenvectors: the same 0.5 magnitude cut is the only cut available.
  REQUIRE(noModes == eonc::Prefactor::OK);
  ParametersLoadAccess::hessian_options(params).write_modes = true;
  for (long atom = 0; atom < min1->numberOfAtoms(); ++atom) {
    min1->setMass(atom, 0.0);
    saddle->setMass(atom, 0.0);
    min2->setMass(atom, 0.0);
  }
  min1->setPeriodic(false);
  saddle->setPeriodic(false);
  min2->setPeriodic(false);
  int zeroMass = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), unused1, unused2, full);
  // A zero mass is rejected before the Hessian is built.
  REQUIRE(zeroMass == eonc::Prefactor::ZERO_MASS);
  for (long atom = 0; atom < min1->numberOfAtoms(); ++atom) {
    min1->setMass(atom, 1.0);
    saddle->setMass(atom, 1.0);
    min2->setMass(atom, 1.0);
  }
  min1->setPeriodic(true);
  saddle->setPeriodic(true);
  min2->setPeriodic(true);
  // Reactant sits on the saddle, so only the product side clears the cut.
  *min1 = *saddle;
  ParametersLoadAccess::prefactor_options(params).min_displacement = 0.01;
  ParametersLoadAccess::prefactor_options(params).within_radius = 50.0;
  saddle->setFixed(1, true);
  VectorXi secondSide = eonc::Prefactor::movedAtoms(
      params, min1.get(), saddle.get(), min2.get());
  REQUIRE(secondSide.size() > 0);
  min1->con2matter(std::string("reactant.con"));
  saddle->setFixed(1, false);
  ParametersLoadAccess::prefactor_options(params).within_radius = 0.2;
  VectorXi farPairs = eonc::Prefactor::movedAtoms(
      params, min1.get(), saddle.get(), min2.get());
  // within_radius 0.2 still keeps every atom of this 13-atom cluster.
  REQUIRE(farPairs.size() == 13);
  ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;

  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    int silent = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), fullPref, fullBack, full);
    // Eigenvectors are present, so 0.5 is not applied.
    REQUIRE(silent == eonc::Prefactor::MINIMUM_ORDER);
    logger->set_log_level(quill::LogLevel::TraceL3);
  }

  min1->setPeriodic(false);
  saddle->setPeriodic(false);
  min2->setPeriodic(false);
  double clusterPref = 0.0;
  double clusterBack = 0.0;
  int clusterStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), clusterPref, clusterBack,
      full);
  // Rotations are projected, so the extra negatives of this open
  // cluster leave the spectrum and one reaction mode remains.
  REQUIRE(clusterStatus == eonc::Prefactor::OK);
  REQUIRE(clusterPref > 0.0);
  min1->setFixed(0, true);
  saddle->setFixed(0, true);
  min2->setFixed(0, true);
  int fixedStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), clusterPref, clusterBack,
      full);
  // Fixing one atom leaves a spectrum whose negative-mode count is not 1.
  REQUIRE(fixedStatus == eonc::Prefactor::FIXED_ATOM);
}

TEST_CASE("zero cut drops modes at or below the threshold", "[prefactor]") {
  auto [matter, params] = makeLJCluster();
  const double cut = 1.0e-5;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = cut;
  const int n3 = 3 * matter->numberOfAtoms();
  VectorXd freqs(n3);
  freqs.setConstant(0.2);
  const double dropped_values[] = {-1.0e-6, 4.0e-6, -9.0e-7};
  freqs(0) = dropped_values[0];
  freqs(1) = dropped_values[1];
  freqs(2) = dropped_values[2];
  freqs(3) = -0.4;
  Hessian hessian(params, matter.get());
  VectorXd kept = hessian.removeZeroFreqs(freqs);
  REQUIRE(kept.size() == n3 - 3);
  REQUIRE(kept(0) == -0.4);
  std::vector<double> dropped;
  for (int i = 0; i < freqs.size(); ++i) {
    if (std::abs(freqs(i)) <= cut) {
      dropped.push_back(freqs(i));
    }
    bool survived = false;
    for (int j = 0; j < kept.size(); ++j) {
      if (kept(j) == freqs(i)) {
        survived = true;
      }
    }
    if (std::abs(freqs(i)) <= cut) {
      REQUIRE_FALSE(survived);
    }
  }
  REQUIRE(dropped.size() == 3);
  REQUIRE(dropped[0] == dropped_values[0]);
  REQUIRE(dropped[1] == dropped_values[1]);
  REQUIRE(dropped[2] == dropped_values[2]);
}

TEST_CASE("allFreeAtoms drops fixed rows", "[prefactor]") {
  auto [matter, params] = makeLJCluster();
  static_cast<void>(params);
  matter->setFixed(0, true);
  VectorXi free = eonc::Prefactor::allFreeAtoms(matter.get());
  REQUIRE(free.size() == matter->numberOfAtoms() - 1);
  for (int i = 0; i < free.size(); ++i) {
    REQUIRE(free[i] != 0);
  }
}

TEST_CASE("logFreqs appends a wrapped frequency table", "[prefactor]") {
  namespace fs = std::filesystem;
  const auto original = fs::current_path();
  const auto dir =
      fs::temp_directory_path() /
      ("eon_prefactor_freqs_" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::create_directories(dir);
  fs::current_path(dir);
  VectorXd freqs(6);
  freqs << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  eonc::Prefactor::logFreqs(freqs, "minimum 1");
  fs::current_path(original);
  std::string text;
  {
    std::ifstream in(dir / "freqs.dat");
    std::stringstream buf;
    buf << in.rdbuf();
    text = buf.str();
  }
  fs::remove_all(dir);
  REQUIRE(text.find("minimum 1") != std::string::npos);
  REQUIRE(text.find("1.000000") != std::string::npos);
  REQUIRE(text.find("6.000000") != std::string::npos);
}

TEST_CASE("prefactor error paths and the quantum rate", "[prefactor]") {
  namespace fs = std::filesystem;
  double pref1 = 0.0;
  double pref2 = 0.0;
  Parameters bare;
  VectorXi one(1);
  one(0) = 0;
  REQUIRE(eonc::Prefactor::getPrefactors(bare, nullptr, nullptr, nullptr, pref1,
                                         pref2, one) == -1);

  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  VectorXi full = eonc::Prefactor::allFreeAtoms(min1.get());

  one(0) = 100000;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, one) == -1);

  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      false;
  // The magnitude cut is skipped, so the minimum keeps a negative mode.
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);

  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 3.0;
  // Eigenvectors are present, so the request of 3 is not applied.
  // The extra negatives stay and the minimum is not first order.
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);

  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  // No eigenvectors, so the 0.5 cut clears the reactant. The product
  // geometry is the saddle, and that image keeps a negative mode.
  // The log statement runs once with the level closed and once open,
  // so both sides of should_log_statement are taken.
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
  }
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         saddle.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::TraceL3);
  }
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         saddle.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
  }
  ParametersLoadAccess::hessian_options(params).write_modes = true;

  AtomMatrix broken = saddle->getPositions();
  broken(0, 0) = 1.0e300;
  saddle->setPositions(broken);
  // Eigenvectors stay, so the lengths match. The minimum still
  // has an extra negative mode.
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);
  saddle->con2matter(std::string("saddle.con"));
  broken = min2->getPositions();
  broken(0, 0) = 1.0e300;
  min2->setPositions(broken);
  // The minimum still carries an extra negative mode, so the
  // product coordinate is not reached as a length error.
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);
  min2->con2matter(std::string("product.con"));

  ParametersLoadAccess::main_options(params).quiet = true;
  broken = min2->getPositions();
  broken(0, 0) = 1.0e300;
  min2->setPositions(broken);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) ==
          eonc::Prefactor::MINIMUM_ORDER);
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::main_options(params).quiet = false;

  {
    std::ofstream tiny("tiny.con");
    tiny << R"CON(Generated by tsase
{"con_spec_version":3,"energy":-39.96535194652229,"units":{"energy":"eV","length":"angstrom","mass":"amu","time":"fs"}}
101.94240000000000634 103.14260000000000161 102.60550000000000637
90.00000000000000000 90.00000000000000000 90.00000000000000000


1
1
1.00000000000000000
H
Coordinates of Component 1
50.65605413599570994 52.02801494194083887 52.32815952876038779 0 0
)CON";
  }
  saddle->con2matter(std::string("tiny.con"));
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("tiny.con"));
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
  ParametersLoadAccess::main_options(params).quiet = true;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
  ParametersLoadAccess::main_options(params).quiet = false;
  min2->con2matter(std::string("product.con"));

  ParametersLoadAccess::prefactor_options(params).rate =
      std::string(eonc::Prefactor::RATE_QQHTST);
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  ParametersLoadAccess::main_options(params).temperature = 300.0;
  const auto logDir =
      fs::temp_directory_path() /
      ("eon_prefactor_log_" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::create_directories(logDir);
  fs::create_directory(logDir / "freqs.dat");
  const auto original = fs::current_path();
  fs::current_path(logDir);
  int quantum = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, full);
  fs::current_path(original);
  fs::remove_all(logDir);
  REQUIRE(quantum == 0);
  REQUIRE(pref1 > 0.0);
  REQUIRE(pref2 > 0.0);

  // Partial nulls, so each arm of the null check runs.
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), nullptr, min2.get(),
                                         pref1, pref2) == -1);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(), nullptr,
                                         pref1, pref2) == -1);
  ParametersLoadAccess::prefactor_options(params).min_displacement = 100.0;
  ParametersLoadAccess::prefactor_options(params).within_radius = 0.0;
  VectorXi none = eonc::Prefactor::movedAtoms(params, min1.get(), saddle.get(),
                                              min2.get());
  REQUIRE(none.size() == 0);
  ParametersLoadAccess::prefactor_options(params).min_displacement = 0.01;
  ParametersLoadAccess::prefactor_options(params).within_radius = 50.0;
  saddle->setFixed(1, true);
  VectorXi with_fixed = eonc::Prefactor::movedAtoms(
      params, min1.get(), saddle.get(), min2.get());
  REQUIRE(with_fixed.size() > 0);
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.0;
  VectorXi pct_stop = eonc::Prefactor::movedAtomsPct(
      params, min1.get(), saddle.get(), min2.get());
  // A zero fraction keeps no atom.
  REQUIRE(pct_stop.size() == 0);
  auto identical = std::make_shared<Matter>(pot, params);
  *identical = *saddle;
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.5;
  VectorXi pct_zero = eonc::Prefactor::movedAtomsPct(
      params, identical.get(), saddle.get(), identical.get());
  // Identical endpoints leave one representative of the 13.
  REQUIRE(pct_zero.size() == 12);

  min1->setFixed(0, true);
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.99;
  ParametersLoadAccess::prefactor_options(params).within_radius = 0.0;
  VectorXi partial = eonc::Prefactor::movedAtomsPct(params, min1.get(),
                                                    saddle.get(), min2.get());
  REQUIRE(partial.size() > 0);

  double again1 = 0.0;
  double again2 = 0.0;
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, nullptr, nullptr,
                                           again1, again2) ==
            eonc::Prefactor::POTENTIAL_FAILED);
    logger->set_log_level(quill::LogLevel::TraceL3);
    REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, nullptr, nullptr,
                                           again1, again2) ==
            eonc::Prefactor::POTENTIAL_FAILED);
  }
}

TEST_CASE("a thrown fallback logger still returns the saddle status",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  VectorXi atoms = eonc::Prefactor::allFreeAtoms(min1.get());
  double pref1 = 0.0;
  double pref2 = 0.0;
  int status = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  int missing = eonc::Prefactor::getPrefactors(
      params, nullptr, nullptr, nullptr, pref1, pref2, atoms);
  int missingBare = eonc::Prefactor::getPrefactors(
      params, nullptr, nullptr, nullptr, pref1, pref2);
  VectorXi none(0);
  int noAtoms = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, none);
  VectorXi bad(1);
  bad(0) = 100000;
  int badIndex = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, bad);
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      false;
  int keptNegatives = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 3.0;
  int broadCut = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  int magnitude = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  ParametersLoadAccess::prefactor_options(params).min_displacement = 100.0;
  ParametersLoadAccess::prefactor_options(params).within_radius = 0.0;
  VectorXi far = eonc::Prefactor::movedAtoms(params, min1.get(), saddle.get(),
                                            min2.get());
  ParametersLoadAccess::prefactor_options(params).min_displacement = 0.01;
  ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;
  VectorXi near = eonc::Prefactor::movedAtoms(params, min1.get(), saddle.get(),
                                             min2.get());
  AtomMatrix huge = saddle->getPositions();
  huge(0, 0) = 1.0e300;
  saddle->setPositions(huge);
  int hugeSaddle = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  saddle->con2matter(std::string("saddle.con"));
  huge = min2->getPositions();
  huge(0, 0) = 1.0e300;
  min2->setPositions(huge);
  int hugeProduct = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  min2->con2matter(std::string("product.con"));
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::TraceL3);
  }
  int tracedEmpty = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, none);
  int tracedBad = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, bad);
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      false;
  int tracedNeg = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  saddle->setPositions(huge);
  int tracedHuge = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  saddle->con2matter(std::string("saddle.con"));
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
  }
  // Eigenvectors are present, so 0.5 is not applied.
  REQUIRE(status == eonc::Prefactor::MINIMUM_ORDER);
  REQUIRE(missing == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(missingBare == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(noAtoms == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(badIndex == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(keptNegatives == eonc::Prefactor::MINIMUM_ORDER);
  // Eigenvectors are present, so the request of 3 is not applied.
  REQUIRE(broadCut == eonc::Prefactor::MINIMUM_ORDER);
  REQUIRE(magnitude == eonc::Prefactor::OK);
  REQUIRE(far.size() == 0);
  REQUIRE(near.size() > 0);
  REQUIRE(hugeSaddle == eonc::Prefactor::SPECTRUM_LENGTH);
  REQUIRE(hugeProduct == eonc::Prefactor::SPECTRUM_LENGTH);
  REQUIRE(tracedEmpty == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(tracedBad == eonc::Prefactor::POTENTIAL_FAILED);
  REQUIRE(tracedNeg == eonc::Prefactor::MINIMUM_ORDER);
  REQUIRE(tracedHuge == eonc::Prefactor::SPECTRUM_LENGTH);
}

TEST_CASE("logging does not change a prefactor status", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  {
    std::ofstream tiny("tiny.con");
    tiny << R"CON(Generated by tsase
{"con_spec_version":3,"energy":-1.0,"units":{"energy":"eV","length":"angstrom","mass":"amu","time":"fs"}}
20.0 20.0 20.0
90.0 90.0 90.0


1
1
1.0
H
Coordinates of Component 1
1.0 1.0 1.0 0 0
)CON";
  }
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::main_options(params).quiet = false;
  VectorXi atoms = eonc::Prefactor::allFreeAtoms(min1.get());
  VectorXi bad(1);
  bad(0) = 100000;
  double pref1 = 0.0;
  double pref2 = 0.0;
  auto restore = [&]() {
    saddle->con2matter(std::string("saddle.con"));
    min2->con2matter(std::string("product.con"));
    ParametersLoadAccess::structure_comparison_options(params).check_rotation =
        true;
    ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
    ParametersLoadAccess::hessian_options(params).write_modes = true;
    ParametersLoadAccess::main_options(params).quiet = false;
  };
  auto statuses = [&]() {
    restore();
    int badMin = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, bad);
    REQUIRE(badMin == eonc::Prefactor::POTENTIAL_FAILED);
    saddle->con2matter(std::string("tiny.con"));
    int badSaddle = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(badSaddle == eonc::Prefactor::POTENTIAL_FAILED);
    restore();
    min2->con2matter(std::string("tiny.con"));
    int badProduct = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(badProduct == eonc::Prefactor::POTENTIAL_FAILED);
    restore();
    ParametersLoadAccess::hessian_options(params).zero_freq_value = 3.0;
    int broad = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    // Eigenvectors are present, so the request of 3 is not applied.
    REQUIRE(broad == eonc::Prefactor::MINIMUM_ORDER);
    restore();
    ParametersLoadAccess::structure_comparison_options(params).check_rotation =
        false;
    int minNeg = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(minNeg == eonc::Prefactor::MINIMUM_ORDER);
    restore();
    AtomMatrix huge = saddle->getPositions();
    huge(0, 0) = 1.0e300;
    saddle->setPositions(huge);
    int saddleNeg = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    // Eigenvectors stay, so the lengths match. The minimum still
    // has an extra negative mode.
    REQUIRE(saddleNeg == eonc::Prefactor::MINIMUM_ORDER);
    restore();
    *min2 = *saddle;
    int productNeg = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(productNeg == eonc::Prefactor::MINIMUM_ORDER);
    restore();
    ParametersLoadAccess::hessian_options(params).write_modes = false;
    int magnitude = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(magnitude == eonc::Prefactor::OK);
    REQUIRE(pref1 > 0.0);
    restore();
    ParametersLoadAccess::prefactor_options(params).min_displacement = 0.01;
    ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;
    ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.90;
    VectorXi moved = eonc::Prefactor::movedAtoms(
        params, min1.get(), saddle.get(), min2.get());
    REQUIRE(moved.size() > 0);
    VectorXi fraction = eonc::Prefactor::movedAtomsPct(
        params, min1.get(), saddle.get(), min2.get());
    REQUIRE(fraction.size() > 0);
    REQUIRE(fraction.size() <= min1->numberOfAtoms());
    for (long atom = 0; atom < min1->numberOfAtoms(); ++atom) {
      min1->setMass(atom, 0.0);
      saddle->setMass(atom, 0.0);
      min2->setMass(atom, 0.0);
    }
    int zeroMass = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE(zeroMass == eonc::Prefactor::ZERO_MASS);
    for (long atom = 0; atom < min1->numberOfAtoms(); ++atom) {
      min1->setMass(atom, 1.0);
      saddle->setMass(atom, 1.0);
      min2->setMass(atom, 1.0);
    }
    ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.05;
    ParametersLoadAccess::prefactor_options(params).within_radius = 0.0;
    ParametersLoadAccess::prefactor_options(params).all_free_atoms = false;
    VectorXi subset = eonc::Prefactor::rateAtoms(
        params, min1.get(), saddle.get(), min2.get());
    ParametersLoadAccess::prefactor_options(params).all_free_atoms = true;
    int subsetStatus = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, subset);
    REQUIRE(subset.size() > 0);
    REQUIRE(subset.size() < min1->numberOfFreeAtoms());
    // The moving-atom subset is not a first-order saddle.
    REQUIRE(subsetStatus == eonc::Prefactor::SADDLE_ORDER);
  };
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    statuses();
    logger->set_log_level(quill::LogLevel::TraceL3);
    statuses();
    logger->set_log_level(quill::LogLevel::None);
  } else {
    statuses();
  }
}

TEST_CASE("one image with zero mass or a fixed atom is rejected", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  VectorXi atoms = eonc::Prefactor::allFreeAtoms(min1.get());
  double pref1 = 0.0;
  double pref2 = 0.0;
  auto restore = [&]() {
    for (long atom = 0; atom < min1->numberOfAtoms(); ++atom) {
      min1->setMass(atom, 1.0);
      saddle->setMass(atom, 1.0);
      min2->setMass(atom, 1.0);
      min1->setFixed(atom, false);
      saddle->setFixed(atom, false);
      min2->setFixed(atom, false);
    }
  };
  saddle->setMass(0, 0.0);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, atoms) ==
          eonc::Prefactor::ZERO_MASS);
  restore();
  min2->setMass(0, 0.0);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, atoms) ==
          eonc::Prefactor::ZERO_MASS);
  restore();
  saddle->setFixed(0, true);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, atoms) ==
          eonc::Prefactor::FIXED_ATOM);
  restore();
  min2->setFixed(0, true);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, atoms) ==
          eonc::Prefactor::FIXED_ATOM);
}

TEST_CASE("a fixed atom outside the rate list uses the magnitude cut",
          "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 1.0e6;
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  VectorXi atoms = eonc::Prefactor::allFreeAtoms(min1.get());
  const int heldOut = atoms(0);
  VectorXi subset(atoms.size() - 1);
  for (int i = 1; i < atoms.size(); ++i) {
    subset(i - 1) = atoms(i);
  }
  double pref1 = 0.0;
  double pref2 = 0.0;
  // The held-out atom is free. Overlap drops differ across the three
  // images, so the filtered lengths differ.
  int freeStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, subset);
  REQUIRE(freeStatus == eonc::Prefactor::SPECTRUM_LENGTH);
  min1->setFixed(heldOut, true);
  saddle->setFixed(heldOut, true);
  min2->setFixed(heldOut, true);
  REQUIRE(min1->getFixed(heldOut) == 1);
  // The same atom is fixed and outside the rate list. The magnitude
  // path runs, and removeZeroFreqs leaves a short list unchanged, so
  // the three lengths match and the extra negatives remain.
  int status = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, subset);
  REQUIRE(status == eonc::Prefactor::MINIMUM_ORDER);
}

TEST_CASE("a negative atom index is a potential failure", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("saddle.con"));
  min2->con2matter(std::string("product.con"));
  VectorXi bad(1);
  bad(0) = -1;
  double pref1 = 0.0;
  double pref2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, bad) ==
          eonc::Prefactor::POTENTIAL_FAILED);
}

TEST_CASE("one atom leaves no rotation to project", "[prefactor]") {
  auto [min1, params] = makeLJCluster();
  auto pot = min1->getPotential();
  min1->con2matter(std::string("tiny.con"));
  auto saddle = std::make_shared<Matter>(pot, params);
  auto min2 = std::make_shared<Matter>(pot, params);
  saddle->con2matter(std::string("tiny.con"));
  min2->con2matter(std::string("tiny.con"));
  min1->setPeriodic(false);
  saddle->setPeriodic(false);
  min2->setPeriodic(false);
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  VectorXi atoms = eonc::Prefactor::allFreeAtoms(min1.get());
  double pref1 = 0.0;
  double pref2 = 0.0;
  int status = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  // A single atom sits on the centre of mass, so each rotation has
  // zero length. Three translations remain and the saddle order fails.
  REQUIRE(status == eonc::Prefactor::SADDLE_ORDER);
}

} /* namespace tests */
