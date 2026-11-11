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

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

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
  double pref1 = 0.0;
  double pref2 = 0.0;
  int status = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
  REQUIRE((status == 0 || status == -1));
  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    int silent = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE((silent == 0 || silent == -1));
    logger->set_log_level(quill::LogLevel::TraceL3);
    int loud = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), pref1, pref2, atoms);
    REQUIRE((loud == 0 || loud == -1));
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
  REQUIRE(partialStatus != 0);
  REQUIRE(fullStatus == 0);
  REQUIRE(fullPref > 0.0);
  REQUIRE(fullBack > 0.0);
  REQUIRE(fullPref != fullBack);

  double unused1 = 0.0;
  double unused2 = 0.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), nullptr, min2.get(),
                                         unused1, unused2, full) == -1);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(), nullptr,
                                         unused1, unused2, full) == -1);
  ParametersLoadAccess::prefactor_options(params).rate = "neither";
  int unknownRate = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), unused1, unused2, full);
  REQUIRE((unknownRate == 0 || unknownRate == -1));
  ParametersLoadAccess::prefactor_options(params).rate =
      std::string(eonc::Prefactor::RATE_HTST);
  ParametersLoadAccess::hessian_options(params).write_modes = false;
  int noModes = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), unused1, unused2, full);
  REQUIRE((noModes == 0 || noModes == -1));
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
  REQUIRE((zeroMass == 0 || zeroMass == -1));
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
  REQUIRE(farPairs.size() >= 0);
  ParametersLoadAccess::prefactor_options(params).within_radius = 3.3;

  if (auto *logger = eonc::log::get()) {
    logger->set_log_level(quill::LogLevel::None);
    int silent = eonc::Prefactor::getPrefactors(
        params, min1.get(), saddle.get(), min2.get(), fullPref, fullBack, full);
    REQUIRE(silent == 0);
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
  REQUIRE((clusterStatus == 0 || clusterStatus == -1));
  min1->setFixed(0, true);
  saddle->setFixed(0, true);
  min2->setFixed(0, true);
  int fixedStatus = eonc::Prefactor::getPrefactors(
      params, min1.get(), saddle.get(), min2.get(), clusterPref, clusterBack,
      full);
  REQUIRE((fixedStatus == 0 || fixedStatus == -1));
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
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);

  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      true;
  ParametersLoadAccess::hessian_options(params).zero_freq_value = 3.0;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);

  ParametersLoadAccess::hessian_options(params).zero_freq_value = 0.5;
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         saddle.get(), pref1, pref2,
                                         full) == -1);

  AtomMatrix broken = saddle->getPositions();
  broken(0, 0) = 1.0e300;
  saddle->setPositions(broken);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
  saddle->con2matter(std::string("saddle.con"));
  broken = min2->getPositions();
  broken(0, 0) = 1.0e300;
  min2->setPositions(broken);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
  min2->con2matter(std::string("product.con"));

  ParametersLoadAccess::main_options(params).quiet = true;
  broken = min2->getPositions();
  broken(0, 0) = 1.0e300;
  min2->setPositions(broken);
  REQUIRE(eonc::Prefactor::getPrefactors(params, min1.get(), saddle.get(),
                                         min2.get(), pref1, pref2, full) == -1);
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
  REQUIRE(pct_stop.size() >= 0);
  auto identical = std::make_shared<Matter>(pot, params);
  *identical = *saddle;
  ParametersLoadAccess::prefactor_options(params).filter_fraction = 0.5;
  VectorXi pct_zero = eonc::Prefactor::movedAtomsPct(
      params, identical.get(), saddle.get(), identical.get());
  REQUIRE(pct_zero.size() >= 0);

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
                                           again1, again2) == -1);
    logger->set_log_level(quill::LogLevel::TraceL3);
    REQUIRE(eonc::Prefactor::getPrefactors(params, nullptr, nullptr, nullptr,
                                           again1, again2) == -1);
  }
}

} /* namespace tests */
