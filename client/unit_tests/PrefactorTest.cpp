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

} /* namespace tests */
