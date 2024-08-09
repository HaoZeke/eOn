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
#include "catch2/catch_amalgamated.hpp"
#include "client/matter/MaskManager.hpp"

using namespace eonc::mat;

TEST_CASE("AtomMask: Basic Functionality", "[AtomMask]") {
  eonc::Vector<size_t> atomIndices(5);
  atomIndices << 0, 1, 2, 3, 4;
  eonc::Vector<bool> freeMask(5);
  freeMask << true, false, true, true, false;

  AtomMask mask(atomIndices, freeMask);

  REQUIRE(mask.atomIndices == atomIndices);
  REQUIRE(mask.freeMask == freeMask);
  REQUIRE(mask.freeIndices == std::vector<int>{0, 2, 3});
}

TEST_CASE("AtomMask: Inconsistent Sizes", "[AtomMask]") {
  eonc::Vector<size_t> atomIndices(5);
  atomIndices << 0, 1, 2, 3, 4;
  eonc::Vector<bool> freeMask(3);
  freeMask << true, false, true;

  REQUIRE_THROWS_AS(AtomMask(atomIndices, freeMask), std::invalid_argument);
}

TEST_CASE("FreeFixer: Filter Positions", "[FreeFixer]") {
  Eigen::MatrixXd positions(5, 3);
  positions << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(5);
  atomIndices << 0, 1, 2, 3, 4;
  eonc::Vector<bool> freeMask(5);
  freeMask << true, false, true, true, false;

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);

  Eigen::MatrixXd filtered_positions = freef(positions);

  REQUIRE(filtered_positions.rows() == 3);
  REQUIRE(filtered_positions.cols() == 3);

  Eigen::MatrixXd expected_positions(3, 3);
  expected_positions << 1, 2, 3, 7, 8, 9, 10, 11, 12;

  REQUIRE(filtered_positions.isApprox(expected_positions));
}

TEST_CASE("FreeFixer: Empty Free Mask", "[FreeFixer]") {
  Eigen::MatrixXd positions(5, 3);
  positions << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(5);
  atomIndices << 0, 1, 2, 3, 4;
  eonc::Vector<bool> freeMask(5);
  freeMask << false, false, false, false, false;

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);

  Eigen::MatrixXd filtered_positions = freef(positions);

  REQUIRE(filtered_positions.rows() == 0);
  REQUIRE(filtered_positions.cols() == positions.cols());
}

TEST_CASE("FreeFixer: All Free Atoms", "[FreeFixer]") {
  Eigen::MatrixXd positions(5, 3);
  positions << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(5);
  atomIndices << 0, 1, 2, 3, 4;
  eonc::Vector<bool> freeMask(5);
  freeMask << true, true, true, true, true;

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);

  Eigen::MatrixXd filtered_positions = freef(positions);

  REQUIRE(filtered_positions.rows() == positions.rows());
  REQUIRE(filtered_positions.cols() == positions.cols());

  REQUIRE(filtered_positions.isApprox(positions));
}

TEST_CASE("FreeFixer: Vector Input", "[FreeFixer]") {
  Eigen::VectorXd velocities(15);
  velocities << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(7);
  atomIndices << 0, 1, 7, 3, 14, 8, 11;
  eonc::Vector<bool> freeMask(7);
  freeMask << true, false, true, true, false, true, true;

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);

  Eigen::VectorXd filtered_velocities = freef(velocities);

  REQUIRE(filtered_velocities.size() == 5);

  Eigen::VectorXd expected_velocities(5);
  expected_velocities << 1, 8, 15, 9, 12;

  REQUIRE(filtered_velocities == expected_velocities);
}
