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

  eonc::Vector<size_t> exp_freeIndices(3);
  exp_freeIndices << 0, 2, 3;
  REQUIRE(mask.freeIndices == exp_freeIndices);
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

  AtomMask mask(atomIndices, freeMask, 3);
  FreeFixer freef(mask);

  Eigen::MatrixXd filtered_positions = freef(positions);

  REQUIRE(filtered_positions.rows() == 3);
  REQUIRE(filtered_positions.cols() == 3);

  Eigen::MatrixXd expected_positions(3, 3);
  expected_positions << 1, 2, 3, 7, 8, 9, 10, 11, 12;

  REQUIRE(filtered_positions == expected_positions);
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

  REQUIRE(filtered_positions == positions);
}

TEST_CASE("FreeFixer: Vector Input", "[FreeFixer]") {
  Eigen::VectorXd velocities(15);
  velocities << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(2);
  atomIndices << 0, 1;
  eonc::Vector<bool> freeMask(2);
  freeMask << true, false;

  Eigen::MatrixXd reshaped_velocities =
      Eigen::Map<Eigen::MatrixXd>(velocities.data(), velocities.size() / 3, 3);

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);
  Eigen::MatrixXd filtered_positions_matrix = freef(reshaped_velocities);

  Eigen::VectorXd filtered_velocities = Eigen::Map<Eigen::VectorXd>(
      filtered_positions_matrix.data(), filtered_positions_matrix.size());

  REQUIRE(filtered_velocities.size() == 3);

  Eigen::VectorXd expected_velocities(3);
  expected_velocities << 1, 2, 3;

  REQUIRE(filtered_velocities == expected_velocities);
}

TEST_CASE("FreeFixer: Out of Bounds Indices", "[FreeFixer]") {
  Eigen::VectorXd velocities(15);
  velocities << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  eonc::Vector<size_t> atomIndices(3);
  atomIndices << 0, 1, 6; // Index 6 is out of bounds for reshaped_velocities
  eonc::Vector<bool> freeMask(3);
  freeMask << true, false, true;

  Eigen::MatrixXd reshaped_velocities =
      Eigen::Map<Eigen::MatrixXd>(velocities.data(), velocities.size() / 3, 3);

  AtomMask mask(atomIndices, freeMask);
  FreeFixer freef(mask);

  // Expecting an out_of_range exception due to out-of-bounds index
  REQUIRE_THROWS_AS(freef(reshaped_velocities), std::out_of_range);
}
