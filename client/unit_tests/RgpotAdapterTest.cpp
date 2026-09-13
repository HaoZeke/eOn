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
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Guards the RgpotAdapter + parameter mapping for every migrated arm: the
// same pinned numbers as the pre-migration in-tree kernels, evaluated
// through helpers::makePotential exactly as jobs do. Kernel numerics
// themselves are pinned upstream in rgpot's own test suite.

TEST_CASE("Adapter LJ matches the pinned reference", "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965351).epsilon(1e-4));
  REQUIRE(matter->getForces().rowwise().norm().maxCoeff() ==
          Catch::Approx(0.004704).epsilon(1e-2));
}

TEST_CASE("Adapter LJCluster matches the pinned reference",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJCLUSTER;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965379).epsilon(1e-4));
}

TEST_CASE("Adapter Morse evaluates and reports shared-safe",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(params);

  REQUIRE(pot->isThreadSafe());
  REQUIRE(pot->isSharedInstanceThreadSafe());
  REQUIRE_FALSE(pot->needsPerImageInstance());

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(matter->getForces().allFinite());
}

TEST_CASE("Adapter ZBL maps config.ini cutoffs through",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::ZBL;
  params.zbl_options.cut_inner = 2.0;
  params.zbl_options.cut_global = 2.5;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

// A batch of one repeated system passes even when an implementation
// evaluates the first system and broadcasts the answer, so the two systems
// here differ. LJ takes the default per-system loop rather than a native
// batch, which is the arm every migrated kernel inherits.
TEST_CASE("Adapter forceBatch agrees with per-system force calls",
          "[pot][rgpot-adapter]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(params);

  // Checked, because the displacement below indexes the first atom: an
  // unread file would leave an empty matrix and turn this into a segfault
  // instead of a failure naming the missing input.
  Matter m0(pot, params);
  REQUIRE(eonc::io::io_ok(m0.con2matter(std::string("reactant.con"))));
  Matter m1(pot, params);
  REQUIRE(eonc::io::io_ok(m1.con2matter(std::string("reactant.con"))));
  AtomMatrix displaced = m1.getPositionsCopy();
  displaced(0, 0) += 0.15;
  m1.setPositions(displaced);

  const long nAtoms = m0.numberOfAtoms();
  const auto nrs0 = m0.getAtomicNrs();
  const auto nrs1 = m1.getAtomicNrs();
  const auto cell0 = m0.getCell();
  const auto cell1 = m1.getCell();
  const AtomMatrix pos0 = m0.getPositionsCopy();
  const AtomMatrix pos1 = m1.getPositionsCopy();

  AtomMatrix singleF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix singleF1 = AtomMatrix::Zero(nAtoms, 3);
  double singleU0{0}, singleU1{0}, var{0};
  pot->force(nAtoms, pos0.data(), nrs0.data(), singleF0.data(), &singleU0, &var,
             cell0.data());
  pot->force(nAtoms, pos1.data(), nrs1.data(), singleF1.data(), &singleU1, &var,
             cell1.data());
  // The displacement has to matter, or the comparison below is vacuous.
  REQUIRE(singleU0 != Catch::Approx(singleU1));

  AtomMatrix batchF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix batchF1 = AtomMatrix::Zero(nAtoms, 3);
  const double *positions[2] = {pos0.data(), pos1.data()};
  const int *atomicNrs[2] = {nrs0.data(), nrs1.data()};
  double *forces[2] = {batchF0.data(), batchF1.data()};
  const double *boxes[2] = {cell0.data(), cell1.data()};
  double energies[2] = {0.0, 0.0};
  double variances[2] = {0.0, 0.0};
  pot->forceBatch(2, nAtoms, positions, atomicNrs, forces, energies, variances,
                  boxes);

  REQUIRE(energies[0] == Catch::Approx(singleU0));
  REQUIRE(energies[1] == Catch::Approx(singleU1));
  REQUIRE((batchF0 - singleF0).cwiseAbs().maxCoeff() < 1e-10);
  REQUIRE((batchF1 - singleF1).cwiseAbs().maxCoeff() < 1e-10);
}

#ifdef RGPOT_HAS_EXPR
TEST_CASE("ExprPot rejects empty expression and unknown child kind",
          "[pot][rgpot-adapter][expr]") {
  Parameters empty;
  empty.potential_options.potential = PotType::EXPR;
  REQUIRE_THROWS_WITH(eonc::helpers::makePotential(empty),
                      Catch::Matchers::ContainsSubstring("expression"));
  Parameters bad;
  bad.potential_options.potential = PotType::EXPR;
  bad.expr_options.expression = "x";
  bad.expr_options.terms = "x:not_a_kernel";
  REQUIRE_THROWS_WITH(eonc::helpers::makePotential(bad),
                      Catch::Matchers::ContainsSubstring("unknown term kind"));
}

TEST_CASE("ExprPot 0.5*lj + morse matches independent sum",
          "[pot][rgpot-adapter][expr]") {
  Parameters params;
  params.potential_options.potential = PotType::EXPR;
  params.expr_options.expression = "0.5*lj + morse";
  params.expr_options.terms = "lj:lj, morse:morse";
  auto expr = eonc::helpers::makePotential(params);

  Parameters lj_p;
  lj_p.potential_options.potential = PotType::LJ;
  auto lj = eonc::helpers::makePotential(lj_p);

  Parameters morse_p;
  morse_p.potential_options.potential = PotType::MORSE_PT;
  auto morse = eonc::helpers::makePotential(morse_p);

  auto matter = std::make_shared<Matter>(expr, params);
  REQUIRE(eonc::io::io_ok(matter->con2matter(std::string("reactant.con"))));
  const long nAtoms = matter->numberOfAtoms();
  const AtomMatrix pos = matter->getPositionsCopy();
  const auto nrs = matter->getAtomicNrs();
  const auto cell = matter->getCell();

  AtomMatrix Fe = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix Flj = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix Fm = AtomMatrix::Zero(nAtoms, 3);
  double Ue{0}, Ulj{0}, Um{0}, var{0};
  expr->force(nAtoms, pos.data(), nrs.data(), Fe.data(), &Ue, &var, cell.data());
  lj->force(nAtoms, pos.data(), nrs.data(), Flj.data(), &Ulj, &var, cell.data());
  morse->force(nAtoms, pos.data(), nrs.data(), Fm.data(), &Um, &var,
               cell.data());

  REQUIRE(Ue == Catch::Approx(0.5 * Ulj + Um).epsilon(1e-12));
  REQUIRE((Fe - 0.5 * Flj - Fm).cwiseAbs().maxCoeff() < 1e-10);
}

#ifdef RGPOT_HAS_DFTD3
TEST_CASE("Adapter D3 evaluates a two-atom Ar pair",
          "[pot][rgpot-adapter][d3]") {
  Parameters params;
  params.potential_options.potential = PotType::DFTD3;
  auto pot = eonc::helpers::makePotential(params);
  const long nAtoms = 2;
  const double R[] = {0.0, 0.0, 0.0, 1.5, 0.0, 0.0};
  const int Z[] = {18, 18};
  const double box[] = {40.0, 0.0, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 40.0};
  double F[6] = {};
  double U = 0.0;
  double var = 0.0;
  pot->force(nAtoms, R, Z, F, &U, &var, box);
  REQUIRE(std::isfinite(U));
  REQUIRE(std::isfinite(F[0]));
  // Implied D3 of rgpot v3.1.2 half_lj_plus_d3 - 0.5 * identity LJ.
  REQUIRE(U == Catch::Approx(-0.023417417876102353).epsilon(1e-10));
}

TEST_CASE("INI ExprPot 0.5*lj + d3 matches independent sum",
          "[pot][rgpot-adapter][expr]") {
  const auto tmp =
      std::filesystem::temp_directory_path() / "eon_expr_lj_d3.ini";
  {
    std::ofstream f(tmp);
    f << "[Potential]\n"
         "potential = expr\n"
         "\n"
         "[ExprPot]\n"
         "expression = 0.5*lj + d3\n"
         "terms = lj:lj, d3:d3\n"
         "\n"
         "[D3Pot]\n"
         "functional = pbe\n"
         "atm = true\n"
         "damping = bj\n";
  }
  Parameters params;
  REQUIRE(params.load(tmp.string()) == 0);
  std::filesystem::remove(tmp);
  REQUIRE(params.potential_options.potential == PotType::EXPR);
  REQUIRE(params.expr_options.expression == "0.5*lj + d3");
  REQUIRE(params.dftd_options.functional == "pbe");

  auto expr = eonc::helpers::makePotential(params);

  Parameters lj_p;
  lj_p.potential_options.potential = PotType::LJ;
  auto lj = eonc::helpers::makePotential(lj_p);

  Parameters d3_p;
  d3_p.potential_options.potential = PotType::DFTD3;
  d3_p.dftd_options = params.dftd_options;
  auto d3 = eonc::helpers::makePotential(d3_p);

  const long nAtoms = 2;
  const double R[] = {0.0, 0.0, 0.0, 1.5, 0.0, 0.0};
  const int Z[] = {18, 18};
  const double box[] = {40.0, 0.0, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 40.0};
  double F_e[6] = {}, F_lj[6] = {}, F_d3[6] = {};
  double U_e = 0.0, U_lj = 0.0, U_d3 = 0.0, var = 0.0;
  expr->force(nAtoms, R, Z, F_e, &U_e, &var, box);
  lj->force(nAtoms, R, Z, F_lj, &U_lj, &var, box);
  d3->force(nAtoms, R, Z, F_d3, &U_d3, &var, box);

  REQUIRE(U_e == Catch::Approx(0.5 * U_lj + U_d3).epsilon(1e-12));
  for (int i = 0; i < 6; ++i) {
    REQUIRE(F_e[i] == Catch::Approx(0.5 * F_lj[i] + F_d3[i]).epsilon(1e-10));
  }
  // rgpot v3.1.2 CppCore/tests/data/expr/half_lj_plus_d3_{energy,forces}.npy
  REQUIRE(U_e == Catch::Approx(-0.18358553943241468).epsilon(1e-12));
  REQUIRE(F_e[0] == Catch::Approx(0.5793081483943868).epsilon(1e-12));
  REQUIRE(F_e[3] == Catch::Approx(-0.5793081483943868).epsilon(1e-12));
}
#endif
#endif

#ifndef RGPOT_HAS_DFTD3
TEST_CASE("DFTD3 without wrap names the rebuild flag",
          "[pot][rgpot-adapter][d3]") {
  Parameters params;
  params.potential_options.potential = PotType::DFTD3;
  REQUIRE_THROWS_WITH(eonc::helpers::makePotential(params),
                      Catch::Matchers::ContainsSubstring("with_dftd3"));
}
#endif

} /* namespace tests */
