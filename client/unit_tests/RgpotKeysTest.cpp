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
#include "eon/Parameters.h"

#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("cpmd section supplies cutOffRy", "[params][rgpot]") {
  Parameters canonical;
  REQUIRE(canonical.load_ini_text("[Main]\n"
                                  "job = point\n"
                                  "[Potential]\n"
                                  "potential = RGPOT\n"
                                  "[RgpotPot]\n"
                                  "backend = cpmdc\n"
                                  "cutOffRy = 10\n"
                                  "charge = 1\n"
                                  "input_block = FROM_SHARED\n"
                                  "[cpmd]\n"
                                  "functional = PBE\n"
                                  "cutOffRy = 55.5\n"
                                  "charge = 4\n"
                                  "input_block = DEMO_BLOCK_SENTINEL\n") == 0);
  REQUIRE(canonical.rgpot_options().functional == "PBE");
  REQUIRE(canonical.rgpot_options().cutoff_ry == Catch::Approx(55.5));
  REQUIRE(canonical.rgpot_options().charge == 4);
  REQUIRE(canonical.rgpot_options().input_block == "DEMO_BLOCK_SENTINEL");

  Parameters legacy;
  REQUIRE(legacy.load_ini_text("[Main]\n"
                               "job = point\n"
                               "[Potential]\n"
                               "potential = RGPOT\n"
                               "[RgpotPot]\n"
                               "backend = cpmdc\n"
                               "cutoff_ry = 40\n") == 0);
  REQUIRE(legacy.rgpot_options().cutoff_ry == Catch::Approx(40.0));

  Parameters older;
  REQUIRE(older.load_ini_text("[Main]\n"
                              "job = point\n"
                              "[Potential]\n"
                              "potential = RGPOT\n"
                              "[RgpotPot]\n"
                              "backend = cpmdc\n"
                              "cpmd_cut_off_ry = 33\n") == 0);
  REQUIRE(older.rgpot_options().cutoff_ry == Catch::Approx(33.0));

  Parameters other_backend;
  REQUIRE(other_backend.load_ini_text("[Main]\n"
                                      "job = point\n"
                                      "[Potential]\n"
                                      "potential = RGPOT\n"
                                      "[RgpotPot]\n"
                                      "backend = nwchemc\n"
                                      "[cpmd]\n"
                                      "cutOffRy = 12.5\n"
                                      "input_block = SHOULD_NOT_APPLY\n") == 0);
  REQUIRE(other_backend.rgpot_options().cutoff_ry == Catch::Approx(70.0));
  REQUIRE(other_backend.rgpot_options().input_block.empty());
}

} // namespace tests
