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
#include "eon/MonteCarloJob.h"
#include "eon/BaseStructures.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/MonteCarlo.h"
#include "eon/PotRegistry.h"
#include "magic_enum/magic_enum.hpp"

#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

std::vector<std::string> MonteCarloJob::run(void) {
  std::string posInFilename = eonc::helpers::getRelevantFile("pos.con");
  std::string posOutFilename("out.con");

  std::vector<std::string> returnFiles;

  auto matter = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(matter->con2matter(posInFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", posInFilename);
    throw std::runtime_error("failed to load " + posInFilename);
  }

  MonteCarlo mc = MonteCarlo(matter, params);
  mc.run(params.monte_carlo_options.steps, params.main_options.temperature,
         params.monte_carlo_options.step_size);

  QUILL_LOG_DEBUG(log, "Saving result to {}", posOutFilename);
  if (eonc::io::io_ok(matter->matter2con(posOutFilename))) {
    returnFiles.push_back(posOutFilename);
  } else {
    QUILL_LOG_ERROR(log, "Failed to write {}", posOutFilename);
  }

  std::string resultsFilename("results.dat");

  std::ofstream out(resultsFilename, std::ios::binary);
  if (!out) {
    QUILL_LOG_CRITICAL(log, "Failed to open {}", resultsFilename);
    throw std::runtime_error("failed to open " + resultsFilename);
  }
  out << std::format("{} termination_reason\n",
                     static_cast<int>(RunStatus::GOOD));
  out << std::format("{} termination_reason_text\n",
                     magic_enum::enum_name<RunStatus>(RunStatus::GOOD));
  out << "monte_carlo job_type\n";
  out << std::format(
      "{} potential_type\n",
      magic_enum::enum_name<PotType>(params.potential_options.potential));
  out << std::format("{} total_force_calls\n",
                     PotRegistry::get().total_force_calls());
  out << std::format("{:f} potential_energy\n", matter->getPotentialEnergy());
  out.close();
  if (!out) {
    QUILL_LOG_CRITICAL(log, "Failed to write {}", resultsFilename);
    throw std::runtime_error("failed to write " + resultsFilename);
  }
  returnFiles.push_back(resultsFilename);

  return returnFiles;
}
