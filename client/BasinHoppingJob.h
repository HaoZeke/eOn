/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
#pragma once
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class BasinHoppingJob : public Job {
public:
  BasinHoppingJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        current{std::make_shared<Matter>(pot, params)},
        trial{std::make_shared<Matter>(pot, params)},
        fcalls{0} {
    log = spdlog::get("combi");
  }
  ~BasinHoppingJob(void) = default;

  std::vector<std::string> run(void) override;

private:
  VectorType calculateDistanceFromCenter(Matter *matter);
  AtomMatrix displaceRandom(double maxDisplacement);
  void randomSwap(Matter *matter);
  std::shared_ptr<Matter> current;
  std::shared_ptr<Matter> trial; // initial configuration
  std::vector<long> getElements(Matter *matter);
  std::vector<std::string> returnFiles;
  int jump_count; // count of jump moves
  int disp_count; // count of displacement moves
  int swap_count; // count of swap moves
  int fcalls;

  std::vector<std::shared_ptr<Matter>> uniqueStructures;
  std::vector<double> uniqueEnergies;
  std::shared_ptr<spdlog::logger> log;
};
