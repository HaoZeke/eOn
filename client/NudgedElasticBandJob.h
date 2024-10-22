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
#include "NudgedElasticBand.h"
#include "Parameters.h"

class NudgedElasticBandJob : public Job {

public:
  NudgedElasticBandJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        fCallsNEB{0} {
    log = spdlog::get("combi");
  }
  ~NudgedElasticBandJob(void) = default;
  std::vector<std::string> run(void);

private:
  // functions
  void printEndState(NudgedElasticBand::NEBStatus status);
  void saveData(NudgedElasticBand::NEBStatus status, NudgedElasticBand *neb);

  // variables
  std::vector<std::string> returnFiles;
  size_t fCallsNEB;
  std::shared_ptr<spdlog::logger> log;
};
