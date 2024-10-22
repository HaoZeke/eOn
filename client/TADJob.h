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
#pragma once

#include "Job.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
namespace eonc {
class TADJob : public Job {
public:
  TADJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~TADJob() = default;
  std::vector<std::string> run(void);

private:
  int dynamics();
  long refine(std::vector<std::shared_ptr<Matter>> buff, long length,
              Matter *reactant);
  bool checkState(Matter *current, Matter *reactant);
  void saveData(int status);
  void dephase();
  bool saddleSearch(std::shared_ptr<Matter> cross);

  std::shared_ptr<Matter> current;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;
  std::shared_ptr<Matter> crossing;
  std::shared_ptr<Matter> final_state;
  std::shared_ptr<Matter> final_tmp;
  std::shared_ptr<Matter> product;
  std::shared_ptr<spdlog::logger> log;

  bool metaStateFlag;
  bool newStateFlag;

  long minimizeFCalls;
  long mdFCalls;
  long dephaseFCalls;
  long refineFCalls;

  long transitionStep;

  double time;
  double minCorrectedTime;
  double transitionTime;
  double barrier;
  double *timeBuffer;
  MinModeSaddleSearch *dimerSearch;

  std::vector<std::string> returnFiles;
};

} // namespace eonc
