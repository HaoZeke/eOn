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
#include "eon/StructureComparisonJob.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Optimizer.h"
#include <stdexcept>

std::vector<std::string> StructureComparisonJob::run() {
  std::vector<std::string> returnFiles;

  auto matter1 = std::make_unique<Matter>(pot, params);
  if (!eonc::io::io_ok(matter1->con2matter("matter1.con"))) {
    EONC_LOG_CRITICAL("Failed to load matter1.con");
    throw std::runtime_error("failed to load matter1.con");
  }

  return returnFiles;
}
