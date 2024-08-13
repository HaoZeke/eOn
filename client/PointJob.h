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

#include "client/matter/Matter.h"
namespace eonc {
class PointJob {
  Matter &_mat;

public:
  PointJob(Matter &mat_a)
      : _mat{mat_a} {
    log = spdlog::get("combi");
  }
  ~PointJob(void) = default;
  bool runImpl(void);

private:
  std::shared_ptr<spdlog::logger> log;
};

} // namespace eonc
