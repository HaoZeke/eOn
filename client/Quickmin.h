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

#include "Optimizer.h"
#include "Parameters.h"

// when changing away from final, remember to mark the destructor as virtual
class Quickmin final : public Optimizer {

public:
  Quickmin(std::shared_ptr<ObjectiveFunction> a_objf,
           std::shared_ptr<Parameters> a_params)
      : Optimizer(a_objf, OptType::QM, a_params),
        m_dt{a_params->optim.timeStep},
        m_dt_max{a_params->optim.maxTimeStep},
        m_max_move{a_params->optim.maxMove},
        m_vel{VectorType::Zero(a_objf->degreesOfFreedom())},
        m_iteration{0},
        m_max_iter{a_params->optim.maxIterations} {
    if (spdlog::get("qm")) {
      m_log = spdlog::get("qm");
    } else {
      m_log = spdlog::basic_logger_mt("qm", "_qm.log", true);
    }
    m_log->set_pattern("[%l] [QM] %v");
  }
  ~Quickmin() = default;

  int step(double a_maxMove) override;
  int run(size_t a_maxIterations, double a_maxMove) override;

private:
  double m_dt, m_dt_max, m_max_move;
  VectorType m_vel;
  size_t m_iteration, m_max_iter;
  std::shared_ptr<spdlog::logger> m_log;
};
