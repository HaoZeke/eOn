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

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>

class Potential {
protected:
  PotType ptype;
  std::shared_ptr<Parameters> m_params;
  std::shared_ptr<spdlog::logger> m_log;

public:
  size_t forceCallCounter;

  // Main Constructor
  Potential(PotType a_ptype, std::shared_ptr<Parameters> a_params)
      : ptype{a_ptype},
        m_params{a_params},
        forceCallCounter{0} {
    SPDLOG_TRACE("CREATED WITH {}", forceCallCounter);
    initializeLogger();
  }

  // Delegating Constructor
  Potential(std::shared_ptr<Parameters> a_params)
      : Potential(a_params->pot.potential, a_params) {}

  virtual ~Potential() {
    SPDLOG_TRACE("DESTROYED AFTER {}\n", forceCallCounter);
    if (m_log) {
      m_log->trace("[{}] destroyed after {} calls",
                   magic_enum::enum_name<PotType>(getType()), forceCallCounter);
    } else {
      std::cerr << "Logger is not initialized\n";
    }
  }

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;

  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix pos, const Vector<int> atmnrs, const Matrix3S box);

  PotType getType() { return this->ptype; }

  // Logger initialization
  void initializeLogger() {
    if (!spdlog::get("_potcalls")) {
      // Create logger if it doesn't exist
      m_log = spdlog::basic_logger_mt("_potcalls", "_potcalls.log", true);
      m_log->set_pattern("[%l] [%Y-%m-%d %H:%M:%S] %v");
    } else {
      // Use existing logger
      m_log = spdlog::get("_potcalls");
    }
    if (m_log) {
      m_log->trace("[{}] created", magic_enum::enum_name<PotType>(getType()));
    }
  }
};

namespace helper_functions {
std::shared_ptr<Potential> makePotential(std::shared_ptr<Parameters> params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         std::shared_ptr<Parameters> params);
} // namespace helper_functions
