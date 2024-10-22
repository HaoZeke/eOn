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
#include "Quickmin.h"
#include "HelperFunctions.h"

int Quickmin::step(double a_maxMove) {
  VectorType force = -m_objf->getGradient();
  if (m_params->optim.QMSteepestDecent) {
    m_vel.setZero();
  } else {
    if (m_vel.dot(force) < 0) {
      m_vel.setZero();
    } else {
      VectorType f_unit = force.normalized();
      m_vel = m_vel.dot(f_unit) * f_unit;
    }
  }

  m_vel += force * m_dt;
  VectorType dr = helper_functions::maxAtomMotionAppliedV(
      m_vel * m_dt, a_maxMove); // used to be m_params->optMaxTimeStep
  SPDLOG_LOGGER_INFO(m_log, "{} M_Vel.norm() is {}", m_iteration, m_vel.norm());
  m_objf->setPositions(m_objf->getPositions() + dr);
  m_iteration++;
  return m_objf->isConverged() ? 1 : 0;
}

int Quickmin::run(size_t a_maxSteps, double a_maxMove) {
  while (!m_objf->isConverged() && m_iteration < a_maxSteps) {
    step(a_maxMove);
  }
  return m_objf->isConverged() ? 1 : 0;
}
