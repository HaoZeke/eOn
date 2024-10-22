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
#include "Matter.h"
#include "SaddleSearchMethod.h"

class BiasedGradientSquaredDescent : public SaddleSearchMethod {
public:
  BiasedGradientSquaredDescent(std::shared_ptr<Matter> matterPassed,
                               double reactantEnergyPassed,
                               std::shared_ptr<Parameters> parametersPassed)
      : SaddleSearchMethod(matterPassed->getPotential(), parametersPassed),
        saddle{matterPassed} {
    reactantEnergy = reactantEnergyPassed;
    saddle = matterPassed;
    eigenvector.resize(saddle->numberOfAtoms(), 3);
    eigenvector.setZero();
    log = spdlog::get("combi");
  }
  ~BiasedGradientSquaredDescent() = default;

  int run(void);
  double getEigenvalue();
  AtomMatrix getEigenvector();

  double eigenvalue;
  AtomMatrix eigenvector;

  std::shared_ptr<Matter> saddle;

  int status;

private:
  double reactantEnergy;
  std::shared_ptr<spdlog::logger> log;
  //        double bgsdAlpha;
};
