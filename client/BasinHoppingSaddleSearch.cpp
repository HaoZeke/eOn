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
#include "eon/BasinHoppingSaddleSearch.h"
#include "eon/HelperFunctions.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/NudgedElasticBand.h"
#include <cmath>
#include <stdexcept>

int BasinHoppingSaddleSearch::run() {
  // minimize "saddle"
  saddle->relax(false, true, false, "displacementmin");
  product = std::make_shared<Matter>(pot, params);
  *product = *saddle;
  // accept or reject based on boltzman
  // exp(-de/(kB*params.main_options.temperature))
  double eproduct, ereactant, de;
  eproduct = product->getPotentialEnergy();
  ereactant = reactant->getPotentialEnergy();
  de = eproduct - ereactant;
  double kB = params.constants.kB;
  double Temperature = params.main_options.temperature;
  if (ereactant < eproduct) {
    if (!(Temperature > 0.0) || !(kB > 0.0)) {
      return 1;
    }
    const double arg = -de / (kB * Temperature);
    const double p = std::exp(arg);
    if (eonc::helpers::random() > p) {
      return 1;
    }
  }
  // NEB reactant to minimized "saddle"
  NudgedElasticBand neb(reactant, product, params, pot);
  if (!eonc::io::io_ok(
          neb.path[0]->matter2con("neb_initial_band.con", false))) {
    QUILL_LOG_WARNING(log, "Failed to write neb_initial_band.con");
  }
  for (int j = 1; j < neb.numImages; j++) {
    if (!eonc::io::io_ok(neb.path[j]->matter2con("neb_initial_band", true))) {
      QUILL_LOG_WARNING(log, "Failed to append neb_initial_band frame");
    }
  }
  neb.compute();
  // pick the maximum energy image along the band
  double Emax = -1e100;
  int HighestImage = 0;

  for (int i = 1; i <= neb.numImages; i++) {
    double Etest = neb.path[i]->getPotentialEnergy();
    QUILL_LOG_DEBUG(log, "i: {} Etest: {:.1f}", i, Etest);
    if (Etest > Emax) {
      Emax = Etest;
      HighestImage = i;
    }
  }
  AtomMatrix r_1 = neb.path[HighestImage - 1]->getPositions();
  AtomMatrix r_3 = neb.path[HighestImage + 1]->getPositions();
  AtomMatrix direction = neb.path[HighestImage]->pbc(r_3 - r_1);
  const double dirNorm = direction.norm();
  if (!(dirNorm > 0.0)) {
    throw std::runtime_error(
        "BasinHoppingSaddleSearch: zero NEB tangent at the highest image");
  }
  direction /= dirNorm;
  MinModeSaddleSearch dim(neb.path[HighestImage], direction, ereactant, params,
                          pot);
  dim.run();
  *saddle = *neb.path[HighestImage];
  eigenvalue = dim.getEigenvalue();
  eigenvector = dim.getEigenvector();
  return 0;
}

double BasinHoppingSaddleSearch::getEigenvalue() { return eigenvalue; }

AtomMatrix BasinHoppingSaddleSearch::getEigenvector() { return eigenvector; }
