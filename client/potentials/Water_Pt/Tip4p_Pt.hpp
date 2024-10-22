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
/** @file
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#ifndef TIP4P_PT
#define TIP4P_PT
#include "../../Potential.h"
#include "zhu_philpott.hpp"

class Tip4p_Pt : public Potential, private forcefields::ZhuPhilpott<> {
public:
  Tip4p_Pt()
      : Potential(PotType::TIP4P_PT),
        // TODO(rg): Expose these
        forcefields::ZhuPhilpott<>(8.5, 1.0){};
  // Functions
  // constructor and destructor

  // To satisfy interface
  void initialize(void) {}
  void cleanMemory(void) {}
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};

#endif
