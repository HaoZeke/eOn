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

#ifndef LENNARD_JONES
#define LENNARD_JONES

// #include "../../system_unit.h" // unit converters
#include "../../Potential.h"

/** Lennard Jones potential.*/
class LJ : public Potential {
private:
  double u0;
  double cuttOffR;
  double psi;
  double cuttOffU;

public:
  LJ(double u0, double cuttOffR, double psi)
      : Potential(PotType::LJ),
        u0{u0},
        cuttOffR{cuttOffR},
        psi{psi} {}

  ~LJ() = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
#endif
