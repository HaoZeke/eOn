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

#ifndef IMD_POTENTIAL
#define IMD_POTENTIAL

#include "../../Potential.h"

class IMD : public Potential {

public:
  IMD(Parameters& a_p)
      : Potential(a_p){};
  ~IMD();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void writeConfIMD(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void readForceIMD(long N, double *F, double *U);
  //        void spawnVASP();
  //        bool vaspRunning();
  //        static bool firstRun;
  //        static long vaspRunCount;
  //        static pid_t vaspPID;
};

#endif
