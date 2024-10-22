/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/

#ifndef __MPI_POTENTIAL__
#define __MPI_POTENTIAL__

#include "../../Parameters.h"
#include "../../Potential.h"

class MPIPot : public Potential {

public:
  MPIPot(std::shared_ptr<Parameters> p);
  ~MPIPot();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  int potentialRank;
  double poll_period;
};

#endif
