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

#include "EDIP.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void EDIP::force(long N, const double *R, const int * /*atomicNrs*/, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  // edipFortran.f90 OpenMP energy reduction races: shared ener is not zeroed
  // under a barrier before partial sums, so OMP_NUM_THREADS>1 yields wrong
  // energies while forces can still look plausible. EDIP is not
  // shared-instance thread-safe; evaluate serially.
#ifdef _OPENMP
  const int prev_threads = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  m_force(&N, R, F, U, &box[0], &box[4], &box[8]);
#ifdef _OPENMP
  omp_set_num_threads(prev_threads);
#endif
}
