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
#include "HessianJob.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"
namespace eonc {
std::vector<std::string> HessianJob::run(void) {
  std::string matter_in("pos.con");

  std::vector<std::string> returnFiles;

  auto matter = std::make_unique<Matter>(pot, params);

  matter->con2matter(matter_in);

  Hessian hessian(params.get(), matter.get());
  long nAtoms = matter->numberOfAtoms();

  Vector<int> moved(nAtoms);
  moved.setConstant(-1);

  int nMoved = 0;
  for (int i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      moved[nMoved] = i;
      nMoved++;
    }
  }
  moved = moved.head(nMoved);
  hessian.getFreqs(matter.get(), moved);

  FILE *fileResults;
  //    FILE *fileMode;

  std::string results_file("results.dat");

  returnFiles.push_back(results_file);

  fileResults = fopen(results_file.c_str(), "wb");

  // fprintf(fileResults, "%d force_calls\n", Potential::fcalls);
  fclose(fileResults);

  return returnFiles;
}

} // namespace eonc
