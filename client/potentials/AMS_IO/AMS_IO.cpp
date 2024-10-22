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

#include "AMS_IO.h"
#include <string>
#include <unistd.h>

AMS_IO::AMS_IO(AMSParams amsp)
    : Potential(PotType::AMS_IO) {
  engine = amsp.engine.c_str();
  forcefield = amsp.forcefield.c_str();
  model = amsp.model.c_str();
  xc = amsp.xc.c_str();
  return;
}

void AMS_IO::cleanMemory(void) { return; }

AMS_IO::~AMS_IO() { cleanMemory(); }

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg",      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr",      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd",      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd",      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf",      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po",      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  NULL};

// guess the atom type from the atomic mass,
std::string mass2atom(double atomicmass) {
  return elementArray[int(atomicmass + .5)];
}

int symbol2atomicNumber(char const *symbol) {
  int i = 0;

  while (elementArray[i] != NULL) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  // invalid symbol
  return -1;
}

char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

void AMS_IO::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  passToSystem(N, R, atomicNrs, box);
  system("chmod +x run_AMS_IO.sh");
  system(
      "./run_AMS_IO.sh >> ams_output"); // Run a single point AMS_IO calculation
                                        // and write the results into ams_output
  recieveFromSystem(N, F, U);
  return;
}

void AMS_IO::passToSystem(long N, const double *R, const int *atomicNrs,
                          const double *box)
// Creating the standard input file that will be read by the AMS_IO driver
{
  FILE *out;
  out = fopen("run_AMS_IO.sh", "w");

  fprintf(out, "#!/bin/sh\n");
  fprintf(out, "ams --delete-old-results <<eor\n");
  fprintf(out, "Task SinglePoint\n");
  fprintf(out, "System\n");
  fprintf(out, " Atoms\n");
  for (int i = 0; i < N; i++) {
    fprintf(out, "  %s\t%.19f\t%.19f\t%.19f\n",
            atomicNumber2symbol(atomicNrs[i]), R[i * 3 + 0], R[i * 3 + 1],
            R[i * 3 + 2]);
  }
  fprintf(out, " End\n");
  if (strlen(model) > 0 || strlen(forcefield)) {
    fprintf(out, " Lattice\n");
    for (int i = 0; i < 3; i++) {
      fprintf(out, "  %.19f\t%.19f\t%.19f\n", box[i * 3 + 0], box[i * 3 + 1],
              box[i * 3 + 2]);
    }
    fprintf(out, " End\n");
  }
  fprintf(out, "End\n");
  fprintf(out, "Engine %s\n", engine);
  if (strlen(forcefield) > 0) {
    fprintf(out, "     Forcefield %s\n", forcefield);
  }
  if (strlen(model) > 0) {
    fprintf(out, "     Model %s\n", model);
  }
  if (strlen(xc) > 0) {
    fprintf(out, "xc %s\n");
    fprintf(out, "     hybrid %s\n",
            xc); // basis set not specified (default = DZ)
    fprintf(out, "end\n");
  }
  fprintf(out, "EndEngine\n");
  fprintf(out, "Properties\n");
  fprintf(out, " Gradients\n");
  fprintf(out, "End\n");
  fprintf(out, "eor");
  fclose(out);
  return;
}

void AMS_IO::recieveFromSystem(long N, double *F, double *U) {

  FILE *in;
  double junkF;
  char junkChar[256];
  double forceX;
  double forceY;
  double forceZ;
  double index;
  char line[256];

  in = fopen("ams_output", "r");

  while (fgets(line, sizeof(line), in)) {

    if (strcmp(line, "     CALCULATION RESULTS\n") ==
        0) { // Finding the Energy in the output file

      fscanf(in, "%s %s %s %lf", junkChar, junkChar, junkChar, U);
      *U = *U * 27.2114; // Energy in hartree to eV
    }
    if (strcmp(line, "  Index   Atom            d/dx            d/dy           "
                     " d/dz\n") == 0) { // Finding the forces
      for (int i = 0; i < N; i++) {
        fscanf(in, "%lf %s %lf %lf %lf", &index, &junkChar, &forceX, &forceY,
               &forceZ);
        F[int(i) * 3 + 0] = -forceX;
        F[int(i) * 3 + 1] = -forceY;
        F[int(i) * 3 + 2] =
            -forceZ; // AMS_IO gives gradients, not forces, hence the change.
      }
    }
  }
  for (int i = 0; i < 3 * N; i++) {
    F[i] = F[i] * 51.4220862; // Forces from hartree/bohr to eV/Angstrom
  }

  fclose(in);
  return;
}
