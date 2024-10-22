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
#include "Matter.h"
#include "Parameters.h"

#include "Eigen.h"
namespace eonc {
/** Functionality relying on the conjugate gradients algorithm. The object is
 * capable of minimizing an Matter object or modified forces being passed in.*/
class BondBoost {

public:
  /** Constructor to be used when a structure is minimized.
  @param[in]   *matter        Pointer to the Matter object to be relaxed.
  @param[in]   *parameters    Pointer to the Parameter object containing the
  runtime parameters.*/
  BondBoost(Matter *matt, Parameters *params);
  ~BondBoost(); ///< Destructor.

  void initialize();
  double boost();

private:
  VectorType Rmdsteps();
  long BondSelect();
  double Booststeps();
  long nAtoms;    ///< Number of free coordinates.
  Matter *matter; ///< Pointer to atom object \b outside the scope of the class.
  Parameters *parameters; ///< Pointer to a structure outside the scope of the
                          ///< class containing runtime parameters.
  long *BAList;
  long *RAList;
  long *TABAList;
  long *BBAList;
  double *Epsr_Q;
  VectorType TABLList; // EquilibriumTaggedAtomInvolvedBondLengthList;
  VectorType EBBLList; // EquilibriumBoostBondLengthList
  VectorType CBBLList; // CurrentBoostBondLengthList
  long nBAs;
  long nRAs;
  long nTABs;
  long nReg;
  long nBBs;
  std::shared_ptr<spdlog::logger> log;
};

class Hyperdynamics {
public:
  //            static const string NONE;
  //            static const string BOND_BOOST;
  static const char NONE[];
  static const char BOND_BOOST[];
};

} // namespace eonc
