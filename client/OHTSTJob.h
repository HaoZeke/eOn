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
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

namespace eonc {

/**
 * @file
 * @ingroup Jobs
 *
 * \brief Optimized hyperplanar transition state theory (OH-TST).
 *
 * Implements the reversible-work variational TST of Johannesson and
 * Jonsson, J. Chem. Phys. 115, 9644 (2001), doi:10.1063/1.1415499.
 * A hyperplanar dividing surface in the 3N-dimensional configuration
 * space, parameterized by a point on a reactant->product guideline and
 * a unit normal, is advanced by damped (velocity-zeroing) Verlet
 * dynamics on the plane position (Eqs 5-7) and orientation (Eqs 8-10),
 * driven by canonical-ensemble averages of the force sampled with
 * thermostatted dynamics constrained to the plane. The free-energy
 * profile follows from the reversible work (Eqs 18-19), the
 * direction-dependent effective mass enters the flux prefactor
 * (Eq 24), and the configurational partition-function ratio is
 * estimated by slice-residence counting on an unconstrained reactant
 * trajectory (Eq 22).
 */
class OHTSTJob : public Job {

public:
  OHTSTJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~OHTSTJob(void) = default;
  std::vector<std::string> run(void);

private:
  //! One thermostatted, plane-constrained sampling block; returns the
  //! averaged normal force <F.n> and the averaged in-plane rotational
  //! force vector <(F.n)(x - gamma)> over the sampling steps.
  void samplePlane(Matter &matter, const VectorXd &gamma,
                   const VectorXd &normal, double &avgFn, VectorXd &avgRot);

  //! Unconstrained reactant-basin trajectory; returns the fraction of
  //! time spent within +/- slabWidth/2 of the reactant plane (Eq 22).
  double reactantResidenceFraction(Matter &matter, const VectorXd &gammaR,
                                   const VectorXd &normal, double slabWidth);

  //! Maxwell-Boltzmann draw on the free DOF, then projection onto the
  //! plane (n.v = 0) when a normal is supplied.
  void drawThermalVelocities(VectorXd &vel, const VectorXd &masses3N,
                             const VectorXd *normal);

  VectorXd m_masses3N;     //!< per-DOF masses of the free atoms (amu)
  double m_dt{0.0};        //!< integration step, internal units
  double m_kbt{0.0};       //!< k_B T (eV)
  double m_andersenProb{0.0}; //!< per-step Andersen collision probability
  long m_seedState{12345}; //!< LCG state for the thermostat draws
  double uniformDraw();    //!< uniform (0,1)
  double gaussDraw();      //!< standard normal (Box-Muller)
  bool m_gaussHave{false};
  double m_gaussSpare{0.0};
};

} // namespace eonc

using eonc::OHTSTJob;
