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
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "EonLogger.h"
#include "HelperFunctions.h"
#include "OHTSTJob.h"

namespace eonc {

// Johannesson-Jonsson OH-TST (JCP 115, 9644 (2001)). The dividing
// surface is the hyperplane n.(x - Gamma(s)) = 0 in the 3N_free
// configuration space, with Gamma(s) = R + s*u on the straight
// reactant->product guideline u = (P - R)/|P - R|. Free energy changes
// are accumulated as reversible work while damped Verlet dynamics on
// (s, n) climb to the maximum-free-energy plane.

double OHTSTJob::uniformDraw() {
  // Deterministic LCG (numerical recipes ranqd1 flavour): the job must
  // be reproducible for a given random_seed across MPI farm workers.
  m_seedState = (1664525L * m_seedState + 1013904223L) & 0x7fffffffL;
  return static_cast<double>(m_seedState) / 2147483648.0;
}

double OHTSTJob::gaussDraw() {
  if (m_gaussHave) {
    m_gaussHave = false;
    return m_gaussSpare;
  }
  double u1 = std::max(uniformDraw(), 1e-12);
  double u2 = uniformDraw();
  double r = std::sqrt(-2.0 * std::log(u1));
  m_gaussSpare = r * std::sin(2.0 * helpers::pi * u2);
  m_gaussHave = true;
  return r * std::cos(2.0 * helpers::pi * u2);
}

void OHTSTJob::drawThermalVelocities(VectorXd &vel, const VectorXd &masses3N,
                                     const VectorXd *normal) {
  for (long k = 0; k < vel.size(); ++k) {
    vel[k] = std::sqrt(m_kbt / masses3N[k]) * gaussDraw();
  }
  if (normal != nullptr) {
    vel -= (*normal) * normal->dot(vel);
  }
}

void OHTSTJob::samplePlane(Matter &matter, const VectorXd &gamma,
                           const VectorXd &normal, double &avgFn,
                           VectorXd &avgRot) {
  const long equilSteps = params.oh_tst_options.equil_steps;
  const long sampleSteps = params.oh_tst_options.sample_steps;

  // Constrain the current geometry exactly onto the plane.
  VectorXd x = matter.getPositionsFreeV();
  x -= normal * normal.dot(x - gamma);
  matter.setPositionsFreeV(x);

  VectorXd v(x.size());
  drawThermalVelocities(v, m_masses3N, &normal);

  VectorXd f = matter.getForcesFreeV();
  VectorXd fPlane = f - normal * normal.dot(f);

  avgFn = 0.0;
  avgRot = VectorXd::Zero(x.size());
  long nAccum = 0;

  for (long step = 0; step < equilSteps + sampleSteps; ++step) {
    // Velocity Verlet on the plane: forces and velocities projected,
    // positions corrected back onto the constraint (RATTLE for a
    // linear constraint is an exact projection).
    VectorXd a = fPlane.cwiseQuotient(m_masses3N);
    x += m_dt * v + 0.5 * m_dt * m_dt * a;
    x -= normal * normal.dot(x - gamma);
    matter.setPositionsFreeV(x);
    f = matter.getForcesFreeV();
    fPlane = f - normal * normal.dot(f);
    VectorXd aNew = fPlane.cwiseQuotient(m_masses3N);
    v += 0.5 * m_dt * (a + aNew);
    v -= normal * normal.dot(v);
    // Andersen collisions keep the constrained ensemble canonical.
    if (uniformDraw() < m_andersenProb) {
      drawThermalVelocities(v, m_masses3N, &normal);
    }
    if (step >= equilSteps) {
      const double fn = normal.dot(f);
      avgFn += fn;
      avgRot.noalias() += fn * (x - gamma);
      ++nAccum;
    }
  }
  if (nAccum > 0) {
    avgFn /= static_cast<double>(nAccum);
    avgRot /= static_cast<double>(nAccum);
  }
}

double OHTSTJob::reactantResidenceFraction(Matter &matter,
                                           const VectorXd &gammaR,
                                           const VectorXd &normal,
                                           double slabWidth) {
  const long steps = params.oh_tst_options.reactant_md_steps;
  const long equilSteps = params.oh_tst_options.equil_steps;
  VectorXd x = matter.getPositionsFreeV();
  VectorXd v(x.size());
  drawThermalVelocities(v, m_masses3N, nullptr);
  VectorXd f = matter.getForcesFreeV();
  long inSlab = 0, counted = 0;
  for (long step = 0; step < equilSteps + steps; ++step) {
    VectorXd a = f.cwiseQuotient(m_masses3N);
    x += m_dt * v + 0.5 * m_dt * m_dt * a;
    matter.setPositionsFreeV(x);
    f = matter.getForcesFreeV();
    VectorXd aNew = f.cwiseQuotient(m_masses3N);
    v += 0.5 * m_dt * (a + aNew);
    if (uniformDraw() < m_andersenProb) {
      drawThermalVelocities(v, m_masses3N, nullptr);
    }
    if (step >= equilSteps) {
      if (std::fabs(normal.dot(x - gammaR)) < 0.5 * slabWidth) {
        ++inSlab;
      }
      ++counted;
    }
  }
  return counted > 0 ? static_cast<double>(inSlab) / counted : 0.0;
}

std::vector<std::string> OHTSTJob::run(void) {
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  if (!reactant->con2matter(params.oh_tst_options.reactant_filename)) {
    EONC_LOG_CRITICAL("OH-TST failed to load {}",
                      params.oh_tst_options.reactant_filename);
    throw std::runtime_error("oh_tst: failed to load reactant");
  }
  if (!product->con2matter(params.oh_tst_options.product_filename)) {
    EONC_LOG_CRITICAL("OH-TST failed to load {}",
                      params.oh_tst_options.product_filename);
    throw std::runtime_error("oh_tst: failed to load product");
  }

  const double temperature = params.main_options.temperature;
  m_kbt = params.constants.kB * temperature;
  m_dt = params.oh_tst_options.time_step / params.constants.timeUnit;
  m_seedState = (params.main_options.randomSeed > 0)
                    ? params.main_options.randomSeed
                    : 12345;
  // Per-step collision probability from the Andersen collision period.
  const double tcol =
      params.thermostat_options.andersen_tcol_input / params.constants.timeUnit;
  m_andersenProb = (tcol > 0.0) ? std::min(1.0, m_dt / tcol) : 0.1;

  // Free-DOF mass vector (amu per coordinate).
  const long nAtoms = reactant->numberOfAtoms();
  auto masses = reactant->getMasses();
  std::vector<double> m3;
  m3.reserve(3 * nAtoms);
  for (long i = 0; i < nAtoms; ++i) {
    if (!reactant->getFixed(i)) {
      for (int j = 0; j < 3; ++j)
        m3.push_back(masses[i]);
    }
  }
  m_masses3N = VectorXd::Map(m3.data(), static_cast<long>(m3.size()));

  // Straight guideline in 3N space, minimum-image on the endpoint
  // difference so the plane never strides a cell boundary.
  const VectorXd xR = reactant->getPositionsFreeV();
  VectorXd diff = product->getPositionsFreeV() - xR;
  {
    AtomMatrix d(AtomMatrix::Map(diff.data(), diff.size() / 3, 3));
    d = reactant->pbc(d);
    diff = VectorXd::Map(d.data(), diff.size());
  }
  const double guideLen = diff.norm();
  if (guideLen < 1e-8) {
    throw std::runtime_error("oh_tst: reactant and product coincide");
  }
  const VectorXd u = diff / guideLen;

  // Plane state: progression coordinate s, normal n, their conjugate
  // velocities. Damped (velocity-zeroing) Verlet climbs A(s, n).
  double s = params.oh_tst_options.s_init * guideLen;
  double vS = 0.0;
  VectorXd n = u;
  VectorXd omega = VectorXd::Zero(n.size());
  const double mS = params.oh_tst_options.plane_mass;
  const double alphaRot = params.oh_tst_options.alpha_rot;
  const double dtPlane = params.oh_tst_options.plane_time_step;
  const double dsMax = params.oh_tst_options.ds_max;
  const double fTol = params.oh_tst_options.force_tol;

  Matter walker(*reactant);

  double aTrans = 0.0, aRot = 0.0, aBest = 0.0, sBest = s;
  VectorXd nBest = n;
  double sPrev = s, fnPrev = 0.0;
  bool havePrev = false;

  FILE *prog = fopen("oh_tst_progression.dat", "w");
  if (prog) {
    fprintf(prog, "# plane  s/L  <F.n> (eV/A)  dA_trans (eV)  dA_rot (eV)  "
                  "A (eV)  n.u\n");
  }

  const long nPlanes = params.oh_tst_options.max_planes;
  long plane = 0;
  bool converged = false;
  for (; plane < nPlanes; ++plane) {
    const VectorXd gamma = xR + s * u;
    double avgFn = 0.0;
    VectorXd avgRot;
    samplePlane(walker, gamma, n, avgFn, avgRot);

    // Translation: dA/ds = -<F.n> (n.u); reversible work accumulated
    // over the distance actually moved (Eq 18).
    const double fS = -avgFn * n.dot(u);
    if (havePrev) {
      aTrans += 0.5 * (fS + fnPrev) * (s - sPrev);
    }
    sPrev = s;
    fnPrev = fS;
    havePrev = true;

    // Rotation: dA/dn = -<(F.n)(x - Gamma)> projected onto the plane
    // (Eq 19 accumulates the corresponding work below).
    VectorXd g = -(avgRot - n * n.dot(avgRot)) / alphaRot;

    const double aTotal = aTrans + aRot;
    if (aTotal > aBest) {
      aBest = aTotal;
      sBest = s;
      nBest = n;
    }
    if (prog) {
      fprintf(prog, "%6ld %10.6f %14.6e %12.6f %12.6f %12.6f %10.6f\n", plane,
              s / guideLen, avgFn, aTrans, aRot, aTotal, n.dot(u));
      fflush(prog);
    }
    EONC_LOG_DEBUG("[oh_tst] plane {} s/L {:.4f} <F.n> {:.4e} A {:.4f} eV",
                   plane, s / guideLen, avgFn, aTotal);

    if (std::fabs(fS) < fTol && g.norm() < fTol && plane > 2) {
      converged = true;
      break;
    }

    // Damped Verlet on s (Eqs 5-7): zero the velocity when it opposes
    // the climbing force so the plane settles at the maximum.
    vS += dtPlane * fS / mS;
    if (vS * fS < 0.0)
      vS = 0.0;
    double ds = dtPlane * vS;
    ds = std::clamp(ds, -dsMax, dsMax);
    s = std::clamp(s + ds, 0.0, guideLen);

    // Damped rotation of the normal (Eqs 8-10): angular velocity in
    // the tangent space of the unit sphere, zeroed on opposing torque,
    // followed by renormalization; the rotational work advances A_rot.
    omega += dtPlane * g;
    omega -= n * n.dot(omega);
    if (omega.dot(g) < 0.0)
      omega.setZero();
    const VectorXd nOld = n;
    n = (n + dtPlane * omega).normalized();
    aRot += g.dot(n - nOld) * alphaRot;
  }
  if (prog)
    fclose(prog);

  // Direction-dependent effective mass (Eq 24) and the one-sided
  // thermal flux factor sqrt(kBT / 2 pi mu).
  const double mu = (m_masses3N.array() * nBest.array().square()).sum();
  const double vFlux = std::sqrt(m_kbt / (2.0 * helpers::pi * mu));

  // Configurational partition ratio by slice residence on an
  // unconstrained reactant trajectory (Eq 22): Q^ZR/Q^R = f / delta.
  Matter rWalker(*reactant);
  const double slab = params.oh_tst_options.slab_width;
  const double fRes =
      reactantResidenceFraction(rWalker, xR, nBest, slab);
  const double qRatio = (slab > 0.0) ? fRes / slab : 0.0;

  // Rate in internal units (1/internal-time), then SI.
  const double kInternal = vFlux * qRatio * std::exp(-aBest / m_kbt);
  const double kSI = kInternal / (params.constants.timeUnit * 1.0e-15);

  std::vector<std::string> returnFiles;
  FILE *out = fopen("results.dat", "w");
  if (out) {
    fprintf(out, "%s job_type\n", "oh_tst");
    fprintf(out, "%d converged\n", converged ? 1 : 0);
    fprintf(out, "%ld planes_used\n", plane);
    fprintf(out, "%.8f free_energy_barrier_eV\n", aBest);
    fprintf(out, "%.8f s_star_over_L\n", sBest / guideLen);
    fprintf(out, "%.8f normal_overlap_with_guideline\n", nBest.dot(u));
    fprintf(out, "%.8e effective_mass_amu\n", mu);
    fprintf(out, "%.8e residence_fraction\n", fRes);
    fprintf(out, "%.8e q_ratio_per_A\n", qRatio);
    fprintf(out, "%.8e rate_ohtst_per_s\n", kSI);
    fprintf(out, "%.4f temperature_K\n", temperature);
    fclose(out);
    returnFiles.push_back("results.dat");
  }
  returnFiles.push_back("oh_tst_progression.dat");
  EONC_LOG_INFO("[oh_tst] {} after {} planes: A = {:.4f} eV at s/L = {:.4f}, "
                "k = {:.4e} 1/s at {:.1f} K",
                converged ? "converged" : "max planes", plane, aBest,
                sBest / guideLen, kSI, temperature);
  return returnFiles;
}

} // namespace eonc
