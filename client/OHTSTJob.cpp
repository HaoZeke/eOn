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
// reactant->product guideline u = (P - R)/|P - R|. Sign conventions
// follow Appendix A: the plane translates along the driving force
// -<F.n> and its normal is driven along +<(n.F) R/(alpha |R|^2)>
// (verified to tilt the normal onto the ridge normal for a tilted-
// ridge model), so both damped Verlet integrations climb the free
// energy and settle where the mean forces vanish.

double OHTSTJob::uniformDraw() {
  // Deterministic LCG: the job must be reproducible for a given
  // random_seed across MPI farm workers.
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

void OHTSTJob::drawThermalVelocities(VectorXd &vel, const VectorXd *normal) {
  for (long k = 0; k < vel.size(); ++k) {
    vel[k] = std::sqrt(m_kbt / m_masses3N[k]) * gaussDraw();
  }
  if (normal != nullptr) {
    vel -= (*normal) * normal->dot(vel);
  }
}

OHTSTJob::PlaneAverages OHTSTJob::samplePlane(Matter &matter,
                                              const VectorXd &gamma,
                                              const VectorXd &normal) {
  const long equilSteps = params.oh_tst_options.equil_steps;
  const long sampleSteps = params.oh_tst_options.sample_steps;
  const double alphaRot = params.oh_tst_options.alpha_rot;

  // Constrain the current geometry exactly onto the plane.
  VectorXd x = matter.getPositionsFreeV();
  x -= normal * normal.dot(x - gamma);
  matter.setPositionsFreeV(x);

  VectorXd v(x.size());
  drawThermalVelocities(v, &normal);

  VectorXd f = matter.getForcesFreeV();
  VectorXd fPlane = f - normal * normal.dot(f);

  PlaneAverages avg;
  avg.rotNorm = VectorXd::Zero(x.size());
  avg.rotRaw = VectorXd::Zero(x.size());
  avg.pos = VectorXd::Zero(x.size());
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
      drawThermalVelocities(v, &normal);
    }
    if (step >= equilSteps) {
      const double fn = normal.dot(f);
      const VectorXd arm = x - gamma;
      const double arm2 = arm.squaredNorm();
      avg.fn += fn;
      if (arm2 > 1e-16) {
        avg.rotNorm.noalias() += (fn / (alphaRot * arm2)) * arm;
      }
      avg.rotRaw.noalias() += fn * arm;
      avg.pos.noalias() += x;
      ++nAccum;
    }
  }
  if (nAccum > 0) {
    const double inv = 1.0 / static_cast<double>(nAccum);
    avg.fn *= inv;
    avg.rotNorm *= inv;
    avg.rotRaw *= inv;
    avg.pos *= inv;
  }
  return avg;
}

double OHTSTJob::reactantQRatio(Matter &matter, const VectorXd &gammaR,
                                const VectorXd &normal) {
  // Eq 22: Q^ZR/Q^R = (dt / t_tot) * sum over plane crossings of
  // 1 / |(r_{i+1} - r_i).n|. The trajectory is unconstrained,
  // thermostatted, and stays in the reactant basin by construction
  // (it starts there and the barrier is >> kT).
  const long steps = params.oh_tst_options.reactant_md_steps;
  const long equilSteps = params.oh_tst_options.equil_steps;
  VectorXd x = matter.getPositionsFreeV();
  VectorXd v(x.size());
  drawThermalVelocities(v, nullptr);
  VectorXd f = matter.getForcesFreeV();
  double side = normal.dot(x - gammaR);
  double crossingSum = 0.0;
  long counted = 0;
  for (long step = 0; step < equilSteps + steps; ++step) {
    const VectorXd xOld = x;
    VectorXd a = f.cwiseQuotient(m_masses3N);
    x += m_dt * v + 0.5 * m_dt * m_dt * a;
    matter.setPositionsFreeV(x);
    f = matter.getForcesFreeV();
    VectorXd aNew = f.cwiseQuotient(m_masses3N);
    v += 0.5 * m_dt * (a + aNew);
    if (uniformDraw() < m_andersenProb) {
      drawThermalVelocities(v, nullptr);
    }
    const double sideNew = normal.dot(x - gammaR);
    if (step >= equilSteps) {
      if (side * sideNew < 0.0) {
        const double proj = std::fabs(normal.dot(x - xOld));
        if (proj > 1e-14) {
          crossingSum += 1.0 / proj;
        }
      }
      ++counted;
    }
    side = sideNew;
  }
  return counted > 0 ? crossingSum / static_cast<double>(counted) : 0.0;
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
  // velocities, and the previous-iteration driving forces for the
  // two-force velocity Verlet updates (Eqs 6-7 and 9-10).
  double s = params.oh_tst_options.s_init * guideLen;
  double vS = 0.0;
  VectorXd n = u;
  VectorXd omega = VectorXd::Zero(n.size());
  const double mS = params.oh_tst_options.plane_mass;
  const double dtPlane = params.oh_tst_options.plane_time_step;
  const double dsMax = params.oh_tst_options.ds_max;
  const double dThetaMax = params.oh_tst_options.dtheta_max;
  const double fTol = params.oh_tst_options.force_tol;

  Matter walker(*reactant);

  // Reversible-work accumulators and the previous plane's averages
  // for the trapezoid rules of Eqs 18-19.
  double aTrans = 0.0, aRot = 0.0, aBest = 0.0, sBest = s;
  VectorXd nBest = n;
  bool havePrev = false;
  double fnPrev = 0.0;
  VectorXd rotRawPrev, posPrev, nPrev;
  double gSPrev = 0.0;
  VectorXd gRotPrev;
  int sideSign = 0; // sign of <F.n> at plane 1 (grooming, Sec IIE)

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
    PlaneAverages avg = samplePlane(walker, gamma, n);

    // Driving forces: translation climbs against <F.n> (Eq 5 with the
    // reversed-force convention); the normal is driven along
    // +<(n.F) R/(alpha |R|^2)> (Appendix A, Eq A1), projected onto
    // the tangent space of the unit sphere.
    const double gS = -avg.fn / mS;
    VectorXd gRot = avg.rotNorm - n * n.dot(avg.rotNorm);

    if (plane == 0) {
      sideSign = (avg.fn < 0.0) ? -1 : 1;
    }

    // Reversible work (Eqs 18-19), trapezoid between consecutive
    // planes: the translational path is the piecewise-linear track of
    // the average configuration <r>, and the rotational work pairs
    // the unnormalized <(n.F) R> with the actual change of the
    // normal. Grooming (Sec IIE): only planes on the reactant side of
    // the ridge (same sign of <F.n> as plane 1) contribute.
    if (havePrev) {
      const bool sameSide = (fnPrev < 0.0 ? -1 : 1) == sideSign &&
                            (avg.fn < 0.0 ? -1 : 1) == sideSign;
      if (sameSide) {
        const VectorXd fParMean = 0.5 * (fnPrev * nPrev + avg.fn * n);
        aTrans += -fParMean.dot(avg.pos - posPrev);
        const VectorXd rotMean = 0.5 * (rotRawPrev + avg.rotRaw);
        aRot += rotMean.dot(n - nPrev);
      }
    }

    const double aTotal = aTrans + aRot;
    if (aTotal > aBest) {
      aBest = aTotal;
      sBest = s;
      nBest = n;
    }
    if (prog) {
      fprintf(prog, "%6ld %10.6f %14.6e %12.6f %12.6f %12.6f %10.6f\n", plane,
              s / guideLen, avg.fn, aTrans, aRot, aTotal, n.dot(u));
      fflush(prog);
    }
    EONC_LOG_DEBUG("[oh_tst] plane {} s/L {:.4f} <F.n> {:.4e} A {:.4f} eV",
                   plane, s / guideLen, avg.fn, aTotal);

    if (plane > 2 && std::fabs(avg.fn) < fTol &&
        gRot.norm() * params.oh_tst_options.alpha_rot < fTol) {
      converged = true;
      // The converged plane is the optimal one even if sampling noise
      // put an earlier plane marginally higher.
      aBest = aTotal;
      sBest = s;
      nBest = n;
      if (prog) {
        fprintf(prog, "# converged at plane %ld\n", plane);
      }
      break;
    }

    // Damped two-force Verlet on s (Eqs 6-7): the velocity is zeroed
    // when it opposes the driving force, so the plane settles at the
    // free-energy maximum instead of oscillating over it ("the
    // velocity of the plane is zeroed if the plane has gone past the
    // maximum free energy position").
    if (havePrev) {
      vS += 0.5 * dtPlane * (gS + gSPrev);
    } else {
      vS += dtPlane * gS;
    }
    if (vS * gS < 0.0)
      vS = 0.0;
    double ds = dtPlane * vS + 0.5 * dtPlane * dtPlane * gS;
    ds = std::clamp(ds, -dsMax, dsMax);
    const double sNew = std::clamp(s + ds, 0.0, guideLen);

    // Damped rotation (Eqs 9-10 with the Appendix A driving): the
    // angular velocity keeps only its projection along the current
    // driving force while aligned with it, and is zeroed otherwise.
    if (havePrev && gRotPrev.size() == gRot.size()) {
      omega += 0.5 * dtPlane * (gRot + gRotPrev);
    } else {
      omega += dtPlane * gRot;
    }
    omega -= n * n.dot(omega);
    const double gNorm = gRot.norm();
    if (gNorm > 1e-14) {
      const VectorXd gHat = gRot / gNorm;
      const double along = omega.dot(gHat);
      if (along > 0.0) {
        omega = gHat * along;
      } else {
        omega.setZero();
      }
    } else {
      omega.setZero();
    }
    VectorXd dn = dtPlane * omega + 0.5 * dtPlane * dtPlane * gRot;
    const double dTheta = dn.norm();
    if (dTheta > dThetaMax) {
      dn *= dThetaMax / dTheta;
    }
    const VectorXd nOld = n;
    n = (n + dn).normalized();
    if (n.dot(u) < 0.0) {
      // Keep the normal pointing towards increasing s.
      n = -n;
    }

    // Eq 11-12 restart: place the walker in the new plane at the
    // rotated image of the average arm, rescaled so rotation does not
    // change the arm length.
    VectorXd arm = avg.pos - gamma;
    VectorXd armNew = arm - nOld * arm.dot(n - nOld);
    const double armLen = arm.norm();
    const double armNewLen = armNew.norm();
    if (armLen > 1e-12 && armNewLen > 1e-12) {
      armNew *= armLen / armNewLen;
    }
    const VectorXd gammaNew = xR + sNew * u;
    VectorXd xStart = gammaNew + armNew;
    xStart -= n * n.dot(xStart - gammaNew);
    walker.setPositionsFreeV(xStart);

    fnPrev = avg.fn;
    nPrev = nOld;
    rotRawPrev = avg.rotRaw;
    posPrev = avg.pos;
    gSPrev = gS;
    gRotPrev = gRot;
    havePrev = true;
    s = sNew;
  }
  if (prog)
    fclose(prog);

  // Direction-dependent effective mass (Eq 24) and the one-sided
  // thermal flux factor sqrt(kBT / 2 pi mu) (Eq 23).
  const double mu = (m_masses3N.array() * nBest.array().square()).sum();
  const double vFlux = std::sqrt(m_kbt / (2.0 * helpers::pi * mu));

  // Q^ZR/Q^R at the FIRST plane of the progression (Z^R), whose
  // reversible work to the optimal plane is what aBest measures.
  Matter rWalker(*reactant);
  const VectorXd gammaR = xR + (params.oh_tst_options.s_init * guideLen) * u;
  const double qRatio = reactantQRatio(rWalker, gammaR, u);

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
    fprintf(out, "%.8f delta_a_trans_eV\n", aTrans);
    fprintf(out, "%.8f delta_a_rot_eV\n", aRot);
    fprintf(out, "%.8f s_star_over_L\n", sBest / guideLen);
    fprintf(out, "%.8f normal_overlap_with_guideline\n", nBest.dot(u));
    fprintf(out, "%.8e effective_mass_amu\n", mu);
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
