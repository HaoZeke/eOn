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
#include "Dimer.h"
#include "HelperFunctions.h"
namespace eonc {
using namespace helper_functions;

Dimer::Dimer(std::shared_ptr<Matter> matter, std::shared_ptr<Parameters> params,
             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  matterCenter = std::make_shared<Matter>(pot, params);
  matterDimer = std::make_shared<Matter>(pot, params);
  *matterCenter = *matter;
  *matterDimer = *matter;
  nAtoms = matter->numberOfAtoms();

  direction.resize(nAtoms, 3);
  rotationalPlane.resize(nAtoms, 3);
  direction.setZero();
  rotationalPlane.setZero();
  totalForceCalls = 0;
  log = spdlog::basic_logger_mt("dimer", "dimer.log", true);
  log->set_pattern("%v");
}

// was estimateLowestEigenmode. rename to compute
void Dimer::compute(std::shared_ptr<Matter> matter,
                    AtomMatrix initialDirection) {
  long rotations = 0;
  long forceCallsCenter;
  long forceCallsDimer;
  double rotationalForce1;
  double rotationalForce2;
  double curvature, rotationalForceChange, forceDimer, rotationAngle;
  double lengthRotationalForceOld;
  double torque = 0;
  bool doneRotating = false;

  *matterCenter = *matter;
  rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
  rotationalForce1 = 0;
  rotationalForce2 = 0;
  AtomMatrix rotationalForce(nAtoms, 3);
  AtomMatrix rotationalForceOld(nAtoms, 3);
  AtomMatrix rotationalPlaneOld(nAtoms, 3);
  rotationalForce.setZero();
  rotationalForceOld.setZero();
  rotationalPlaneOld.setZero();
  initialDirection.normalize();
  direction = initialDirection;

  statsAngle = 0;
  lengthRotationalForceOld = 0;
  forceCallsCenter = matterCenter->getForceCalls();
  forceCallsDimer = matterDimer->getForceCalls();

  // uses two force calls per rotation
  while (!doneRotating) {
    // calculate the rotational force and curvature
    curvature = calcRotationalForceReturnCurvature(rotationalForce);

    // determine the new rotational plane
    determineRotationalPlane(rotationalForce, rotationalForceOld,
                             rotationalPlaneOld, &lengthRotationalForceOld);

    // calculate the torque on the dimer
    // GH        torque = rotationalForce.squaredNorm();
    torque = rotationalForce.norm();
    assert(std::isnormal(torque));

    // convergence scheme
    if ((torque > params->dimer.torqueMax &&
         rotations >= params->dimer.rotationsMax) ||
        (torque < params->dimer.torqueMax &&
         torque >= params->dimer.torqueMin &&
         rotations >= params->dimer.rotationsMin) ||
        (torque < params->dimer.torqueMin)) {
      /*            cout << "torque: "<<torque<<endl;
                  cout << "params->dimer.torqueMax:
         "<<params->dimer.torqueMax<<endl; cout <<
         "params->dimer.torqueMin: "<<params->dimer.torqueMin<<endl; cout
         << "rotations: "<<rotations<<endl; cout <<
         "params->dimer.rotationsMax: "<<params->dimer.rotationsMax<<endl;
                  cout << "params->dimer.rotationsMin:
         "<<params->dimer.rotationsMin<<endl; */
      doneRotating = true;
    }

    // rotational force along the rotational planes normal
    rotationalForce1 =
        (rotationalForce.array() * rotationalPlane.array()).sum();

    rotate(params->dimer.rotationAngle);

    if (!doneRotating) {
      // rotated dimer
      curvature = calcRotationalForceReturnCurvature(rotationalForce);

      rotationalForce2 =
          (rotationalForce.array() * rotationalPlane.array()).sum();

      rotationalForceChange =
          ((rotationalForce1 - rotationalForce2) / params->dimer.rotationAngle);

      forceDimer = (rotationalForce1 + rotationalForce2) / 2.0;

      rotationAngle = (atan(2.0 * forceDimer / rotationalForceChange) / 2.0 -
                       params->dimer.rotationAngle / 2.0);

      //            std::cout << "Rotation Angle: " <<rotationAngle <<endl;
      //            //debug

      if (rotationalForceChange < 0) {
        rotationAngle = rotationAngle + M_PI / 2.0;
      }

      rotate(rotationAngle);
      rotationalPlaneOld = rotationalPlane; // XXX: Is this copying correctly???
      rotations++;
    }
    SPDLOG_LOGGER_DEBUG(log,
                        "[DimerRot]   -----   ---------   ----------------   "
                        "---------  {:9.3e}  {:9.3e}  {:9.3e}   ---------\n",
                        curvature, torque, rotationAngle * (180.0 / M_PI));
  }

  statsTorque = torque;
  statsCurvature = curvature;
  direction.normalize();
  statsAngle = acos((direction.array() * initialDirection.array()).sum());
  statsAngle *= (180.0 / M_PI);
  statsRotations = rotations;

  eigenvalue = curvature;

  forceCallsCenter = matterCenter->getForceCalls() - forceCallsCenter;
  forceCallsDimer = matterDimer->getForceCalls() - forceCallsDimer;

  totalForceCalls += forceCallsCenter + forceCallsDimer;

  return;
}

double Dimer::getEigenvalue() { return eigenvalue; }

AtomMatrix Dimer::getEigenvector() { return direction; }

double Dimer::calcRotationalForceReturnCurvature(AtomMatrix &rotationalForce) {
  double projectedForceA, projectedForceB;
  AtomMatrix posCenter(nAtoms, 3);
  AtomMatrix posDimer(nAtoms, 3);
  AtomMatrix forceCenter(nAtoms, 3);
  AtomMatrix forceA(nAtoms, 3);
  AtomMatrix forceB(nAtoms, 3);

  posCenter = matterCenter->getPositions();

  // displace to get the dimer configuration A
  posDimer = posCenter + direction * params->main.finiteDifference;

  // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
  // http://doi.org/10.1021/ct501155k
  if (params->dimer.removeRotation) {
    matterDimer->setPositions(posDimer);
    rotationRemove(matterCenter, matterDimer);
    posDimer = matterDimer->getPositions();
    direction = posDimer - posCenter;
    direction.normalize();
    posDimer = posCenter + direction * params->main.finiteDifference;
  }

  // obtain the force for the dimer configuration
  matterDimer->setPositions(posDimer);
  forceA = matterDimer->getForces();

  // use forward difference to obtain the force for configuration B
  forceCenter = matterCenter->getForces();
  forceB = 2.0 * forceCenter - forceA;

  projectedForceA = (direction.array() * forceA.array()).sum();
  projectedForceB = (direction.array() * forceB.array()).sum();

  // remove force component parallel to dimer
  forceA = makeOrthogonal(forceA, direction);
  forceB = makeOrthogonal(forceB, direction);

  // determine difference in force orthogonal to dimer
  rotationalForce = (forceA - forceB) / (2.0 * params->main.finiteDifference);

  // curvature along the dimer
  return (projectedForceB - projectedForceA) /
         (2.0 * params->main.finiteDifference);
}

void Dimer::determineRotationalPlane(AtomMatrix rotationalForce,
                                     AtomMatrix &rotationalForceOld,
                                     AtomMatrix rotationalPlaneOld,
                                     double *lengthRotationalForceOld) {
  double a, b, gamma = 0.0;

  a = fabs((rotationalForce.array() * rotationalForceOld.array()).sum());
  b = rotationalForceOld.squaredNorm();
  if (a < 0.5 * b) {
    gamma = (rotationalForce.array() *
             (rotationalForce - rotationalForceOld).array())
                .sum() /
            b;
  } else
    gamma = 0.0;

  // new rotational plane based on the current rotational force and the previous
  // rotational plane force
  rotationalPlane = rotationalForce +
                    rotationalPlaneOld * (*(lengthRotationalForceOld)) * gamma;

  // plane normal is made orthogonal to the dimer direction and normalized
  *lengthRotationalForceOld = rotationalPlane.norm();
  rotationalPlane = makeOrthogonal(rotationalPlane, direction);
  rotationalPlane.normalize();

  rotationalForceOld = rotationalForce;

  return;
}

void Dimer::rotate(double rotationAngle) {
  double cosAngle, sinAngle;

  statsAngle += rotationAngle;

  cosAngle = cos(rotationAngle);
  sinAngle = sin(rotationAngle);

  direction = direction * cosAngle + rotationalPlane * sinAngle;
  rotationalPlane = rotationalPlane * cosAngle - direction * sinAngle;

  direction.normalize();
  rotationalPlane.normalize();

  // remove component from rotationalPlane parallel to direction
  rotationalPlane = makeOrthogonal(rotationalPlane, direction);
  rotationalPlane.normalize();

  return;
}

} // namespace eonc
