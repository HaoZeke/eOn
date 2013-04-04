//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MinModeSaddleSearch.h"
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "Lanczos.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "EpiCenters.h"
#include "ObjectiveFunction.h"
#include "Log.h"

#include <cstdlib>

using namespace helper_functions;

class MinModeObjectiveFunction : public ObjectiveFunction
{
    public:
        MinModeObjectiveFunction(Matter *matterPassed, LowestEigenmode *minModeMethodPassed,
                                 AtomMatrix modePassed, Parameters *parametersPassed)
        {
            matter = matterPassed;
            minModeMethod = minModeMethodPassed;
            eigenvector = modePassed;
            parameters = parametersPassed;
        }
        ~MinModeObjectiveFunction(void){};

        VectorXd getGradient(bool fdstep=false)
        {
            AtomMatrix proj;
            AtomMatrix force = matter->getForces();

            if (!fdstep) {
                minModeMethod->compute(matter, eigenvector);
            }

            eigenvector = minModeMethod->getEigenvector();
            double eigenvalue = minModeMethod->getEigenvalue();

            proj = (force.cwise() * eigenvector).sum() * eigenvector.normalized();

            if (0 < eigenvalue) {
                if (parameters->saddlePerpForceRatio > 0.0) {
                    // reverse force parallel to eigenvector, and reduce perpendicular force
                    double const d = parameters->saddlePerpForceRatio;
                    force = d*force - (1.+d)*proj;

                // zero out the smallest forces to keep displacement confined
                }else if(parameters->saddleConfinePositive) {
                    int sufficientForce = 0;
                    double minForce = parameters->saddleConfinePositiveMinForce;
                    while(sufficientForce < parameters->saddleConfinePositiveMinActive){
                        sufficientForce = 0;
                        force = matter->getForces();
                        for (int i=0; i<3*matter->numberOfAtoms(); i++) {
                            if (fabs(force[i]) < minForce)
                                force[i] = 0;
                            else{
                                sufficientForce = sufficientForce + 1;
                                force[i] = -parameters->saddleConfinePositiveBoost*proj[i];
                            }
                        }
                        minForce *= parameters->saddleConfinePositiveScaleRatio;
                    }
                }else{
                    // follow eigenmode
                    force = -proj;
                }
            }else{
                // reversing force parallel to eigenmode
                force += -2.*proj;
            }

            VectorXd forceV = VectorXd::Map(force.data(), 3*matter->numberOfAtoms());
            return -forceV;
        }
        double getEnergy() { return matter->getPotentialEnergy(); }
        void setPositions(VectorXd x) { matter->setPositionsV(x); }
        VectorXd getPositions() { return matter->getPositionsV(); }
        int degreesOfFreedom() { return 3*matter->numberOfAtoms(); }
        bool isConverged() { return getConvergence() < parameters->saddleConvergedForce; }

        double getConvergence() {
            if (parameters->optConvergenceMetric == "norm") {
                return matter->getForcesFreeV().norm(); 
            } else if (parameters->optConvergenceMetric == "max_atom") {
                return matter->maxForce(); 
            } else if (parameters->optConvergenceMetric == "max_component") {
                return matter->getForces().maxCoeff(); 
            } else {
                log("[MinModeSaddleSearch] unknown opt_convergence_metric: %s\n", parameters->optConvergenceMetric.c_str());
                exit(1);
            }
        }

        VectorXd difference(VectorXd a, VectorXd b) {
            return matter->pbcV(a-b);
        }

    private:
        AtomMatrix eigenvector;
        LowestEigenmode *minModeMethod;
        Matter *matter;
        Parameters *parameters;
};

MinModeSaddleSearch::MinModeSaddleSearch(Matter *matterPassed, AtomMatrix modePassed,
                          double reactantEnergyPassed, Parameters *parametersPassed)
{
    reactantEnergy = reactantEnergyPassed;
    matter = matterPassed;
    mode = modePassed;
    parameters = parametersPassed;
    status = STATUS_GOOD;
    iteration = 0;

    if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        if (parameters->dimerImproved) {
            minModeMethod = new ImprovedDimer(matter, parameters);
        }else{
            minModeMethod = new Dimer(matter, parameters);
        }
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        minModeMethod = new Lanczos(matter, parameters);
    }
}

MinModeSaddleSearch::~MinModeSaddleSearch()
{
    delete minModeMethod;
}

int MinModeSaddleSearch::run()
{
    log("Saddle point search started from reactant with energy %f eV.\n", reactantEnergy);

    if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        log("[Dimer]  %9s   %9s   %10s   %9s   %9s   %7s   %6s   %4s\n", 
            "Step", "Step Size", "Delta E", "Force", "Curvature", 
            "Torque", "Angle", "Rots");
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        log("[Lanczos]  %9s  %9s  %10s  %9s  %9s\n", 
            "Step", "Step Size", "Delta E", "Force", "Curvature");
    }

    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        matter->matter2con(climb.str(), false);
    }
    
    AtomMatrix initialPosition = matter->getPositions();

    MinModeObjectiveFunction objf(matter, minModeMethod, mode, parameters);
    if (parameters->saddleNonnegativeDisplacementAbort) {
        objf.getGradient();
        if (minModeMethod->getEigenvalue() > 0) {
            printf("%f\n", minModeMethod->getEigenvalue());
            return STATUS_NONNEGATIVE_ABORT;
        }
    }

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    while (!objf.isConverged() || iteration == 0) {
        
        if(parameters->saddleNonlocalCountAbort != 0) {
            long nm = numAtomsMoved(initialPosition - matter->getPositions(), 
                                    parameters->saddleNonlocalDistanceAbort);
            if(nm >= parameters->saddleNonlocalCountAbort) {
                status = STATUS_NONLOCAL_ABORT;
                break;
            }
        }
        
        if (iteration >= parameters->saddleMaxIterations) {
            status = STATUS_BAD_MAX_ITERATIONS;
            break;
        }

        AtomMatrix pos = matter->getPositions();

        optimizer->step(parameters->optMaxMove);

        double de = objf.getEnergy()-reactantEnergy;
        double stepSize = helper_functions::maxAtomMotion(matter->pbc(matter->getPositions() - pos));

        if (de > parameters->saddleMaxEnergy) {
            status = STATUS_BAD_HIGH_ENERGY;
            break;
        }
        
        iteration++;

        if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER)
        {
            log("[Dimer]  %9ld   %9.7f   %10.4f   %9.5f   %9.4f   %7.3f   %6.3f   %4ld\n",
                        iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                        matter->maxForce(),
                        minModeMethod->getEigenvalue(),
                        minModeMethod->statsTorque,
                        minModeMethod->statsAngle,
                        minModeMethod->statsRotations);
        }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
            log("[Lanczos]  %9ld  % 9.6f   %10.4f  %9.5f  %9.5f\n", 
                iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                matter->maxForce(),
                minModeMethod->getEigenvalue());
        }

        if (parameters->writeMovies) {
            matter->matter2con(climb.str(), true);
        }

        if (parameters->checkpoint) {
            matter->matter2con("displacement_cp.con", false);
            FILE *fileMode = fopen("mode_cp.dat", "wb");
            helper_functions::saveMode(fileMode, matter, 
                                       minModeMethod->getEigenvector());
            fclose(fileMode);
        }
    }

    if (iteration == 0) minModeMethod->compute(matter, mode);

    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
        log("[MinModeSaddleSearch] eigenvalue not negative\n");
        status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }

    delete optimizer;

    return status;
}

double MinModeSaddleSearch::getEigenvalue()
{
    return minModeMethod->getEigenvalue();
}

AtomMatrix MinModeSaddleSearch::getEigenvector()
{
    return minModeMethod->getEigenvector();
}
