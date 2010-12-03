//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef SADDLE_POINT_H
#define SADDLE_POINT_H

#include "Matter.h"
#include "LowestEigenmodeInterface.h" 

#include <string>

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

using namespace std;

class Matter;
class Parameters;
class LowestEigenmodeInterface;

/** Rely on the dimer method for saddle point determination. The object rely on an object being able to determine the lowest eigenmode and a function to determine where a displacement prior the saddle point search should be centered.*/
class SaddlePoint {
public:

    // Return codes passed from server to client to indicate calculation status
    enum{
        statusGood=0,
        statusInit,
        statusBadNoConvex,
        statusBadHighEnergy,
        statusBadMaxConcaveIterations,
        statusBadMaxIterations,
        statusBadNotConnected,
        statusBadPrefactor,
        statusBadHighBarrier,
        statusBadMinima
    };

    // Constants used to displace atoms before a saddle search
    enum{
        dispNone=0,
        dispNotFccOrHcp,
        dispMinCoordinated,
        dispLastAtom
    };

    // Methods for finding the minimum mode
    enum{
        MINMODE_DIMER,
        MINMODE_LANCZOS
    };
 
    SaddlePoint(); // The object shall be initialized later with SaddlePoint::initialize
 
    /** Constructor where object is initialized.
    @param[in]  initial      Pointer to where the initial state (minimum) shall be stored.
    @param[in]  saddle       Pointer to where the conformation where to start the saddle point search. It will also be used to return the saddle point.
    @param[in]  *parameters  Pointer to the Parameter object containing the runtime parameters */
    SaddlePoint(Matter *initial, Matter *saddle, Parameters *parameters);

    ~SaddlePoint(); // destructor

    /** Initialized the object.
    @param[in]  *initial     Initial state (minimum)
    @param[in]  *saddle      Where to start the saddle point search (may be equal to initial)
    @param[in]  *parameters  Pointer to the Parameter object containing the runtime parameters
    */
    void initialize(Matter *initial, Matter *saddle, Parameters *parameters);

    /** Try to determine a nearby saddle point
    @param[out]  *min1     Pointer to Matter object containing one of the minima connected to the saddle point
    @param[out]  *min2     Pointer to Matter object containing the other minima connected to the saddle point
    The values returned is true if the calculation converged.*/
    long locate(Matter *min1, Matter *min2);
    LowestEigenmodeInterface const * getLowestEigenmode() const;
    long getnFreeCoord() const;
    Matrix<double, Eigen::Dynamic, 3> getEigenMode();

    Matrix<double, Eigen::Dynamic, 3> mode;
    void loadMode(string filename);
    void loadMode(FILE * modeFile);
    void saveMode(FILE * modeFile);

    long forceCallsSaddlePointConcave;
    long forceCallsSaddlePointConvex;
    long forceCallsMinimization;

private:
    Matter *initial;
    Matter *saddle; // pointer to atom object outside the scope of the class
    Parameters *parameters; // pointer to a structure outside the scope of the class containing runtime parameters
    LowestEigenmodeInterface *lowestEigenmode; // pointer to the method used to determine the lowest eigenmode
 
    double eigenValue; // containing an estimate for the lowest eigenvalue
    Matrix<double, Eigen::Dynamic, 3> eigenMode; // double array for the lowest eigenmode
    Matrix<double, Eigen::Dynamic, 3> initialDisplacement; //RT: used to keep track of the initial displacement
    long nFreeCoord; // number of free coordinates
    long status; // keep track of where problems occured

    void clean(); // clean up dynamical allocated memory

    void displaceStateAndSetMode(Matter *matter); // displacing the atomic positions in @matter according to values in Parameters, being centered on the atom determined by one of the EpiCenter functions and set the initial mode accordingly  
    Matrix<double, Eigen::Dynamic, 3> correctingForces(Matrix<double, Eigen::Dynamic, 3> force); // projected min-mode force

    /** Determine the two minima connected to the saddle point, by displacing the positions in the saddle point by either adding or subtracting a part of the lowest eigenmode
    @param[out]  *min1   Matter object containing one of the minima connected to the saddle point
    @param[out]  *min2   Matter object containing other of the minima connected to the saddle point */
    void relaxFromSaddle(Matter *min1, Matter *min2);

    void jumpToConvexRegion();
    void displaceInConcaveRegion();

    void searchForSaddlePoint(double initialEnergy);
    void addForceCallsSaddlePoint(long fcalls, double eigenvalue);

};

#endif
