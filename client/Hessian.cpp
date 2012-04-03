//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <math.h>
#include "Hessian.h"
#include "Log.h"

Hessian::Hessian(Parameters* params, Matter *matter)
{
    parameters = params;
    hessian.resize(0,0);
    freqs.resize(0);
}

Hessian::~Hessian()
{
}

Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian::getHessian(Matter *matterIn, VectorXi atomsIn)
{
    if((matter != matterIn) || (atoms != atomsIn) || (hessian.rows() == 0))
    {
        hessian.resize(0,0);
        matter = matterIn;
        atoms = atomsIn;

        if (!calculate())
        {
            hessian.resize(0,0);
        }
    }
    return hessian;
}

VectorXd Hessian::getFreqs(Matter *matterIn, VectorXi atomsIn)
{
    if((matter != matterIn) || (atoms != atomsIn) || (hessian.rows() == 0))
    {
        hessian.resize(0,0);
        matter = matterIn;
        atoms = atomsIn;

        if (!calculate())
        {
            freqs.resize(0);
        }
    }
    return freqs;
}

bool Hessian::calculate(void)
{
    int nAtoms = matter->numberOfAtoms();

    // Determine the Hessian size
    int size = 0;
    size = atoms.rows()*3;
    cout<<"Hessian size: "<<size<<endl;
    assert(size > 0);

    // Build the hessian 
    Matter matterTemp(parameters);
    matterTemp = *matter;

    double dr = parameters->finiteDifference;

    AtomMatrix pos = matter->getPositions();
    AtomMatrix posDisplace(nAtoms, 3);
    AtomMatrix posTemp(nAtoms, 3);
    AtomMatrix force1(nAtoms, 3);
    AtomMatrix force2(nAtoms, 3);

//    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian(size, size);
    hessian.resize(size, size);

    int i, j; 
    force1 = matterTemp.getForces();
    for(i = 0; i<size; i++)
    {
        posDisplace.setZero(); 

        // Displacing one coordinate
        posDisplace(atoms(i/3), i%3) = dr;

        posTemp = pos + posDisplace; 
        matterTemp.setPositions(posTemp);
        force2 = matterTemp.getForces();

        // To use central difference estimate of the hessian uncomment following (and divide by 2*dr) in the following.
        // This does use an additional 'size' forcecalls and will generally not lead to very different results.
        // In most cases, the additional accuracy is not worth the computation time.

        /*
        posTemp = pos - posDisplace; 
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        */

        for(j=0; j<size; j++)
        {
            hessian(i,j) = -(force2(atoms(j/3), j%3)-force1(atoms(j/3), j%3))/dr;
            // get the effective mass of the moving atoms 
            double effMass = sqrt(matter->getMass(atoms(j/3))*matter->getMass(atoms(i/3)));
            hessian(i,j) /= effMass;
        }
    }

    // Force hessian to be symmetric
    
    // hessian = (hessian + hessian.transpose())/2;  
    // cannot be used, messes up the lower trianguler 
    // transpose does not seem to be a hardcopy, rather just an index manipulation 

    for(i=0; i<size; i++) {
        for(j=0; j<i; j++) {
            hessian(i,j) = (hessian(i,j) + hessian(j,i)) / 2;
            hessian(j,i) = hessian(i,j);    
        }
    }

    if(!parameters->quiet)
    {
        cout <<"writing hessian"<<endl;
        ofstream hessfile;
        hessfile.open("hessian.dat");
        hessfile <<hessian;
    }

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(hessian);
    freqs = es.eigenvalues();

    return true;
}

// If we are checking for rotation, then the system has no frozen atoms and
// can rotate and translate. This gives effectively zero eigenvalues. We
// need to remove them from the prefactor calculation. 
// the condition requires that every atom moves. Otherwise, we don't 
// get the 6 rotational and translational modes.
// XXX: what happens if the entire particle rotates about one atom or a line of atoms?

bool Hessian::removeZeroFreqs(VectorXd freqs)
{
    int size = freqs.size();
    if(size != 3*matter->numberOfAtoms()) {
        return false;
    }
    VectorXd newfreqs(size);
    int nremoved = 0;
    for(int i=0; i<size; i++)
    {
        if(abs(freqs(i)) > parameters->hessianZeroFreqValue)
        {
            newfreqs(i-nremoved) = freqs(i);
        }
        else
        {
            nremoved++;
        }
    }
    if(nremoved != 6)
    {
        log("Error: Found %i trivial eigenmodes instead of 6\n", nremoved);
    }
    freqs = newfreqs;
    return true;
}

