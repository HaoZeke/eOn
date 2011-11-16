//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <time.h>

#include "Potential.h"
#include "Parameters.h"
#include "Log.h"
#include "HelperFunctions.h"
#include "potentials/NewPotential/NewPotential.h"

#include "potentials/IMD/IMD.h"
#include "potentials/EDIP/EDIP.h"
#include "potentials/EMT/EffectiveMediumTheory.h"
#include "potentials/Morse/Morse.h"
#include "potentials/LennardJones/LJ.h"
#include "potentials/SW/SW.h"
#include "potentials/Tersoff/Tersoff.h"
#include "potentials/Aluminum/Aluminum.h"
#include "potentials/EAM/EAM.h"
#include "potentials/Lenosky/Lenosky.h"
#include "potentials/QSC/QSC.h"
//#include "potentials/water/water_for_eon.hpp"
#include "potentials/Water/Water.hpp"
#include "potentials/Water_Pt/Tip4p_Pt.hpp"
#include "potentials/Water_H/Tip4P_H.h"

#ifndef WIN32
    #include "potentials/VASP/VASP.h"
#endif
#include "potentials/bopfox/bopfox.h"
#ifdef BOPFOX
    #include "potentials/bop/bop.h"
#endif
#ifdef LAMMPS_POT
    #include "potentials/LAMMPS/LAMMPS_EON.h"
#endif
#ifdef EONMPI
    #include "potentials/MPIPot/MPIPot.h"
#endif

#include <cstdlib>

const char Potential::POT_LJ[] =          "lj";
const char Potential::POT_IMD[] =         "imd";
const char Potential::POT_EAM_AL[] =      "eam_al";
const char Potential::POT_MORSE_PT[] =    "morse_pt";
const char Potential::POT_EMT[] =         "emt";
const char Potential::POT_QSC[] =         "qsc";
const char Potential::POT_TIP4P[] =       "tip4p";
const char Potential::POT_TIP4P_PT[] =    "tip4p_pt";
const char Potential::POT_TIP4P_H[] =     "tip4p_h";
const char Potential::POT_SPCE[] =        "spce";
const char Potential::POT_LENOSKY_SI[] =  "lenosky_si";
const char Potential::POT_SW_SI[] =       "sw_si";
const char Potential::POT_TERSOFF_SI[] =  "tersoff_si";
const char Potential::POT_EDIP[] =        "edip";
const char Potential::POT_VASP[] =        "vasp";
const char Potential::POT_BOPFOX[] =      "bopfox";
const char Potential::POT_BOP[] =         "bop";
const char Potential::POT_LAMMPS[] =      "lammps";
const char Potential::POT_MPI[] =         "mpi";

Potential* Potential::pot = NULL;

Potential *Potential::getPotential(Parameters *parameters)
{
    if(pot) {
        return pot;
    }
    if(parameters->potential == POT_LJ)
        pot = new LJ();
    else if(parameters->potential == POT_IMD)
        pot = new IMD();
    else if(parameters->potential == POT_MORSE_PT)
        pot = new Morse();
    else if(parameters->potential == POT_EMT)
        pot = new EffectiveMediumTheory();
    else if(parameters->potential == POT_QSC)
        pot = new QSC();
    else if(parameters->potential == POT_TIP4P)
        pot = new Tip4p();
    else if(parameters->potential == POT_TIP4P_PT)
        pot = new Tip4p_Pt();
    else if(parameters->potential == POT_SPCE)
        pot = new SpceCcl();
    #ifndef NO_FORTRAN
    else if(parameters->potential == POT_EAM_AL)
        pot = new Aluminum();
    else if(parameters->potential == POT_LENOSKY_SI)
        pot = new Lenosky();
    else if(parameters->potential == POT_SW_SI)
        pot = new SW();
    else if(parameters->potential == POT_TERSOFF_SI)
        pot = new Tersoff();
    else if(parameters->potential == POT_EDIP)
        pot = new EDIP();
    else if(parameters->potential == POT_TIP4P_H)
        pot = new Tip4p_H();
    
    #ifndef WIN32
    else if(parameters->potential == POT_VASP)
        pot = new VASP();
    #endif
    else if(parameters->potential == POT_BOPFOX)
        pot = new bopfox();
    #endif
    #ifdef BOPFOX
    else if(parameters->potential == POT_BOP)
        pot = new bop();
    #endif
    #ifdef LAMMPS_POT
    else if(parameters->potential == POT_LAMMPS)
        pot = new lammps_eon(parameters);
    #endif
    #ifdef EONMPI
    else if(parameters->potential == POT_MPI)
        pot = new MPIPot(parameters);
    #endif
    else {
        printf("Unknown Potential: %s\n", parameters->potential.c_str());
        std::exit(1);
    }
    pot->initialize();
    return pot;
};

int Potential::fcalls = 0;
int Potential::fcallsTotal = 0;

AtomMatrix Potential::force(long nAtoms, AtomMatrix positions,
                            VectorXi atomicNrs, double *energy, Matrix3d box) 
{
    AtomMatrix forces(nAtoms,3);

    double start, userStart, sysStart;
    helper_functions::getTime(&start, &userStart, &sysStart);

    force(nAtoms, positions.data(), atomicNrs.data(), forces.data(), energy, 
          box.data());

    double finish, userFinish, sysFinish;
    helper_functions::getTime(&finish, &userFinish, &sysFinish);

    log_file("[Potential] fcall#: %4d  real: %8.2g  user: %8.2g  sys: %8.2g seconds\n", 
             fcalls, finish - start, userFinish - userStart, sysFinish - sysStart);

    fcalls+=1;
    fcallsTotal += 1;

    return forces;
};
