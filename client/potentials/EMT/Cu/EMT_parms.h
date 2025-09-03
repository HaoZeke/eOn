/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
// all information for the EMT potential

#ifndef EMT_PARMS_H
#define EMT_PARMS_H

// the number of different spicies that will be used in simulation
int const EMT_NumType = 1;
// int const EMT_NumType = 2;

// array where number of elements should corresponds to EMT_NumType
// the value of element [0] correspond to the number of atoms of the first type
// the value of element [1] correspond to the number of atoms of the second type
// and so on
int const EMT_NumAtomPerType[] = {1300};
// int const EMT_NumAtomPerType[] = {1470, 1};

// array where number of elements should corresponds to EMT_NumType
// the value of element [0] correspond to the atomic number of the atoms of the
// first type the value of element [1] correspond to the atomic number of the
// atoms of the second type and so on the last type is the type of atoms that
// eventually will be deposited!
int const EMT_AtomType[] = {29};
// int const EMT_AtomType[] = {12,1};

#endif
