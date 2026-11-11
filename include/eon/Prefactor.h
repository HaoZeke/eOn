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
#include "Matter.h"

#include "Parameters.h"

#include <string_view>

namespace eonc {

namespace Prefactor {
inline constexpr int OK = 0;
// The potential or the Hessian did not return a usable spectrum.
inline constexpr int POTENTIAL_FAILED = -1;
// The spectrum is usable and is not a first-order saddle.
inline constexpr int SADDLE_ORDER = -2;
inline constexpr std::string_view RATE_HTST{"htst"};
inline constexpr std::string_view RATE_QQHTST{"qqhtst"};
inline constexpr std::string_view FILTER_CUTOFF{"cutoff"};
inline constexpr std::string_view FILTER_FRACTION{"fraction"};

VectorXi rateAtoms(const Parameters &parameters, Matter *min1, Matter *saddle,
                   Matter *min2);
int getPrefactors(const Parameters &parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2);
int getPrefactors(const Parameters &parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2,
                  const VectorXi &atoms);
VectorXi movedAtoms(const Parameters &parameters, Matter *min1, Matter *saddle,
                    Matter *min2);
VectorXi movedAtomsPct(const Parameters &parameters, Matter *min1,
                       Matter *saddle, Matter *min2);
VectorXi allFreeAtoms(Matter *matter);
void logFreqs(const VectorXd &freqs, const char *name);
} // namespace Prefactor

} // namespace eonc
