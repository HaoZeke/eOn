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
#define NPARAMS 3

const EAM::element_parameters EAM::el_params[NPARAMS] = {
    /*
     * These are the definitions for element parameters and look like:
     * { Z, Dm, alphaM, Rm, beta1, beta2, r_cut, func_coeff[] }
     */

    /* Aluminum -- 13 */
    {13,
     3.63262,
     1.49303,
     2.12396,
     3.45762,
     3.54208,
     5.555,
     {-0.0847305, 17.0887, 60.7909, -75.8412, 39.9948, -6.21563, -3.07993,
      1.51365, -0.196034}},
    /* Palladium -- 46 */
    {46,
     1.6054,
     1.55083,
     2.35604,
     3.35148,
     3.34715,
     5.4120,
     {0.0992627, 4.8356, -0.368647, 17.2739, -18.450, 9.5018, -2.69876,
      0.405341, -0.0251768}},
    /* Gold -- 79 */
    {79,
     0.669134,
     1.89707,
     2.57021,
     3.69785,
     3.69674,
     5.5155,
     {-0.198791, -27.7663, 224.068, -594.101, 1316.31, -2070.62, 1973.14,
      -1012.02, 214.373}}};
