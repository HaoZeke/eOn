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

#include "CatLearnPot.h"

CatLearnPot::CatLearnPot(std::shared_ptr<Parameters> a_params)
    : SurrogatePotential(PotType::CatLearn, a_params) {
  py::module_ sys = py::module_::import("sys");
  py::exec(fmt::format("sys.path.insert(0, {})", a_params->catl.path));

  py::module_ gp_module = py::module_::import(
      "catlearn.regression.gaussianprocess.calculator.mlmodel");

  // Import the required modules
  // GP Model
  this->m_gpmod =
      gp_module.attr("get_default_model")("model"_a = a_params->catl.model);
};

void CatLearnPot::train_optimize(MatrixType features, MatrixType targets) {
  m_gpmod.attr("optimize")(features, targets, py::arg("retrain") = true);
  return;
}

void CatLearnPot::force(long nAtoms, const double *positions,
                        const int *atomicNrs, double *forces, double *energy,
                        double *variance, const double *box) {
  // TODO(rg): Test after matrixtype map
  MatrixType features = MatrixType::Map(positions, 1, nAtoms * 3);
  py::tuple ef_and_unc = (this->m_gpmod.attr("predict")(
      features, "get_variance"_a = true, "get_derivatives"_a = true));
  auto ef_dat = ef_and_unc[0].cast<MatrixType>();
  auto vari = ef_and_unc[1].cast<MatrixType>();
  auto gradients = ef_dat.block(0, 1, 1, nAtoms * 3);
  for (int idx = 0; idx < nAtoms; idx++) {
    forces[3 * idx] = gradients(0, 3 * idx) * -1;
    forces[3 * idx + 1] = gradients(0, 3 * idx + 1) * -1;
    forces[3 * idx + 2] = gradients(0, 3 * idx + 2) * -1;
  }
  *variance = vari(0, 0); // energy variance only
  *energy = ef_dat(0, 0);
  return;
}
