/*
** Min-mode solvers: Dimer, ImprovedDimer, Lanczos, Davidson.
** First-class surface (not Job.run workdir wrappers).
*/
#include "Davidson.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "eigen_numpy.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <memory>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

/// Flatten AtomMatrix (n,3) row-major into VectorXd of length 3n.
VectorXd flatten_atom_matrix(const AtomMatrix &am) {
  VectorXd v(am.rows() * 3);
  for (Eigen::Index i = 0; i < am.rows(); ++i)
    for (Eigen::Index j = 0; j < 3; ++j)
      v(i * 3 + j) = am(i, j);
  return v;
}

template <typename Class>
nb::class_<Class> &bind_minmode_common(nb::class_<Class> &cls) {
  using eonc::Matter;
  using eonc::Parameters;
  using eonc::Potential;

  cls.def(
         "__init__",
         [](Class *self, std::shared_ptr<Matter> matter,
            const Parameters &params, std::shared_ptr<Potential> pot) {
           new (self) Class(std::move(matter), params, std::move(pot));
         },
         nb::arg("matter"), nb::arg("parameters"), nb::arg("potential"),
         nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
         nb::keep_alive<1, 4>(),
         "Matter + Parameters + Potential (params must outlive the solver)")
      .def(
          "compute",
          [](Class &self, std::shared_ptr<Matter> matter, const NpF64 &dir) {
            AtomMatrix direction = atom_matrix_from_numpy(dir);
            nb::gil_scoped_release release;
            self.compute(std::move(matter), direction);
          },
          nb::arg("matter"), nb::arg("direction"),
          "Converge lowest eigenmode. direction: float64 (n_atoms, 3)")
      .def_prop_ro(
          "eigenvalue", [](Class &self) { return self.getEigenvalue(); },
          "Lowest curvature after compute()")
      .def_prop_ro(
          "eigenvector",
          [](Class &self) { return matrix_to_numpy(self.getEigenvector()); },
          nb::rv_policy::move, "Mode as float64 (n_atoms, 3)")
      .def_prop_ro("total_force_calls",
                   [](const Class &self) { return self.totalForceCalls; })
      .def_prop_ro("total_iterations",
                   [](const Class &self) { return self.totalIterations; })
      .def_prop_ro("stats_torque",
                   [](const Class &self) { return self.statsTorque; })
      .def_prop_ro("stats_curvature",
                   [](const Class &self) { return self.statsCurvature; })
      .def_prop_ro("stats_angle",
                   [](const Class &self) { return self.statsAngle; })
      .def_prop_ro("stats_rotations",
                   [](const Class &self) { return self.statsRotations; });
  return cls;
}

} // namespace

void bind_eigenmode(nb::module_ &m) {
  {
    auto cls = nb::class_<eonc::Dimer>(
        m, "Dimer",
        "Classic finite-difference dimer min-mode (LowestEigenmode)");
    bind_minmode_common(cls);
  }
  {
    auto cls = nb::class_<eonc::ImprovedDimer>(
        m, "ImprovedDimer",
        "Improved dimer min-mode (default eOn min-mode; CG/LBFGS rotation)");
    bind_minmode_common(cls);
    cls.def(
           "set_reference_mode",
           [](eonc::ImprovedDimer &self, const NpF64 &mode) {
             if (mode.ndim() == 2)
               self.setReferenceMode(
                   flatten_atom_matrix(atom_matrix_from_numpy(mode)));
             else
               self.setReferenceMode(vector_from_numpy(mode));
           },
           nb::arg("mode"),
           "Lock rotation reference mode (OCI / mode tracking); (n,3) or 3n")
        .def("clear_reference_mode",
             &eonc::ImprovedDimer::clearReferenceMode)
        .def_prop_ro("rotation_did_converge",
                     [](const eonc::ImprovedDimer &s) {
                       return s.rotationDidConverge;
                     })
        .def_prop_ro("found_negative_curvature",
                     [](const eonc::ImprovedDimer &s) {
                       return s.foundNegativeCurvature;
                     });
  }
  {
    auto cls = nb::class_<eonc::Lanczos>(
        m, "Lanczos",
        "Lanczos iteration for lowest Hessian curvature mode");
    bind_minmode_common(cls);
  }
  {
    auto cls = nb::class_<eonc::Davidson>(
        m, "Davidson",
        "Davidson residual-correction subspace min-mode solver");
    bind_minmode_common(cls);
  }
}

} // namespace eonc::pybind
