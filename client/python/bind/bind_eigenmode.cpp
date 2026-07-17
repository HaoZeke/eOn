/*
** Min-mode: chemist-facing Dimer(method=..., accelerant=...) plus low-level classes.
** Default method is "improved". accelerant="gp" routes to AtomicGPDimer (WITH_GPRD).
*/
#include "Davidson.h"
#include "Dimer.h"
#include "EigenmodeStrategy.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "eigen_numpy.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

VectorXd flatten_atom_matrix(const AtomMatrix &am) {
  VectorXd v(am.rows() * 3);
  for (Eigen::Index i = 0; i < am.rows(); ++i)
    for (Eigen::Index j = 0; j < 3; ++j)
      v(i * 3 + j) = am(i, j);
  return v;
}

std::string lower_ascii(std::string s) {
  for (char &c : s)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

/// Route Parameters for buildEigenmodeStrategy from chemist kwargs.
eonc::Parameters route_minmode_params(eonc::Parameters params,
                                      const std::string &method_in,
                                      const std::string &accelerant_in) {
  const std::string method = lower_ascii(method_in);
  const std::string acc = lower_ascii(accelerant_in);

  if (!acc.empty() && acc != "none" && acc != "gp") {
    throw std::runtime_error(
        "Dimer: accelerant must be None/\"\" or \"gp\", got \"" + accelerant_in +
        "\"");
  }
  const bool want_gp = (acc == "gp");

  if (want_gp) {
#ifndef WITH_GPRD
    throw std::runtime_error(
        "Dimer(accelerant=\"gp\") requires build with -Dwith_gprd=true "
        "(WITH_GPRD); this extension has built_with_gprd()==False");
#else
    // GP accelerant is always the improved-dimer + GP stack (AtomicGPDimer).
    if (method != "improved" && method != "dimer" && !method.empty()) {
      throw std::runtime_error(
          "Dimer(accelerant=\"gp\") only applies to method=\"improved\" "
          "(got \"" +
          method_in +
          "\"). Use DimerSpec / method=\"improved\" — not classic/lanczos/"
          "davidson.");
    }
    params.saddle_search_options.minmode_method =
        eonc::LowestEigenmode::MINMODE_GPRDIMER;
    params.dimer_options.improved = true;
    return params;
#endif
  }

  if (method == "classic") {
    params.saddle_search_options.minmode_method =
        eonc::LowestEigenmode::MINMODE_DIMER;
    params.dimer_options.improved = false;
  } else if (method == "improved" || method == "dimer" || method.empty()) {
    // Default and "dimer" synonym: improved dimer
    params.saddle_search_options.minmode_method =
        eonc::LowestEigenmode::MINMODE_DIMER;
    params.dimer_options.improved = true;
  } else if (method == "lanczos") {
    params.saddle_search_options.minmode_method =
        eonc::LowestEigenmode::MINMODE_LANCZOS;
  } else if (method == "davidson") {
    params.saddle_search_options.minmode_method =
        eonc::LowestEigenmode::MINMODE_DAVIDSON;
  } else {
    throw std::runtime_error(
        "Dimer: method must be improved|classic|lanczos|davidson, got \"" +
        method_in + "\"");
  }
  return params;
}

/// Chemist-facing dimer: method defaults to improved; accelerant="gp" optional.
struct PyDimer {
  eonc::Parameters params;
  std::shared_ptr<eonc::Potential> pot;
  std::shared_ptr<eonc::Matter> seed;
  std::shared_ptr<eonc::EigenmodeStrategy> strategy;
  std::string method;
  std::string accelerant;

  PyDimer(std::shared_ptr<eonc::Matter> matter, const eonc::Parameters &p_in,
          std::shared_ptr<eonc::Potential> pot_in, std::string method_in,
          std::string accelerant_in)
      : params(route_minmode_params(p_in, method_in, accelerant_in)),
        pot(std::move(pot_in)),
        seed(std::move(matter)),
        method(std::move(method_in)),
        accelerant(std::move(accelerant_in)) {
    if (!seed)
      throw std::runtime_error("Dimer: matter is required");
    if (!pot)
      pot = seed->getPotential();
    if (!pot)
      throw std::runtime_error("Dimer: potential is required");
    strategy = eonc::buildEigenmodeStrategy(seed, params, pot);
  }

  void compute(std::shared_ptr<eonc::Matter> matter, const NpF64 &dir) {
    AtomMatrix direction = atom_matrix_from_numpy(dir);
    auto m = matter ? matter : seed;
    nb::gil_scoped_release release;
    eonc::eigenmodeCompute(*strategy, std::move(m), direction);
  }

  double eigenvalue() const {
    return eonc::eigenmodeGetEigenvalue(*strategy);
  }

  auto eigenvector() const {
    return matrix_to_numpy(eonc::eigenmodeGetEigenvector(*strategy));
  }

  long total_force_calls() const {
    return std::visit([](const auto &impl) { return impl.totalForceCalls; },
                      *strategy);
  }

  long total_iterations() const {
    return std::visit([](const auto &impl) { return impl.totalIterations; },
                      *strategy);
  }

  double stats_torque() const {
    return std::visit([](const auto &impl) { return impl.statsTorque; },
                      *strategy);
  }

  double stats_curvature() const {
    return std::visit([](const auto &impl) { return impl.statsCurvature; },
                      *strategy);
  }

  double stats_angle() const {
    return std::visit([](const auto &impl) { return impl.statsAngle; },
                      *strategy);
  }

  long stats_rotations() const {
    return std::visit([](const auto &impl) { return impl.statsRotations; },
                      *strategy);
  }
};

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
         "Low-level solver: Matter + Parameters + Potential")
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
  m.def(
      "built_with_gprd",
      []() {
#ifdef WITH_GPRD
        return true;
#else
        return false;
#endif
      },
      "True if compiled with -Dwith_gprd=true (AtomicGPDimer / GP accelerant)");

  // --- Chemist entry: Dimer(method="improved", accelerant=None|"gp") ---
  nb::class_<PyDimer>(
      m, "Dimer",
      "Min-mode entry. Default method=\"improved\". "
      "accelerant=\"gp\" selects GP-dimer (WITH_GPRD) without a separate type.")
      .def(
          "__init__",
          [](PyDimer *self, std::shared_ptr<eonc::Matter> matter,
             const eonc::Parameters &params,
             std::shared_ptr<eonc::Potential> pot, const std::string &method,
             const std::string &accelerant) {
            new (self) PyDimer(std::move(matter), params, std::move(pot),
                               method, accelerant);
          },
          nb::arg("matter"), nb::arg("parameters"), nb::arg("potential"),
          nb::arg("method") = "improved", nb::arg("accelerant") = "",
          nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
          nb::keep_alive<1, 4>(),
          "method: improved (default) | classic | lanczos | davidson. "
          "accelerant: \"\" or \"gp\".")
      .def("compute", &PyDimer::compute, nb::arg("matter"),
           nb::arg("direction"),
           "Converge lowest eigenmode. direction: float64 (n_atoms, 3)")
      .def_prop_ro("eigenvalue", &PyDimer::eigenvalue)
      .def_prop_ro("eigenvector", &PyDimer::eigenvector, nb::rv_policy::move)
      .def_prop_ro("method",
                   [](const PyDimer &d) { return d.method; })
      .def_prop_ro("accelerant",
                   [](const PyDimer &d) { return d.accelerant; })
      .def_prop_ro("total_force_calls", &PyDimer::total_force_calls)
      .def_prop_ro("total_iterations", &PyDimer::total_iterations)
      .def_prop_ro("stats_torque", &PyDimer::stats_torque)
      .def_prop_ro("stats_curvature", &PyDimer::stats_curvature)
      .def_prop_ro("stats_angle", &PyDimer::stats_angle)
      .def_prop_ro("stats_rotations", &PyDimer::stats_rotations);

  // Low-level engines (power users)
  {
    auto cls = nb::class_<eonc::Dimer>(
        m, "ClassicDimer",
        "Classic finite-difference dimer (low-level). Prefer Dimer(method="
        "\"classic\").");
    bind_minmode_common(cls);
  }
  {
    auto cls = nb::class_<eonc::ImprovedDimer>(
        m, "ImprovedDimer",
        "Improved dimer engine (low-level). Prefer Dimer() / method="
        "\"improved\".");
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
        .def("clear_reference_mode", &eonc::ImprovedDimer::clearReferenceMode)
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
        "Lanczos min-mode (low-level). Prefer Dimer(method=\"lanczos\").");
    bind_minmode_common(cls);
  }
  {
    auto cls = nb::class_<eonc::Davidson>(
        m, "Davidson",
        "Davidson min-mode (low-level). Prefer Dimer(method=\"davidson\").");
    bind_minmode_common(cls);
  }
#ifdef WITH_GPRD
  {
    auto cls = nb::class_<eonc::AtomicGPDimer>(
        m, "AtomicGPDimer",
        "GP-dimer engine (low-level). Prefer Dimer(accelerant=\"gp\").");
    bind_minmode_common(cls);
  }
#endif
}

} // namespace eonc::pybind
