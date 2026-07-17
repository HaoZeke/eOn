/*
** MinModeSaddleSearch — first-class single-ended saddle search.
*/
#include "Matter.h"
#include "MinModeSaddleSearch.h"
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

void bind_saddle(nb::module_ &m) {
  using eonc::Matter;
  using eonc::MinModeSaddleSearch;
  using eonc::Parameters;
  using eonc::Potential;

  nb::enum_<MinModeSaddleSearch::Status>(m, "SaddleStatus")
      .value("GOOD", MinModeSaddleSearch::STATUS_GOOD)
      .value("INIT", MinModeSaddleSearch::STATUS_INIT)
      .value("BAD_NO_CONVEX", MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
      .value("BAD_HIGH_ENERGY", MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
      .value("BAD_MAX_CONCAVE_ITERATIONS",
             MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
      .value("BAD_MAX_ITERATIONS",
             MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
      .value("BAD_NOT_CONNECTED",
             MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
      .value("BAD_PREFACTOR", MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
      .value("BAD_HIGH_BARRIER", MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
      .value("BAD_MINIMA", MinModeSaddleSearch::STATUS_BAD_MINIMA)
      .value("FAILED_PREFACTOR",
             MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
      .value("POTENTIAL_FAILED",
             MinModeSaddleSearch::STATUS_POTENTIAL_FAILED)
      .value("NONNEGATIVE_ABORT",
             MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
      .value("NONLOCAL_ABORT", MinModeSaddleSearch::STATUS_NONLOCAL_ABORT)
      .value("NEGATIVE_BARRIER",
             MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
      .value("BAD_MD_TRAJECTORY_TOO_SHORT",
             MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
      .value("BAD_NO_NEGATIVE_MODE_AT_SADDLE",
             MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
      .value("BAD_NO_BARRIER", MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
      .value("ZEROMODE_ABORT", MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
      .value("OPTIMIZER_ERROR", MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
      .value("DIMER_LOST_MODE", MinModeSaddleSearch::STATUS_DIMER_LOST_MODE)
      .value("DIMER_RESTORED_BEST",
             MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST)
      .export_values();

  m.def(
      "saddle_status_message",
      [](int status) {
        return std::string(MinModeSaddleSearch::statusMessage(status));
      },
      nb::arg("status"), "Human-readable MinModeSaddleSearch status string");

  nb::class_<MinModeSaddleSearch>(
      m, "MinModeSaddleSearch",
      "Single-ended min-mode following saddle search. Mutates the Matter "
      "geometry to the saddle. Min-mode method follows "
      "Parameters.saddle_minmode_method (dimer / lanczos / davidson).")
      .def(
          "__init__",
          [](MinModeSaddleSearch *self, std::shared_ptr<Matter> matter,
             const NpF64 &mode, double reactant_energy,
             const Parameters &params, std::shared_ptr<Potential> pot) {
            AtomMatrix mode_am = atom_matrix_from_numpy(mode);
            new (self) MinModeSaddleSearch(std::move(matter), mode_am,
                                           reactant_energy, params,
                                           std::move(pot));
          },
          nb::arg("matter"), nb::arg("mode"), nb::arg("reactant_energy"),
          nb::arg("parameters"), nb::arg("potential"), nb::keep_alive<1, 2>(),
          nb::keep_alive<1, 5>(), nb::keep_alive<1, 6>(),
          "matter: starting (displaced) geometry; mode: float64 (n,3) initial "
          "min-mode; reactant_energy: energy of the und placed reactant for "
          "barrier checks")
      .def(
          "run",
          [](MinModeSaddleSearch &self) {
            nb::gil_scoped_release release;
            return self.run();
          },
          "Run saddle search; returns SaddleStatus int. On GOOD, matter is at "
          "the saddle.")
      .def(
          "run",
          [](MinModeSaddleSearch &self, long max_iter) {
            nb::gil_scoped_release release;
            return self.run(max_iter);
          },
          nb::arg("max_iterations"),
          "Run with iteration cap override")
      .def_prop_ro(
          "status",
          [](const MinModeSaddleSearch &self) { return self.getStatus(); })
      .def_prop_ro(
          "status_message",
          [](const MinModeSaddleSearch &self) {
            return std::string(
                MinModeSaddleSearch::statusMessage(self.getStatus()));
          })
      .def_prop_ro(
          "iteration",
          [](const MinModeSaddleSearch &self) {
            return self.getIterationCount();
          })
      .def_prop_ro(
          "force_calls",
          [](const MinModeSaddleSearch &self) {
            return self.getForceCalls();
          })
      .def_prop_ro(
          "eigenvalue",
          [](MinModeSaddleSearch &self) { return self.getEigenvalue(); })
      .def_prop_ro(
          "eigenvector",
          [](MinModeSaddleSearch &self) {
            return matrix_to_numpy(self.getEigenvector());
          },
          nb::rv_policy::move)
      .def("__repr__", [](const MinModeSaddleSearch &self) {
        return "<MinModeSaddleSearch status=" +
               std::to_string(self.getStatus()) + " iter=" +
               std::to_string(self.getIterationCount()) + ">";
      });
}

} // namespace eonc::pybind
