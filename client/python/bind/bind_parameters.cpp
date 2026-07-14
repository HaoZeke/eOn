#include "Parameters.h"

#include <magic_enum/magic_enum.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_parameters(nb::module_ &m) {
  nb::class_<eonc::Parameters>(m, "Parameters",
                               "Client runtime parameters (config.ini / JSON)")
      .def(nb::init<>())
      .def(
          "load",
          [](eonc::Parameters &self, const std::string &path) {
            if (self.load(path))
              throw std::runtime_error("Parameters.load failed for " + path);
          },
          nb::arg("path"))
      .def(
          "load_json",
          [](eonc::Parameters &self, const std::string &json_str) {
            if (self.load_json(json_str))
              throw std::runtime_error("Parameters.load_json failed");
          },
          nb::arg("json_str"))
      .def("to_json", &eonc::Parameters::to_json)
      // --- Main ---
      .def_prop_rw(
          "job",
          [](const eonc::Parameters &s) { return s.main_options.job; },
          [](eonc::Parameters &s, eonc::JobType j) { s.main_options.job = j; })
      .def_prop_rw(
          "potential",
          [](const eonc::Parameters &s) {
            return s.potential_options.potential;
          },
          [](eonc::Parameters &s, eonc::PotType p) {
            s.potential_options.potential = p;
          })
      .def_prop_rw(
          "temperature",
          [](const eonc::Parameters &s) { return s.main_options.temperature; },
          [](eonc::Parameters &s, double t) { s.main_options.temperature = t; })
      .def_prop_rw(
          "random_seed",
          [](const eonc::Parameters &s) { return s.main_options.randomSeed; },
          [](eonc::Parameters &s, long v) { s.main_options.randomSeed = v; })
      .def_prop_rw(
          "quiet",
          [](const eonc::Parameters &s) { return s.main_options.quiet; },
          [](eonc::Parameters &s, bool v) { s.main_options.quiet = v; })
      .def_prop_rw(
          "write_log",
          [](const eonc::Parameters &s) { return s.main_options.writeLog; },
          [](eonc::Parameters &s, bool v) { s.main_options.writeLog = v; })
      .def_prop_rw(
          "checkpoint",
          [](const eonc::Parameters &s) { return s.main_options.checkpoint; },
          [](eonc::Parameters &s, bool v) { s.main_options.checkpoint = v; })
      .def_prop_rw(
          "remove_net_force",
          [](const eonc::Parameters &s) {
            return s.main_options.removeNetForce;
          },
          [](eonc::Parameters &s, bool v) {
            s.main_options.removeNetForce = v;
          })
      .def_prop_rw(
          "write_con_forces",
          [](const eonc::Parameters &s) {
            return s.main_options.writeConForces;
          },
          [](eonc::Parameters &s, bool v) {
            s.main_options.writeConForces = v;
          })
      .def_prop_rw(
          "con_filename",
          [](const eonc::Parameters &s) { return s.main_options.conFilename; },
          [](eonc::Parameters &s, const std::string &v) {
            s.main_options.conFilename = v;
          })
      .def_prop_rw(
          "ini_filename",
          [](const eonc::Parameters &s) { return s.main_options.iniFilename; },
          [](eonc::Parameters &s, const std::string &v) {
            s.main_options.iniFilename = v;
          })
      // --- Potential ---
      .def_prop_rw(
          "potentials_path",
          [](const eonc::Parameters &s) {
            return s.potential_options.potentialsPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.potential_options.potentialsPath = v;
          })
      .def_prop_rw(
          "ext_pot_path",
          [](const eonc::Parameters &s) {
            return s.potential_options.extPotPath;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.potential_options.extPotPath = v;
          })
      // --- Optimizer ---
      .def_prop_rw(
          "opt_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.optimizer_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            s.optimizer_options.max_iterations = static_cast<size_t>(v);
          })
      .def_prop_rw(
          "opt_converged_force",
          [](const eonc::Parameters &s) {
            return s.optimizer_options.converged_force;
          },
          [](eonc::Parameters &s, double v) {
            s.optimizer_options.converged_force = v;
          })
      .def_prop_rw(
          "opt_max_move",
          [](const eonc::Parameters &s) { return s.optimizer_options.max_move; },
          [](eonc::Parameters &s, double v) {
            s.optimizer_options.max_move = v;
          })
      // --- Structure comparison ---
      .def_prop_rw(
          "comp_eps_r",
          [](const eonc::Parameters &s) {
            return s.structure_comparison_options.distance_difference;
          },
          [](eonc::Parameters &s, double v) {
            s.structure_comparison_options.distance_difference = v;
          })
      .def("__repr__", [](const eonc::Parameters &self) {
        return "<Parameters job=" +
               std::string(magic_enum::enum_name(self.main_options.job)) +
               " pot=" +
               std::string(
                   magic_enum::enum_name(self.potential_options.potential)) +
               ">";
      });
}

} // namespace eonc::pybind
