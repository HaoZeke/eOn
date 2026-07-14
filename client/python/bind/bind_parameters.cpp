#include "Parameters.h"
#include <magic_enum/magic_enum.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <memory>
#include <stdexcept>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_parameters(nb::module_ &m) {
  nb::class_<eonc::Parameters>(m, "Parameters",
                               "Client runtime parameters (config.ini / JSON)")
      .def(nb::init<>())
      .def(
          "load",
          [](eonc::Parameters &self, const std::string &path) {
            int err = self.load(path);
            if (err) {
              throw std::runtime_error("Parameters.load failed for " + path);
            }
          },
          nb::arg("path"), "Load from config.ini path")
      .def(
          "load_json",
          [](eonc::Parameters &self, const std::string &json_str) {
            int err = self.load_json(json_str);
            if (err) {
              throw std::runtime_error("Parameters.load_json failed");
            }
          },
          nb::arg("json_str"))
      .def("to_json", &eonc::Parameters::to_json)
      .def_prop_rw(
          "job",
          [](const eonc::Parameters &self) { return self.main_options.job; },
          [](eonc::Parameters &self, eonc::JobType j) {
            self.main_options.job = j;
          })
      .def_prop_rw(
          "potential",
          [](const eonc::Parameters &self) {
            return self.potential_options.potential;
          },
          [](eonc::Parameters &self, eonc::PotType p) {
            self.potential_options.potential = p;
          })
      .def_prop_rw(
          "temperature",
          [](const eonc::Parameters &self) {
            return self.main_options.temperature;
          },
          [](eonc::Parameters &self, double t) {
            self.main_options.temperature = t;
          })
      .def_prop_rw(
          "random_seed",
          [](const eonc::Parameters &self) {
            return self.main_options.randomSeed;
          },
          [](eonc::Parameters &self, long s) {
            self.main_options.randomSeed = s;
          })
      .def_prop_rw(
          "quiet",
          [](const eonc::Parameters &self) { return self.main_options.quiet; },
          [](eonc::Parameters &self, bool q) { self.main_options.quiet = q; })
      .def_prop_rw(
          "write_log",
          [](const eonc::Parameters &self) {
            return self.main_options.writeLog;
          },
          [](eonc::Parameters &self, bool w) {
            self.main_options.writeLog = w;
          })
      .def_prop_rw(
          "remove_net_force",
          [](const eonc::Parameters &self) {
            return self.main_options.removeNetForce;
          },
          [](eonc::Parameters &self, bool v) {
            self.main_options.removeNetForce = v;
          })
      .def_prop_rw(
          "write_con_forces",
          [](const eonc::Parameters &self) {
            return self.main_options.writeConForces;
          },
          [](eonc::Parameters &self, bool v) {
            self.main_options.writeConForces = v;
          })
      .def_prop_rw(
          "potentials_path",
          [](const eonc::Parameters &self) {
            return self.potential_options.potentialsPath;
          },
          [](eonc::Parameters &self, const std::string &p) {
            self.potential_options.potentialsPath = p;
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
