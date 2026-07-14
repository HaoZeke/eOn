#include "Parameters.h"
#include "Potential.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <magic_enum/magic_enum.hpp>
#include <memory>
#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_potential(nb::module_ &m) {
  nb::class_<eonc::Potential>(m, "Potential",
                              "Abstract potential (force + energy)")
      .def_prop_ro("type", &eonc::Potential::getType)
      .def_prop_ro(
          "force_call_counter",
          [](const eonc::Potential &self) {
            return self.forceCallCounter.load();
          })
      .def("is_surrogate", &eonc::Potential::isSurrogate)
      .def("requires_isolated_molecule_layout",
           &eonc::Potential::requiresIsolatedMoleculeLayout)
      .def("is_thread_safe", &eonc::Potential::isThreadSafe)
      .def("__repr__", [](const eonc::Potential &self) {
        return "<Potential " +
               std::string(magic_enum::enum_name(self.getType())) + ">";
      });

  m.def(
      "make_potential",
      [](eonc::PotType ptype, eonc::Parameters &params) {
        auto pot = eonc::helpers::makePotential(ptype, params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null");
        }
        return pot;
      },
      nb::arg("pot_type"), nb::arg("parameters"),
      "Construct a Potential from PotType + Parameters");

  m.def(
      "make_potential",
      [](const std::string &name, eonc::Parameters &params) {
        auto v = magic_enum::enum_cast<eonc::PotType>(
            name, magic_enum::case_insensitive);
        if (!v) {
          throw std::invalid_argument("unknown PotType: " + name);
        }
        params.potential_options.potential = *v;
        auto pot = eonc::helpers::makePotential(*v, params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null for " + name);
        }
        return pot;
      },
      nb::arg("name"), nb::arg("parameters"),
      "Construct a Potential from config-style name + Parameters");

  m.def(
      "make_potential",
      [](eonc::Parameters &params) {
        auto pot = eonc::helpers::makePotential(params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null");
        }
        return pot;
      },
      nb::arg("parameters"),
      "Construct a Potential from Parameters.potential");
}

} // namespace eonc::pybind
