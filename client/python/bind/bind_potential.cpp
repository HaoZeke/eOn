#include "Parameters.h"
#include "Potential.h"
#include "Eigen.h"
#include <nanobind/ndarray.h>
#include <algorithm>

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
      .def(
          "get_ef",
          [](eonc::Potential &self,
             nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu> pos,
             nb::ndarray<nb::numpy, int64_t, nb::c_contig, nb::device::cpu> z,
             nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu> box) {
            if (pos.ndim() != 2 || pos.shape(1) != 3)
              throw std::invalid_argument("pos must be (n,3)");
            if (z.ndim() != 1 || z.shape(0) != pos.shape(0))
              throw std::invalid_argument("z must be (n,)");
            if (box.ndim() != 2 || box.shape(0) != 3 || box.shape(1) != 3)
              throw std::invalid_argument("box must be (3,3)");
            const auto n = static_cast<Eigen::Index>(pos.shape(0));
            AtomMatrix R(n, 3);
            std::copy(pos.data(), pos.data() + n * 3, R.data());
            VectorXi atmnrs(n);
            for (Eigen::Index i = 0; i < n; ++i)
              atmnrs(i) = static_cast<int>(z.data()[i]);
            Matrix3d cell;
            std::copy(box.data(), box.data() + 9, cell.data());
            auto [energy, forces] = self.get_ef(R, atmnrs, cell);
            // return (energy, forces ndarray)
            const size_t rows = static_cast<size_t>(forces.rows());
            const size_t cols = 3;
            double *buf = new double[rows * cols];
            std::copy(forces.data(), forces.data() + rows * cols, buf);
            nb::capsule owner(buf, [](void *p) noexcept {
              delete[] static_cast<double *>(p);
            });
            auto f_np = nb::ndarray<nb::numpy, double, nb::c_contig>(
                buf, {rows, cols}, owner);
            return nb::make_tuple(energy, f_np);
          },
          nb::arg("positions"), nb::arg("atomic_numbers"), nb::arg("cell"),
          "Evaluate energy and forces: returns (energy, forces (n,3))")
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
