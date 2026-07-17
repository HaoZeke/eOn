#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "eigen_numpy.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/array.h>

#include <memory>
#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_matter(nb::module_ &m) {
  using eonc::Matter;
  using eonc::Parameters;
  using eonc::Potential;

  nb::class_<Matter>(m, "Matter",
                     "Atomic structure + potential (eOn C++ client core)")
      .def(nb::init<std::shared_ptr<Potential>, const Parameters &>(),
           nb::arg("potential"), nb::arg("parameters"),
           nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(),
           "Create empty Matter bound to a potential and parameters "
           "(Parameters must outlive Matter)")
      .def(nb::init<const Matter &>(), nb::arg("other"), "Copy construct")
      .def("resize", &Matter::resize, nb::arg("n_atoms"))
      .def_prop_ro("n_atoms", &Matter::numberOfAtoms)
      .def_prop_ro("n_free", &Matter::numberOfFreeAtoms)
      .def_prop_ro("n_fixed", &Matter::numberOfFixedAtoms)

      // --- positions (zero-copy view of Matter storage) ---
      // Getter: NumPy view into Eigen row-major positions (n,3).
      // Setter: bulk memcpy + PBC + dirty flag (no temporary AtomMatrix).
      // Write-through on the view: call mark_geometry_dirty() after mutate.
      .def_prop_rw(
          "positions",
          [](Matter &self) {
            return view_n3(self.positionsData(), self.numberOfAtoms());
          },
          [](Matter &self, const NpF64 &arr) {
            const long n = self.numberOfAtoms();
            const double *p = require_n3(arr, n);
            self.assignPositions(p);
          },
          nb::rv_policy::reference_internal,
          "Cartesian positions (n,3) float64 — zero-copy view of C++ storage")
      .def("mark_geometry_dirty", &Matter::markGeometryDirty,
           "After in-place mutation of positions/cell views, force recompute")
      .def_prop_ro(
          "positions_free",
          [](const Matter &self) {
            return matrix_to_numpy(self.getPositionsFree());
          },
          nb::rv_policy::move)
      .def("set_positions_free",
           [](Matter &self, const NpF64 &arr) {
             self.setPositionsFree(atom_matrix_from_numpy(arr));
           },
           nb::arg("positions"))

      // --- cell (zero-copy 3x3 view) ---
      .def_prop_rw(
          "cell",
          [](Matter &self) { return view_33(self.cellData()); },
          [](Matter &self, const NpF64 &arr) {
            if (arr.ndim() != 2 || arr.shape(0) != 3 || arr.shape(1) != 3) {
              throw std::invalid_argument("cell must be float64 (3,3)");
            }
            self.assignCell(arr.data());
          },
          nb::rv_policy::reference_internal,
          "Cell matrix (3,3) row-major — zero-copy view")

      // --- velocities ---
      .def_prop_rw(
          "velocities",
          [](Matter &self) {
            return view_n3(self.velocitiesData(), self.numberOfAtoms());
          },
          [](Matter &self, const NpF64 &arr) {
            const long n = self.numberOfAtoms();
            self.assignVelocities(require_n3(arr, n));
          },
          nb::rv_policy::reference_internal)

      // --- forces (zero-copy into cached force storage after evaluate) ---
      .def_prop_ro(
          "forces",
          [](Matter &self) {
            const AtomMatrix &f = self.getForces();
            return view_n3(const_cast<double *>(f.data()),
                           self.numberOfAtoms());
          },
          nb::rv_policy::reference_internal,
          "Forces with fixed atoms zeroed (zero-copy view of cache)")
      .def_prop_ro(
          "forces_raw",
          [](Matter &self) {
            const AtomMatrix &f = self.getForcesRaw();
            return view_n3(const_cast<double *>(f.data()),
                           self.numberOfAtoms());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "forces_free",
          [](const Matter &self) {
            return matrix_to_numpy(self.getForcesFree());
          },
          nb::rv_policy::move)
      .def("set_forces",
           [](Matter &self, const NpF64 &arr) {
             self.setForces(atom_matrix_from_numpy(arr));
           },
           nb::arg("forces"))

      // --- masses / Z / fixed (zero-copy where layout allows) ---
      .def_prop_rw(
          "masses",
          [](Matter &self) {
            return view_n(self.massesData(), self.numberOfAtoms());
          },
          [](Matter &self, const NpF64 &arr) {
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "masses must be float64 length n_atoms");
            }
            self.assignMasses(arr.data());
          },
          nb::rv_policy::reference_internal)
      .def_prop_rw(
          "atomic_numbers",
          [](Matter &self) {
            // int32 view of VectorXi storage (portable on LP64)
            return view_n_i32(self.atomicNrsData(), self.numberOfAtoms());
          },
          [](Matter &self, nb::object obj) {
            // Accept int32 or int64 from NumPy
            auto np = nb::module_::import_("numpy");
            nb::object cont = np.attr("ascontiguousarray")(
                obj, nb::arg("dtype") = "int32");
            auto arr = nb::cast<NpI32>(cont);
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "atomic_numbers must be length n_atoms");
            }
            self.assignAtomicNrs(
                reinterpret_cast<const int *>(arr.data()));
          },
          nb::rv_policy::reference_internal,
          "Atomic numbers (Z) — zero-copy int32 view of C++ storage")
      .def(
          "get_fixed",
          [](const Matter &self, long atom) { return self.getFixed(atom); },
          nb::arg("atom"))
      .def(
          "set_fixed",
          [](Matter &self, long atom, int fixed) {
            self.setFixed(atom, fixed);
          },
          nb::arg("atom"), nb::arg("fixed"))
      .def_prop_rw(
          "fixed",
          [](Matter &self) {
            return view_n_i32(self.isFixedData(), self.numberOfAtoms());
          },
          [](Matter &self, nb::object obj) {
            auto np = nb::module_::import_("numpy");
            nb::object cont = np.attr("ascontiguousarray")(
                obj, nb::arg("dtype") = "int32");
            auto arr = nb::cast<NpI32>(cont);
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "fixed mask must be 1-d int array of length n_atoms");
            }
            self.assignIsFixed(reinterpret_cast<const int *>(arr.data()));
          },
          nb::rv_policy::reference_internal,
          "Fixed flags (1=fixed, 0=free) — zero-copy int32 view")

      // --- free mask ---
      .def_prop_ro(
          "free_mask",
          [](const Matter &self) {
            return matrix_to_numpy(self.getFree());
          },
          nb::rv_policy::move,
          "Nx3 free mask (1 free, 0 fixed)")

      // --- energy / mechanics ---
      .def_prop_ro("potential_energy", &Matter::getPotentialEnergy)
      .def_prop_ro("kinetic_energy", &Matter::getKineticEnergy)
      .def_prop_ro("mechanical_energy", &Matter::getMechanicalEnergy)
      .def_prop_ro("energy_variance", &Matter::getEnergyVariance)
      .def_prop_ro("max_force", &Matter::maxForce)
      .def_prop_ro("force_calls", &Matter::getForceCalls)
      .def("reset_force_calls", &Matter::resetForceCalls)
      .def_prop_ro("needs_force_update", &Matter::needsForceUpdate)
      .def_prop_ro("potential_calls", &Matter::getPotentialCalls)

      // --- PBC ---
      .def_prop_rw("periodic", &Matter::getPeriodic, &Matter::setPeriodic)
      .def_prop_rw("pbc_convention", &Matter::getPbcConvention,
                   &Matter::setPbcConvention)
      .def(
          "pbc",
          [](const Matter &self, const NpF64 &diff) {
            return matrix_to_numpy(self.pbc(atom_matrix_from_numpy(diff)));
          },
          nb::arg("diff"), nb::rv_policy::move,
          "Minimum-image map for displacement array (n,3)")
      .def("apply_periodic_boundary_if_enabled",
           &Matter::applyPeriodicBoundaryIfEnabled)

      // --- potential handle ---
      .def("set_potential", &Matter::setPotential, nb::arg("potential"))
      .def("get_potential", &Matter::getPotential)

      // --- geometry helpers ---
      .def("distance",
           nb::overload_cast<long, long>(&Matter::distance, nb::const_),
           nb::arg("i"), nb::arg("j"))
      .def(
          "distance_to_matter",
          [](const Matter &self, const Matter &other, long index) {
            return self.distance(other, index);
          },
          nb::arg("other"), nb::arg("index"))
      .def("distance_to", &Matter::distanceTo, nb::arg("other"))
      .def("per_atom_norm", &Matter::perAtomNorm, nb::arg("other"))
      .def("compare", &Matter::compare, nb::arg("other"),
           nb::arg("indistinguishable") = false)

      // --- relax (default non-mutating; inplace=True mutates self) ---
      .def(
          "relax",
          [](Matter &self, bool inplace, bool quiet, bool write_movie,
             bool checkpoint, const std::string &prefix_movie,
             const std::string &prefix_checkpoint) -> nb::object {
            if (inplace) {
              bool ok = false;
              {
                nb::gil_scoped_release release;
                ok = self.relax(quiet, write_movie, checkpoint, prefix_movie,
                                prefix_checkpoint);
              }
              return nb::make_tuple(
                  nb::cast(self, nb::rv_policy::reference), ok);
            }
            Matter out(self);
            bool ok = false;
            {
              nb::gil_scoped_release release;
              ok = out.relax(quiet, write_movie, checkpoint, prefix_movie,
                             prefix_checkpoint);
            }
            return nb::make_tuple(std::move(out), ok);
          },
          nb::arg("inplace") = false, nb::arg("quiet") = false,
          nb::arg("write_movie") = false, nb::arg("checkpoint") = false,
          nb::arg("prefix_movie") = "", nb::arg("prefix_checkpoint") = "",
          "Minimize. Default: copy then relax (self unchanged). "
          "inplace=True mutates self. Returns (Matter, converged: bool).")

      // --- I/O (same ConFileIO / readcon path as CLI client) ---
      .def(
          "con2matter",
          [](Matter &self, const std::string &path) {
            return self.con2matter(path);
          },
          nb::arg("path"), "Load .con file into this Matter")
      .def(
          "matter2con",
          [](Matter &self, const std::string &path, bool append) {
            return self.matter2con(path, append, nullptr);
          },
          nb::arg("path"), nb::arg("append") = false, "Write Matter to .con")
      .def(
          "convel2matter",
          [](Matter &self, const std::string &path) {
            return self.convel2matter(path);
          },
          nb::arg("path"))
      .def(
          "matter2convel",
          [](Matter &self, const std::string &path) {
            return self.matter2convel(path);
          },
          nb::arg("path"))
      .def(
          "matter2xyz",
          [](Matter &self, const std::string &path, bool append) {
            return self.matter2xyz(path, append);
          },
          nb::arg("path"), nb::arg("append") = false)
      .def(
          "write_tibble",
          [](Matter &self, const std::string &path) {
            return self.writeTibble(path);
          },
          nb::arg("path"))

      // --- atom index / header (bindings-friendly) ---
      .def("get_atom_index", &Matter::getAtomIndex, nb::arg("atom"))
      .def("set_atom_index", &Matter::setAtomIndex, nb::arg("atom"),
           nb::arg("index"))
      .def_prop_ro(
          "header_con",
          [](const Matter &self) { return self.getHeaderCon(); })
      .def(
          "set_header_con_line",
          [](Matter &self, size_t i, const std::string &line) {
            self.setHeaderConLine(i, line);
          },
          nb::arg("i"), nb::arg("line"))

      .def("__len__", &Matter::numberOfAtoms)
      .def("__repr__", [](const Matter &self) {
        return "<Matter n=" + std::to_string(self.numberOfAtoms()) +
               " E=" + std::to_string(self.getPotentialEnergy()) + ">";
      });
}

} // namespace eonc::pybind
