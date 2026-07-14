#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "eigen_numpy.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
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

      // --- positions ---
      .def_prop_rw(
          "positions",
          [](const Matter &self) {
            return matrix_to_numpy(self.getPositions());
          },
          [](Matter &self, const NpF64 &arr) {
            self.setPositions(atom_matrix_from_numpy(arr));
          },
          "Cartesian positions, shape (n, 3) float64")
      .def_prop_ro(
          "positions_free",
          [](const Matter &self) {
            return matrix_to_numpy(self.getPositionsFree());
          })
      .def("set_positions_free",
           [](Matter &self, const NpF64 &arr) {
             self.setPositionsFree(atom_matrix_from_numpy(arr));
           },
           nb::arg("positions"))

      // --- cell ---
      .def_prop_rw(
          "cell",
          [](const Matter &self) { return matrix_to_numpy(self.getCell()); },
          [](Matter &self, const NpF64 &arr) {
            self.setCell(matrix3_from_numpy(arr));
          },
          "Cell matrix (3, 3) row-vector lattice")

      // --- velocities ---
      .def_prop_rw(
          "velocities",
          [](const Matter &self) {
            return matrix_to_numpy(self.getVelocities());
          },
          [](Matter &self, const NpF64 &arr) {
            self.setVelocities(atom_matrix_from_numpy(arr));
          })

      // --- forces (const getters; may evaluate potential) ---
      .def_prop_ro(
          "forces",
          [](const Matter &self) {
            return matrix_to_numpy(self.getForces());
          },
          "Forces with fixed atoms zeroed as appropriate")
      .def_prop_ro(
          "forces_raw",
          [](const Matter &self) {
            return matrix_to_numpy(self.getForcesRaw());
          })
      .def_prop_ro(
          "forces_free",
          [](const Matter &self) {
            return matrix_to_numpy(self.getForcesFree());
          })
      .def("set_forces",
           [](Matter &self, const NpF64 &arr) {
             self.setForces(atom_matrix_from_numpy(arr));
           },
           nb::arg("forces"))

      // --- masses / Z / fixed ---
      .def_prop_rw(
          "masses",
          [](const Matter &self) {
            return vector_to_numpy(self.getMasses());
          },
          [](Matter &self, const NpF64 &arr) {
            self.setMasses(vector_from_numpy(arr));
          })
      .def_prop_rw(
          "atomic_numbers",
          [](const Matter &self) {
            return vectori_to_numpy(self.getAtomicNrs());
          },
          [](Matter &self, const NpI64 &arr) {
            self.setAtomicNrs(vectori_from_numpy_i64(arr));
          },
          "Atomic numbers (Z), int64 array length n")
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
          [](const Matter &self) {
            const long n = self.numberOfAtoms();
            VectorXi v(n);
            for (long i = 0; i < n; ++i) {
              v(i) = self.getFixed(i);
            }
            return vectori_to_numpy(v);
          },
          [](Matter &self, const NpI64 &arr) {
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "fixed mask must be 1-d int array of length n_atoms");
            }
            for (size_t i = 0; i < arr.shape(0); ++i) {
              self.setFixed(static_cast<long>(i),
                            static_cast<int>(arr.data()[i]) ? 1 : 0);
            }
          },
          "Fixed flags (1=fixed, 0=free), length n")

      // --- free mask ---
      .def_prop_ro(
          "free_mask",
          [](const Matter &self) {
            return matrix_to_numpy(self.getFree());
          },
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
          nb::arg("diff"), "Minimum-image map for displacement array (n,3)")
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

      // --- relax ---
      .def(
          "relax",
          [](Matter &self, bool quiet, bool write_movie, bool checkpoint,
             const std::string &prefix_movie,
             const std::string &prefix_checkpoint) {
            return self.relax(quiet, write_movie, checkpoint, prefix_movie,
                              prefix_checkpoint);
          },
          nb::arg("quiet") = false, nb::arg("write_movie") = false,
          nb::arg("checkpoint") = false, nb::arg("prefix_movie") = "",
          nb::arg("prefix_checkpoint") = "",
          "Run client minimizer on this Matter; returns True if converged")

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
