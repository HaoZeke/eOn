#include "Parameters.h"
#include "Potential.h"
#include "BaseStructures.h"
#include "Eigen.h"
#include <nanobind/ndarray.h>
#include <algorithm>
#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <magic_enum/magic_enum.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;


namespace {

/// Live ASE calculator as Potential (seamless Python path; no ASE_POT script).
///
/// Hot path design (force is called every MD/opt step):
/// - Reuse one ``ase.Atoms`` + the same calculator (no per-call construction).
/// - Wrap C++ ``R`` / ``box`` / ``Z`` as NumPy **views** (zero alloc / zero
///   Python list materialization). ASE then bulk-copies positions into its
///   own buffer via ``positions[:] = view`` (one memcpy, not O(N) Python).
/// - One ``calc.calculate(..., ['energy','forces'])`` (not energy then forces).
/// - Forces leave via ``memcpy`` from a contiguous float64 ndarray.
class AseCalcPotential final : public eonc::Potential {
  nb::object calc_;
  nb::object atoms_;          // cached ase.Atoms
  nb::object np_;             // numpy module
  nb::object props_;          // ['energy', 'forces']
  nb::object all_changes_;    // ase.calculators.calculator.all_changes
  std::vector<int> z_cache_;  // last atomic numbers (detect composition change)
  long n_cached_{-1};
  bool modules_ready_{false};

  void ensure_modules() {
    if (modules_ready_)
      return;
    np_ = nb::module_::import_("numpy");
    auto calc_mod = nb::module_::import_("ase.calculators.calculator");
    all_changes_ = calc_mod.attr("all_changes");
    props_ = nb::list();
    props_.attr("append")("energy");
    props_.attr("append")("forces");
    modules_ready_ = true;
  }

  /// NumPy view of external C++ buffer (no copy; valid only for this force()).
  static nb::object view_f64(const double *data, size_t d0, size_t d1) {
    // ndarray without capsule owner → view; caller must not outlive data.
    nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu> arr(
        const_cast<double *>(data), {d0, d1});
    return nb::cast(arr, nb::rv_policy::reference);
  }

  static nb::object view_i32(const int *data, size_t n) {
    static_assert(sizeof(int) == 4, "AseCalcPotential Z view assumes 32-bit int");
    nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu> arr(
        const_cast<int32_t *>(reinterpret_cast<const int32_t *>(data)), {n});
    return nb::cast(arr, nb::rv_policy::reference);
  }

  void ensure_atoms(long nAtoms, const int *atomicNrs) {
    const size_t n = static_cast<size_t>(nAtoms);
    const bool size_changed = (n_cached_ != nAtoms);
    bool z_changed = size_changed;
    if (!z_changed) {
      z_changed = !std::equal(z_cache_.begin(), z_cache_.end(), atomicNrs);
    }
    if (!size_changed && !z_changed && atoms_.is_valid() && !atoms_.is_none())
      return;

    z_cache_.assign(atomicNrs, atomicNrs + n);
    n_cached_ = nAtoms;

    // Numbers: view into z_cache_ (stable for lifetime of this Potential call
    // sequence; reallocated only on composition change).
    nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu> z_arr(
        z_cache_.data(), {n});
    nb::object numbers = nb::cast(z_arr, nb::rv_policy::reference);

    // Placeholder geometry; force() overwrites via zero-copy views each call.
    nb::object zeros = np_.attr("zeros");
    nb::object positions =
        zeros(nb::make_tuple(nAtoms, 3), nb::arg("dtype") = "float64");
    nb::object cell = np_.attr("eye")(3, nb::arg("dtype") = "float64");

    nb::object Atoms = nb::module_::import_("ase").attr("Atoms");
    atoms_ =
        Atoms(nb::arg("numbers") = numbers, nb::arg("positions") = positions,
              nb::arg("cell") = cell, nb::arg("pbc") = true);
    atoms_.attr("calc") = calc_;
  }

  static void copy_forces_out(nb::object forces_obj, long nAtoms, double *F) {
    const size_t n3 = static_cast<size_t>(nAtoms) * 3u;
    // Prefer a direct buffer view of the calculator result (no Python loop).
    try {
      auto farr =
          nb::cast<nb::ndarray<nb::numpy, double, nb::device::cpu>>(forces_obj);
      if (farr.ndim() == 2 &&
          static_cast<long>(farr.shape(0)) == nAtoms && farr.shape(1) == 3) {
        // C-contiguous (n,3): stride0==3, stride1==1 in elements
        if (farr.stride(0) == 3 && farr.stride(1) == 1) {
          std::memcpy(F, farr.data(), n3 * sizeof(double));
          return;
        }
        // Generic strided read (still bulk C++, no Python scalars)
        const double *base = farr.data();
        const auto s0 = static_cast<size_t>(farr.stride(0));
        const auto s1 = static_cast<size_t>(farr.stride(1));
        for (long i = 0; i < nAtoms; ++i)
          for (long j = 0; j < 3; ++j)
            F[i * 3 + j] =
                base[static_cast<size_t>(i) * s0 + static_cast<size_t>(j) * s1];
        return;
      }
    } catch (const nb::cast_error &) {
      // fall through to ascontiguousarray
    }
    auto np = nb::module_::import_("numpy");
    nb::object cont = np.attr("ascontiguousarray")(
        forces_obj, nb::arg("dtype") = "float64");
    auto farr =
        nb::cast<nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>>(
            cont);
    if (farr.ndim() != 2 || static_cast<long>(farr.shape(0)) != nAtoms ||
        farr.shape(1) != 3) {
      throw std::runtime_error(
          "ASE forces must be shape (n_atoms, 3) float64");
    }
    std::memcpy(F, farr.data(), n3 * sizeof(double));
  }

public:
  explicit AseCalcPotential(nb::object calc)
      : eonc::Potential(eonc::PotType::ASE_POT), calc_(std::move(calc)) {
    if (!calc_.is_valid() || calc_.is_none()) {
      throw std::invalid_argument(
          "make_potential_from_ase: calculator is None");
    }
  }

  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double * /*variance*/, const double *box) override {
    if (nAtoms <= 0) {
      *U = 0.0;
      return;
    }
    if (!R || !atomicNrs || !F || !U || !box) {
      throw std::invalid_argument("make_potential_from_ase force: null buffer");
    }
    nb::gil_scoped_acquire gil;
    try {
      ensure_modules();
      ensure_atoms(nAtoms, atomicNrs);

      // Zero-copy views of caller-owned buffers (valid for this call only).
      nb::object pos_view =
          view_f64(R, static_cast<size_t>(nAtoms), size_t{3});
      nb::object cell_view = view_f64(box, size_t{3}, size_t{3});

      // Bulk assign into ASE arrays (one C-level memcpy each via NumPy, not
      // per-atom Python). Views above avoid building list/asarray trees.
      atoms_.attr("set_positions")(pos_view);
      atoms_.attr("set_cell")(cell_view, nb::arg("scale_atoms") = false);

      // Single calculator evaluation for energy + forces.
      nb::object forces_obj;
      try {
        calc_.attr("calculate")(atoms_, props_, all_changes_);
        nb::object results = calc_.attr("results");
        *U = nb::cast<double>(results["energy"]);
        forces_obj = results["forces"];
      } catch (const nb::python_error &) {
        // Fallback for calculators without calculate()/results.
        *U = nb::cast<double>(atoms_.attr("get_potential_energy")());
        forces_obj = atoms_.attr("get_forces")();
      }
      copy_forces_out(forces_obj, nAtoms, F);
    } catch (const nb::python_error &e) {
      throw std::runtime_error(std::string("ASE calculator Python error: ") +
                               e.what());
    } catch (const std::exception &e) {
      throw std::runtime_error(std::string("ASE calculator error: ") + e.what());
    }
  }

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }
};

} // namespace

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

  // Seamless ASE Calculator → Potential (runtime ASE import; no -Dwith_ase).
  m.def(
      "make_potential_from_ase",
      [](nb::object calculator) -> std::shared_ptr<eonc::Potential> {
        return std::make_shared<AseCalcPotential>(std::move(calculator));
      },
      nb::arg("calculator"),
      "Wrap a live ASE Calculator as eOn Potential (no ASE_POT script path).");
}

} // namespace eonc::pybind
