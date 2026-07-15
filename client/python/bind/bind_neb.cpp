/*
** NudgedElasticBand + path init helpers — first-class NEB surface for pyeonclient.
*/
#include "Matter.h"
#include "NEBInitialPaths.hpp"
#include "NEBSplineExtrema.h"
#include "NudgedElasticBand.h"
#include "Parameters.h"
#include "Potential.h"
#include "PotRegistry.h"
#include "eigen_numpy.hpp"

#include <magic_enum/magic_enum.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_neb(nb::module_ &m) {
  using eonc::Matter;
  using eonc::NudgedElasticBand;
  using eonc::Parameters;
  using eonc::Potential;

  nb::enum_<NudgedElasticBand::NEBStatus>(m, "NEBStatus")
      .value("GOOD", NudgedElasticBand::NEBStatus::GOOD)
      .value("INIT", NudgedElasticBand::NEBStatus::INIT)
      .value("BAD_MAX_ITERATIONS",
             NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS)
      .value("RUNNING", NudgedElasticBand::NEBStatus::RUNNING)
      .value("MAX_UNCERTAINTY", NudgedElasticBand::NEBStatus::MAX_UNCERTAINTY)
      .export_values();

  nb::class_<NudgedElasticBand>(m, "NudgedElasticBand",
                                "NEB band: path images, compute, forces, extrema")
      // --- construct from endpoints (linear / IDPP / FILE init via Parameters) ---
      .def(
          "__init__",
          [](NudgedElasticBand *self, std::shared_ptr<Matter> initial,
             std::shared_ptr<Matter> final_state, const Parameters &params,
             std::shared_ptr<Potential> pot) {
            new (self)
                NudgedElasticBand(std::move(initial), std::move(final_state),
                                  params, std::move(pot));
          },
          nb::arg("initial"), nb::arg("final"), nb::arg("parameters"),
          nb::arg("potential"),
          nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
          nb::keep_alive<1, 4>(), nb::keep_alive<1, 5>(),
          "Build NEB from reactant/product Matter; path init follows "
          "Parameters.neb_options.initialization")
      // --- construct from explicit path vector (copied into band) ---
      .def(
          "__init__",
          [](NudgedElasticBand *self, const std::vector<Matter> &path,
             const Parameters &params, std::shared_ptr<Potential> pot) {
            new (self) NudgedElasticBand(path, params, std::move(pot));
          },
          nb::arg("path"), nb::arg("parameters"), nb::arg("potential"),
          nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>(),
          "Build NEB from a list of Matter frames (length = intermediates + 2)")

      // --- algorithm steps ---
      .def("compute", &NudgedElasticBand::compute,
           "Run NEB optimization (energy-weighted, CI, OCI per Parameters)")
      .def("update_forces",
           nb::overload_cast<bool>(&NudgedElasticBand::updateForces),
           nb::arg("ci_active"), "Recompute projected forces; ci_active toggles CI")
      .def("update_forces", nb::overload_cast<>(&NudgedElasticBand::updateForces),
           "Recompute projected forces with current CI flag")
      .def("set_ci_enabled", &NudgedElasticBand::setCIEnabled, nb::arg("enabled"))
      .def("convergence_force", &NudgedElasticBand::convergenceForce)
      .def("find_extrema", &NudgedElasticBand::findExtrema)
      .def("print_image_data", &NudgedElasticBand::printImageData,
           nb::arg("write_to_file") = false, nb::arg("idx") = 0,
           "Print (and optionally write neb.dat) image reaction coordinates")

      // --- status / geometry ---
      .def_prop_ro("status", [](NudgedElasticBand &n) { return n.getStatus(); })
      .def_prop_ro("num_images",
                   [](const NudgedElasticBand &n) { return n.numImages; })
      .def_prop_ro("climbing_image",
                   [](const NudgedElasticBand &n) { return n.climbingImage; })
      .def_prop_ro("max_energy_image",
                   [](const NudgedElasticBand &n) {
                     return static_cast<long>(n.maxEnergyImage);
                   })
      .def_prop_ro("num_extrema",
                   [](const NudgedElasticBand &n) { return n.numExtrema; })
      .def_prop_ro("energy_reference",
                   [](const NudgedElasticBand &n) { return n.E_ref; })
      .def_prop_ro("ksp", [](const NudgedElasticBand &n) { return n.ksp; })
      .def_prop_ro(
          "n_path",
          [](const NudgedElasticBand &n) {
            return static_cast<long>(n.path.size());
          },
          "len(path) including endpoints (= num_images + 2)")

      // --- path images (shared_ptr Matter) ---
      .def(
          "image",
          [](NudgedElasticBand &n, long i) {
            if (i < 0 || static_cast<size_t>(i) >= n.path.size())
              throw std::out_of_range("NEB image index out of range");
            return n.path[static_cast<size_t>(i)];
          },
          nb::arg("i"), nb::rv_policy::reference_internal,
          "Shared Matter at path index i (0=reactant, -1/last=product)")
      .def(
          "path_images",
          [](NudgedElasticBand &n) { return n.path; },
          nb::rv_policy::reference_internal,
          "Full path as list of shared Matter (endpoints + intermediates)")
      .def(
          "image_energy",
          [](const NudgedElasticBand &n, long i) {
            if (i < 0 || static_cast<size_t>(i) >= n.path.size())
              throw std::out_of_range("NEB image index out of range");
            return n.path[static_cast<size_t>(i)]->getPotentialEnergy();
          },
          nb::arg("i"))
      .def(
          "image_force_norm",
          [](const NudgedElasticBand &n, long i) {
            if (i < 0 || static_cast<size_t>(i) >= n.path.size())
              throw std::out_of_range("NEB image index out of range");
            return n.path[static_cast<size_t>(i)]->getForces().norm();
          },
          nb::arg("i"))
      .def(
          "projected_force_norm",
          [](const NudgedElasticBand &n, long i) {
            if (i < 1 || i > n.numImages)
              return 0.0;
            if (static_cast<size_t>(i) >= n.projectedForce.size() ||
                !n.projectedForce[static_cast<size_t>(i)])
              return 0.0;
            return n.projectedForce[static_cast<size_t>(i)]->norm();
          },
          nb::arg("i"))
      .def_prop_ro(
          "extremum_positions",
          [](const NudgedElasticBand &n) { return n.extremumPosition; })
      .def_prop_ro(
          "extremum_energies",
          [](const NudgedElasticBand &n) { return n.extremumEnergy; })
      .def_prop_ro(
          "extremum_curvatures",
          [](const NudgedElasticBand &n) { return n.extremumCurvature; })
      .def("__repr__", [](NudgedElasticBand &n) {
        return "<NudgedElasticBand images=" + std::to_string(n.numImages) +
               " status=" +
               std::to_string(static_cast<int>(n.getStatus())) + ">";
      });

  // --- path init helpers (same as NEBJob / NudgedElasticBand internals) ---
  m.def(
      "neb_read_file_paths",
      [](const std::string &list_file) {
        auto paths = eonc::helpers::neb_paths::readFilePaths(list_file);
        std::vector<std::string> out;
        out.reserve(paths.size());
        for (const auto &p : paths)
          out.push_back(p.string());
        return out;
      },
      nb::arg("list_file"),
      "Read idppPath.dat-style list of .con paths (one path per line)");

  m.def(
      "neb_load_path_from_files",
      [](const std::vector<std::string> &files,
         std::shared_ptr<Potential> pot, const Parameters &params) {
        if (files.size() < 2)
          throw std::invalid_argument(
              "neb_load_path_from_files needs >= 2 frames");
        std::vector<Matter> path;
        path.reserve(files.size());
        for (const auto &f : files) {
          Matter m(pot, params);
          auto st = m.con2matter(f);
          if (!eonc::io::io_ok(st))
            throw std::runtime_error("failed to load NEB frame: " + f);
          path.push_back(std::move(m));
        }
        return path;
      },
      nb::arg("files"), nb::arg("potential"), nb::arg("parameters"),
      "Load list of .con files into Matter frames for NEB path constructor");

  m.def(
      "neb_linear_path",
      [](const Matter &initial, const Matter &final_state, long n_intermediate) {
        return eonc::helpers::neb_paths::linearPath(
            initial, final_state, static_cast<size_t>(n_intermediate));
      },
      nb::arg("initial"), nb::arg("final"), nb::arg("n_intermediate"),
      "Linear interpolate n_intermediate images between endpoints");

  // Full NEBJob::saveData artifact set (callable as a Python step).
  m.def(
      "neb_write_results",
      [](NudgedElasticBand &neb, const Parameters &params,
         size_t force_calls_neb) {
        std::vector<std::string> returnFiles;
        {
          std::ofstream out("results.dat", std::ios::binary);
          if (!out)
            throw std::runtime_error("neb_write_results: cannot open results.dat");
          auto status = neb.getStatus();
          out << std::format("{} termination_reason\n",
                             static_cast<int>(status));
          out << std::format("{} termination_reason_text\n",
                             magic_enum::enum_name(status));
          out << std::format(
              "{} potential_type\n",
              magic_enum::enum_name(params.potential_options.potential));
          out << std::format("{} total_force_calls\n",
                             PotRegistry::get().total_force_calls());
          out << std::format("{} force_calls_neb\n", force_calls_neb);
          out << std::format("{:f} energy_reference\n",
                             neb.path[0]->getPotentialEnergy());
          out << std::format("{} number_of_images\n", neb.numImages);
          for (long i = 0; i <= neb.numImages + 1; i++) {
            out << std::format("{:f} image{}_energy\n",
                               neb.path[static_cast<size_t>(i)]
                                       ->getPotentialEnergy() -
                                   neb.path[0]->getPotentialEnergy(),
                               i);
            out << std::format(
                "{:f} image{}_force\n",
                neb.path[static_cast<size_t>(i)]->getForces().norm(), i);
            double proj = 0.0;
            if (i >= 1 && i <= neb.numImages &&
                static_cast<size_t>(i) < neb.projectedForce.size() &&
                neb.projectedForce[static_cast<size_t>(i)])
              proj = neb.projectedForce[static_cast<size_t>(i)]->norm();
            out << std::format("{:f} image{}_projected_force\n", proj, i);
          }
          long safe_ext = std::min(
              neb.numExtrema, static_cast<long>(neb.extremumPosition.size()));
          out << std::format("{} number_of_extrema\n", safe_ext);
          for (long i = 0; i < safe_ext; i++) {
            out << std::format("{:f} extremum{}_position\n",
                               neb.extremumPosition[static_cast<size_t>(i)], i);
            out << std::format("{:f} extremum{}_energy\n",
                               neb.extremumEnergy[static_cast<size_t>(i)], i);
          }
        }
        returnFiles.emplace_back("results.dat");

        // Multi-frame path (neb.con) — required by cookbook plots
        const std::string nebFilename = "neb.con";
        if (!eonc::io::io_ok(eonc::neb::writePathCon(
                neb.path, neb.tangent, neb.eigenmode_solvers, neb.numImages,
                params.debug_options.estimate_neb_eigenvalues, nebFilename))) {
          throw std::runtime_error("neb_write_results: failed neb.con");
        }
        returnFiles.push_back(nebFilename);

        // Discrete saddle
        const std::string sp = "sp.con";
        if (!eonc::io::io_ok(
                neb.path[neb.maxEnergyImage]->matter2con(sp)))
          throw std::runtime_error("neb_write_results: failed sp.con");
        returnFiles.push_back(sp);

        // MMF peaks (same filters as NEBJob::saveData)
        if (params.neb_options.mmf_peaks.enabled && neb.numExtrema > 0) {
          int peakCount = 0;
          for (long i = 0; i < neb.numExtrema; i++) {
            double relativeEnergy =
                neb.extremumEnergy[static_cast<size_t>(i)] -
                neb.path[0]->getPotentialEnergy();
            if (!(neb.extremumCurvature[static_cast<size_t>(i)] < 0 &&
                  relativeEnergy > params.neb_options.mmf_peaks.tolerance))
              continue;
            double posFraction = neb.extremumPosition[static_cast<size_t>(i)];
            int leftIdx = static_cast<int>(std::floor(posFraction));
            double f = posFraction - leftIdx;
            if (leftIdx < 0 || leftIdx >= neb.numImages + 1)
              continue;
            Matter peakPos = eonc::helpers::neb_paths::interpolateImage(
                *neb.path[static_cast<size_t>(leftIdx)],
                *neb.path[static_cast<size_t>(leftIdx + 1)], f);
            std::string peakPosFile =
                std::format("peak{:02d}_pos.con", peakCount);
            if (!eonc::io::io_ok(peakPos.matter2con(peakPosFile)))
              throw std::runtime_error("neb_write_results: " + peakPosFile);
            returnFiles.push_back(peakPosFile);
            AtomMatrix peakMode =
                (1.0 - f) * (*neb.tangent[static_cast<size_t>(leftIdx)]) +
                f * (*neb.tangent[static_cast<size_t>(leftIdx + 1)]);
            peakMode.normalize();
            std::string peakModeFile =
                std::format("peak{:02d}_mode.dat", peakCount);
            {
              std::ofstream modeOut(peakModeFile);
              if (!modeOut)
                throw std::runtime_error("neb_write_results: " + peakModeFile);
              for (long row = 0; row < peakMode.rows(); ++row) {
                modeOut << std::format("{:12.6f} {:12.6f} {:12.6f}\n",
                                       peakMode(row, 0), peakMode(row, 1),
                                       peakMode(row, 2));
              }
            }
            returnFiles.push_back(peakModeFile);
            peakCount++;
          }
        }

        // neb.dat (+ per-iteration neb_###.dat already written during compute)
        returnFiles.emplace_back("neb.dat");
        neb.printImageData(true, std::numeric_limits<size_t>::max());
        return returnFiles;
      },
      nb::arg("neb"), nb::arg("parameters"),
      nb::arg("force_calls_neb") = 0,
      "Write full NEBJob::saveData artifact set: results.dat, neb.con, "
      "sp.con, peaks, neb.dat");

  m.def(
      "pot_registry_total_force_calls",
      []() { return PotRegistry::get().total_force_calls(); },
      "Global PotRegistry force-call counter (for force_calls_neb deltas)");
}

} // namespace eonc::pybind
