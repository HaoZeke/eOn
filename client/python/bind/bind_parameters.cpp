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
      // --- Debug ---
      .def_prop_rw(
          "write_movies",
          [](const eonc::Parameters &s) { return s.debug_options.write_movies; },
          [](eonc::Parameters &s, bool v) { s.debug_options.write_movies = v; })
      // --- Metatomic ---
      .def_prop_rw(
          "metatomic_model_path",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.model_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.model_path = v;
          })
      .def_prop_rw(
          "metatomic_device",
          [](const eonc::Parameters &s) { return s.metatomic_options.device; },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.device = v;
          })
      .def_prop_rw(
          "metatomic_length_unit",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.length_unit;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.metatomic_options.length_unit = v;
          })
      .def_prop_rw(
          "metatomic_deterministic",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.deterministic;
          },
          [](eonc::Parameters &s, bool v) {
            s.metatomic_options.deterministic = v;
          })
      .def_prop_rw(
          "metatomic_deterministic_strict",
          [](const eonc::Parameters &s) {
            return s.metatomic_options.deterministic_strict;
          },
          [](eonc::Parameters &s, bool v) {
            s.metatomic_options.deterministic_strict = v;
          })
      // --- NEB ---
      .def_prop_rw(
          "neb_images",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.neb_options.image_count);
          },
          [](eonc::Parameters &s, long v) {
            s.neb_options.image_count = v;
          })
      .def_prop_rw(
          "neb_max_iterations",
          [](const eonc::Parameters &s) {
            return static_cast<long>(s.neb_options.max_iterations);
          },
          [](eonc::Parameters &s, long v) {
            s.neb_options.max_iterations = v;
          })
      .def_prop_rw(
          "neb_force_tolerance",
          [](const eonc::Parameters &s) {
            return s.neb_options.force_tolerance;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.force_tolerance = v;
          })
      .def_prop_rw(
          "neb_init_method",
          [](const eonc::Parameters &s) {
            return s.neb_options.initialization.method;
          },
          [](eonc::Parameters &s, eonc::NEBInit v) {
            s.neb_options.initialization.method = v;
          })
      .def_prop_rw(
          "neb_initial_path",
          [](const eonc::Parameters &s) {
            return s.neb_options.initialization.input_path;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.neb_options.initialization.input_path = v;
          })
      .def_prop_rw(
          "neb_minimize_endpoints",
          [](const eonc::Parameters &s) {
            return s.neb_options.endpoints.minimize;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.endpoints.minimize = v;
          })
      .def_prop_rw(
          "neb_climbing_image",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.enabled = v;
          })
      .def_prop_rw(
          "neb_climbing_converged_only",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.converged_only;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.converged_only = v;
          })
      .def_prop_rw(
          "neb_ci_after",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_after_rel",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_energy_weighted",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.spring.weighting.enabled = v;
          })
      .def_prop_rw(
          "neb_ew_ksp_min",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.k_min;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.spring.weighting.k_min = v;
          })
      .def_prop_rw(
          "neb_ew_ksp_max",
          [](const eonc::Parameters &s) {
            return s.neb_options.spring.weighting.k_max;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.spring.weighting.k_max = v;
          })
      .def_prop_rw(
          "neb_ci_mmf",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.use_mmf;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.climbing_image.ocineb.use_mmf = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.trigger_force;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.trigger_force = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_after_rel",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.trigger_factor;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.trigger_factor = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_nsteps",
          [](const eonc::Parameters &s) {
            return static_cast<long>(
                s.neb_options.climbing_image.ocineb.max_steps);
          },
          [](eonc::Parameters &s, long v) {
            s.neb_options.climbing_image.ocineb.max_steps = v;
          })
      .def_prop_rw(
          "neb_ci_mmf_angle",
          [](const eonc::Parameters &s) {
            return s.neb_options.climbing_image.ocineb.angle_tol;
          },
          [](eonc::Parameters &s, double v) {
            s.neb_options.climbing_image.ocineb.angle_tol = v;
          })
      .def_prop_rw(
          "neb_mmf_peaks",
          [](const eonc::Parameters &s) {
            return s.neb_options.mmf_peaks.enabled;
          },
          [](eonc::Parameters &s, bool v) {
            s.neb_options.mmf_peaks.enabled = v;
          })
      // --- RgpotPot (incl. backend=metatomic → libmetatomic_engine) ---
      .def_prop_rw(
          "rgpot_backend",
          [](const eonc::Parameters &s) { return s.rgpot_options.backend; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.backend = v;
          })
      .def_prop_rw(
          "rgpot_engine_path",
          [](const eonc::Parameters &s) { return s.rgpot_options.engine_path; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.engine_path = v;
          })
      .def_prop_rw(
          "rgpot_engine_library",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.engine_library;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.engine_library = v;
          })
      .def_prop_rw(
          "rgpot_model_path",
          [](const eonc::Parameters &s) { return s.rgpot_options.model_path; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.model_path = v;
          })
      .def_prop_rw(
          "rgpot_device",
          [](const eonc::Parameters &s) { return s.rgpot_options.device; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.device = v;
          })
      .def_prop_rw(
          "rgpot_length_unit",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.length_unit;
          },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.length_unit = v;
          })
      .def_prop_rw(
          "rgpot_check_consistency",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.check_consistency;
          },
          [](eonc::Parameters &s, bool v) {
            s.rgpot_options.check_consistency = v;
          })
      .def_prop_rw(
          "rgpot_uncertainty_threshold",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.uncertainty_threshold;
          },
          [](eonc::Parameters &s, double v) {
            s.rgpot_options.uncertainty_threshold = v;
          })
      .def_prop_rw(
          "rgpot_torch_determinism_strict",
          [](const eonc::Parameters &s) {
            return s.rgpot_options.torch_determinism_strict;
          },
          [](eonc::Parameters &s, bool v) {
            s.rgpot_options.torch_determinism_strict = v;
          })
      .def_prop_rw(
          "rgpot_basis",
          [](const eonc::Parameters &s) { return s.rgpot_options.basis; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.basis = v;
          })
      .def_prop_rw(
          "rgpot_theory",
          [](const eonc::Parameters &s) { return s.rgpot_options.theory; },
          [](eonc::Parameters &s, const std::string &v) {
            s.rgpot_options.theory = v;
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
