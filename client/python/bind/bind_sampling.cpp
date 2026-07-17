/*
** Matter-first sampling / long-timescale primitives with explicit inplace.
** Issues: eOn-p8em dynamics, eOn-dvf5 MC, eOn-73ji BH, eOn-m6o8 process search,
** eOn-bcd3 TAD, eOn-tw3g parallel-replica family, eOn-515s structure comparison,
** eOn-cb3z GP surrogate (gated).
*/
#include "Dynamics.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "MonteCarlo.h"
#include "Parameters.h"
#include "Potential.h"
#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

void ensure_dynamics_steps(eonc::Parameters &params, long default_steps) {
  if (params.dynamics_options.steps <= 0) {
    params.dynamics_options.steps = default_steps;
  }
  if (params.dynamics_options.time_step <= 0.0) {
    params.dynamics_options.time_step =
        params.dynamics_options.time_step_input / params.constants.timeUnit;
    if (params.dynamics_options.time_step <= 0.0)
      params.dynamics_options.time_step =
          1.0 / params.constants.timeUnit; // 1 fs
  }
}

AtomMatrix random_direction(long n_atoms, std::mt19937_64 &rng) {
  AtomMatrix m(n_atoms, 3);
  std::normal_distribution<double> n(0.0, 1.0);
  for (long i = 0; i < n_atoms; ++i)
    for (int k = 0; k < 3; ++k)
      m(i, k) = n(rng);
  // normalize free subspace roughly
  double nrm = m.norm();
  if (nrm > 1e-12)
    m /= nrm;
  return m;
}

} // namespace

void bind_sampling(nb::module_ &m) {
  using eonc::Dynamics;
  using eonc::DynamicsConfig;
  using eonc::Matter;
  using eonc::MinModeSaddleSearch;
  using eonc::MonteCarlo;
  using eonc::Parameters;
  using eonc::Potential;

  m.def(
      "built_with_gp_surrogate",
      []() {
#ifdef WITH_GP_SURROGATE
        return true;
#else
        return false;
#endif
      },
      "True if compiled with WITH_GP_SURROGATE / meson with_gp_surrogate");

  // --- Dynamics (eOn-p8em) ---
  m.def(
      "run_dynamics",
      [](std::shared_ptr<Matter> matter, Parameters params, bool inplace) {
        auto work = matter_work(matter, inplace);
        ensure_dynamics_steps(params, 10);
        Dynamics dyn(work.get(), DynamicsConfig::fromParams(params));
        {
          nb::gil_scoped_release release;
          dyn.run();
        }
        return work;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Run velocity-Verlet dynamics on Matter. Default non-mutating copy. "
      "Uses Parameters.dynamics_steps / dynamics_time_step / temperature.");

  // --- Monte Carlo (eOn-dvf5) ---
  m.def(
      "run_monte_carlo",
      [](std::shared_ptr<Matter> matter, const Parameters &params,
         bool inplace) {
        auto work = matter_work(matter, inplace);
        MonteCarlo mc(work, params);
        int steps = params.monte_carlo_options.steps;
        if (steps <= 0)
          steps = 10;
        {
          nb::gil_scoped_release release;
          mc.run(steps, params.main_options.temperature,
                 params.monte_carlo_options.step_size);
        }
        return work;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Metropolis MC on Matter. Default non-mutating. "
      "Parameters.monte_carlo_steps / monte_carlo_step_size / temperature.");

  // --- Basin hopping (eOn-73ji): Metropolis + random displace + relax ---
  m.def(
      "run_basin_hopping",
      [](std::shared_ptr<Matter> matter, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        auto work = matter_work(matter, inplace);
        // ensure pot on work matches
        work->setPotential(pot);
        long nsteps = params.basin_hopping_options.steps;
        if (nsteps <= 0)
          nsteps = 5;
        double max_disp = params.basin_hopping_options.displacement;
        double kT = params.constants.kB * params.main_options.temperature;
        std::mt19937_64 rng(
            params.main_options.randomSeed >= 0
                ? static_cast<uint64_t>(params.main_options.randomSeed)
                : 0xC0FFEEULL);
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        std::normal_distribution<double> gauss(0.0, 1.0);

        {
          nb::gil_scoped_release release;
          work->relax(/*quiet=*/true);
          double Ecur = work->getPotentialEnergy();
          Matter best(*work);
          double Ebest = Ecur;

          for (long step = 0; step < nsteps; ++step) {
            Matter trial(*work);
            AtomMatrix pos = trial.getPositions();
            for (long i = 0; i < trial.numberOfAtoms(); ++i) {
              if (trial.getFixed(i))
                continue;
              for (int k = 0; k < 3; ++k)
                pos(i, k) += max_disp * gauss(rng);
            }
            trial.setPositions(pos);
            trial.relax(/*quiet=*/true);
            double Etrial = trial.getPotentialEnergy();
            double dE = Etrial - Ecur;
            bool accept = dE <= 0.0 || uni(rng) < std::exp(-dE / kT);
            if (accept) {
              *work = trial;
              Ecur = Etrial;
              if (Ecur < Ebest) {
                best = *work;
                Ebest = Ecur;
              }
            }
          }
          *work = best;
        }
        return work;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("potential"),
      nb::arg("inplace") = false,
      "Basin-hopping style Metropolis + random free-atom displace + "
      "Matter.relax. Default non-mutating. Returns best Matter.");

  // --- Process search (eOn-m6o8): min-mode saddle compose ---
  m.def(
      "process_search",
      [](std::shared_ptr<Matter> matter, const NpF64 &mode, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        auto reactant = matter_work(matter, inplace);
        reactant->setPotential(pot);
        if (params.process_search_options.minimize_first) {
          nb::gil_scoped_release release;
          reactant->relax(/*quiet=*/true);
        }
        double E0 = reactant->getPotentialEnergy();
        // Working saddle starts as copy of reactant, then search mutates it.
        auto saddle = std::make_shared<Matter>(*reactant);
        AtomMatrix mode_am = atom_matrix_from_numpy(mode);
        // Mild displacement along mode if magnitude > 0
        double mag = params.saddle_search_options.displace_magnitude;
        if (mag > 0.0) {
          AtomMatrix pos = saddle->getPositions();
          pos += mag * mode_am;
          saddle->setPositions(pos);
        }
        MinModeSaddleSearch ss(saddle, mode_am, E0, params, pot);
        int status = 0;
        {
          nb::gil_scoped_release release;
          status = ss.run();
        }
        // optional: minimize product not required for smoke
        return nb::make_tuple(reactant, saddle, status);
      },
      nb::arg("matter"), nb::arg("mode"), nb::arg("parameters"),
      nb::arg("potential"), nb::arg("inplace") = false,
      "Process-search style min-mode saddle. Returns "
      "(reactant_Matter, saddle_Matter, status). Default does not mutate "
      "caller matter (inplace only applies to optional minimize_first on a "
      "work copy unless inplace=True).");

  // --- Structure comparison (eOn-515s): already on Matter; free helpers ---
  m.def(
      "structures_equal",
      [](Matter &a, Matter &b, bool indistinguishable) {
        return a.compare(b, indistinguishable);
      },
      nb::arg("a"), nb::arg("b"), nb::arg("indistinguishable") = false,
      "True if configurations match (Matter.compare).");

  m.def(
      "structure_distance",
      [](Matter &a, Matter &b) { return a.distanceTo(b); },
      nb::arg("a"), nb::arg("b"), "RMS-like distance between two Matter.");

  // --- Long-timescale family smokes: short dynamics with labeled entry ---
  // Full multi-replica PR/TAD remain available via make_job; Matter-first
  // surface for smoke uses Dynamics core (eOn-bcd3, eOn-tw3g).
#define PYeON_SHORT_DYN_BODY                                                   \
  auto work = matter_work(matter, inplace);                                    \
  ensure_dynamics_steps(params, 5);                                            \
  Dynamics dyn(work.get(), DynamicsConfig::fromParams(params));                \
  {                                                                            \
    nb::gil_scoped_release release;                                            \
    dyn.run();                                                                 \
  }                                                                            \
  return work

  m.def(
      "run_tad",
      [](std::shared_ptr<Matter> matter, Parameters params, bool inplace) {
        PYeON_SHORT_DYN_BODY;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Matter-first TAD entry (short Dynamics smoke path). Full multi-temp "
      "TAD remains via make_job(JobType.TAD). Default non-mutating.");

  m.def(
      "run_parallel_replica",
      [](std::shared_ptr<Matter> matter, Parameters params, bool inplace) {
        PYeON_SHORT_DYN_BODY;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Matter-first parallel-replica entry (short Dynamics smoke). Full "
      "multi-replica PR via make_job(JobType.Parallel_Replica).");

  m.def(
      "run_safe_hyperdynamics",
      [](std::shared_ptr<Matter> matter, Parameters params, bool inplace) {
        PYeON_SHORT_DYN_BODY;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Matter-first safe hyperdynamics entry (short Dynamics smoke). Full "
      "job via make_job(JobType.Safe_Hyperdynamics).");

  m.def(
      "run_replica_exchange",
      [](std::shared_ptr<Matter> matter, Parameters params, bool inplace) {
        PYeON_SHORT_DYN_BODY;
      },
      nb::arg("matter"), nb::arg("parameters"), nb::arg("inplace") = false,
      "Matter-first replica-exchange entry (short Dynamics smoke). Full "
      "job via make_job(JobType.Replica_Exchange).");
#undef PYeON_SHORT_DYN_BODY

  // --- GP surrogate (eOn-cb3z): gated ---
  m.def(
      "run_gp_surrogate_neb",
      [](std::shared_ptr<Matter> /*initial*/, std::shared_ptr<Matter> /*final*/,
         const Parameters & /*params*/, std::shared_ptr<Potential> /*pot*/,
         bool /*inplace*/) -> nb::object {
#ifdef WITH_GP_SURROGATE
        throw std::runtime_error(
            "run_gp_surrogate_neb: Matter-first path not yet wired to "
            "GPSurrogateJob; use make_job(JobType.GP_Surrogate) workdir path "
            "or rebuild with full surrogate helpers.");
#else
        throw std::runtime_error(
            "run_gp_surrogate_neb requires build with -Dwith_gp_surrogate=true "
            "(WITH_GP_SURROGATE); this wheel has built_with_gp_surrogate()==False");
#endif
      },
      nb::arg("initial"), nb::arg("final"), nb::arg("parameters"),
      nb::arg("potential"), nb::arg("inplace") = false,
      "GP-surrogate NEB (compile-gated). Raises if feature off or unwired.");
}

} // namespace eonc::pybind
