/*
** Matter-first sampling / long-timescale primitives with explicit inplace.
** Issues: eOn-p8em dynamics, eOn-n6oi MC, eOn-fmxm BH, eOn-412j process search,
** eOn-du2q TAD, eOn-v1d0 PR/SafeHyper/REX, eOn-nmcu structure comparison,
** eOn-hyjn GP surrogate (gated), eOn-rczn inplace contract.
*/
#include "Dynamics.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "MonteCarlo.h"
#include "Parameters.h"
#include "ParallelReplicaJob.h"
#include "Potential.h"
#include "ReplicaExchangeJob.h"
#include "SafeHyperJob.h"
#include "TADJob.h"
#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"

#ifdef WITH_GP_SURROGATE
#include "GPSurrogateJob.h"
#include "NudgedElasticBand.h"
#endif

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
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

/// Ensure PR/TAD/REX interval fields are finite and compatible with a short
/// Matter-first smoke run (avoids zero-length buffers / div-by-zero).
void ensure_long_timescale_params(eonc::Parameters &params,
                                  long default_steps) {
  ensure_dynamics_steps(params, default_steps);
  auto &dyn = params.dynamics_options;
  const double dt = dyn.time_step;
  const double horizon = dt * static_cast<double>(dyn.steps);

  auto &pr = params.parallel_replica_options;
  if (pr.state_check_interval <= 0.0)
    pr.state_check_interval = std::max(dt, horizon);
  if (pr.record_interval <= 0.0)
    pr.record_interval = std::max(dt, pr.state_check_interval / 2.0);
  if (pr.record_interval > pr.state_check_interval)
    pr.record_interval = pr.state_check_interval;
  // Ensure at least one MD step per interval (integer division in TAD/PR).
  if (pr.state_check_interval < dt)
    pr.state_check_interval = dt;
  if (pr.record_interval < dt)
    pr.record_interval = dt;
  if (pr.dephase_time <= 0.0)
    pr.dephase_time = dt; // one step dephase for smoke
  if (pr.corr_time <= 0.0)
    pr.corr_time = dt;
  // Short Matter-first smokes: skip refine (needs filled snapshot buffers)
  // and always cap dephase so LJ dimers cannot hang forever.
  if (dyn.steps > 0 && dyn.steps < 50) {
    pr.refine_transition = false;
    pr.dephase_loop_stop = true;
    if (pr.dephase_loop_max <= 0 || pr.dephase_loop_max > 5)
      pr.dephase_loop_max = 2;
  }

  auto &rex = params.replica_exchange_options;
  if (rex.replicas < 2)
    rex.replicas = 2;
  // Prefer explicit dynamics budget over stale sampling_time_input conversion.
  if (rex.sampling_time <= 0.0 ||
      (dyn.steps > 0 && dyn.steps < 50 &&
       rex.sampling_time > horizon * 2.0)) {
    rex.sampling_time = horizon > 0.0 ? horizon : dt;
  }
  if (rex.exchange_period <= 0.0 ||
      rex.exchange_period > rex.sampling_time)
    rex.exchange_period = std::max(dt, rex.sampling_time / 2.0);
  if (rex.temperature_low <= 0.0)
    rex.temperature_low = params.main_options.temperature > 0.0
                              ? params.main_options.temperature
                              : 300.0;
  if (rex.temperature_high <= rex.temperature_low)
    rex.temperature_high = rex.temperature_low * 1.5;
  if (rex.exchange_trials <= 0)
    rex.exchange_trials = std::max(1L, rex.replicas - 1);
  // Short Matter-first smoke: two replicas is enough to exercise Metropolis.
  if (dyn.steps > 0 && dyn.steps < 50 && rex.replicas > 3)
    rex.replicas = 2;
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

  // --- Long-timescale jobs via real Job classes (eOn-du2q, eOn-v1d0) ---
  // Never Dynamics aliases. Default non-mutating; inplace=True writes back.
  m.def(
      "run_tad",
      [](std::shared_ptr<Matter> matter, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        ensure_long_timescale_params(params, 5);
        auto work = matter_work(matter, inplace);
        if (pot)
          work->setPotential(pot);
        auto job_pot = pot ? pot : work->getPotential();
        TADJob job(job_pot, params);
        std::shared_ptr<Matter> result;
        {
          nb::gil_scoped_release release;
          result = job.runFromMatter(work);
        }
        return matter_result(matter, result, inplace);
      },
      nb::arg("matter"), nb::arg("parameters"),
      nb::arg("potential") = nb::none(), nb::arg("inplace") = false,
      "Matter-first Temperature Accelerated Dynamics (TADJob). "
      "Default non-mutating copy; inplace=True updates matter.");

  m.def(
      "run_parallel_replica",
      [](std::shared_ptr<Matter> matter, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        ensure_long_timescale_params(params, 5);
        auto work = matter_work(matter, inplace);
        if (pot)
          work->setPotential(pot);
        auto job_pot = pot ? pot : work->getPotential();
        ParallelReplicaJob job(job_pot, params);
        std::shared_ptr<Matter> result;
        {
          nb::gil_scoped_release release;
          result = job.runFromMatter(work);
        }
        return matter_result(matter, result, inplace);
      },
      nb::arg("matter"), nb::arg("parameters"),
      nb::arg("potential") = nb::none(), nb::arg("inplace") = false,
      "Matter-first Parallel Replica Dynamics (ParallelReplicaJob). "
      "Default non-mutating; inplace=True updates matter with final trajectory.");

  m.def(
      "run_safe_hyperdynamics",
      [](std::shared_ptr<Matter> matter, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        ensure_long_timescale_params(params, 5);
        auto work = matter_work(matter, inplace);
        if (pot)
          work->setPotential(pot);
        auto job_pot = pot ? pot : work->getPotential();
        SafeHyperJob job(job_pot, params);
        std::shared_ptr<Matter> result;
        {
          nb::gil_scoped_release release;
          result = job.runFromMatter(work);
        }
        return matter_result(matter, result, inplace);
      },
      nb::arg("matter"), nb::arg("parameters"),
      nb::arg("potential") = nb::none(), nb::arg("inplace") = false,
      "Matter-first Safe Hyperdynamics (SafeHyperJob). "
      "Default non-mutating; inplace=True updates matter.");

  m.def(
      "run_replica_exchange",
      [](std::shared_ptr<Matter> matter, Parameters params,
         std::shared_ptr<Potential> pot, bool inplace) {
        ensure_long_timescale_params(params, 5);
        auto work = matter_work(matter, inplace);
        if (pot)
          work->setPotential(pot);
        auto job_pot = pot ? pot : work->getPotential();
        ReplicaExchangeJob job(job_pot, params);
        std::shared_ptr<Matter> result;
        {
          nb::gil_scoped_release release;
          result = job.runFromMatter(work);
        }
        return matter_result(matter, result, inplace);
      },
      nb::arg("matter"), nb::arg("parameters"),
      nb::arg("potential") = nb::none(), nb::arg("inplace") = false,
      "Matter-first Replica Exchange (ReplicaExchangeJob). "
      "Default non-mutating; returns replica-0 Matter after sampling.");

  // --- GP surrogate (eOn-hyjn): gated ---
  m.def(
      "run_gp_surrogate_neb",
      [](std::shared_ptr<Matter> initial, std::shared_ptr<Matter> final_state,
         Parameters params, std::shared_ptr<Potential> pot,
         bool inplace) -> nb::object {
#ifdef WITH_GP_SURROGATE
        if (!initial || !final_state)
          throw std::runtime_error(
              "run_gp_surrogate_neb: initial and final Matter required");
        auto work_i = matter_work(initial, inplace);
        // final is never mutated by the surrogate path; always copy
        auto work_f = std::make_shared<Matter>(*final_state);
        if (pot) {
          work_i->setPotential(pot);
          work_f->setPotential(pot);
        }
        auto job_pot = pot ? pot : work_i->getPotential();
        GPSurrogateJob job(job_pot, params);
        std::shared_ptr<NudgedElasticBand> neb;
        {
          nb::gil_scoped_release release;
          neb = job.runFromMatter(work_i, work_f);
        }
        return nb::cast(neb);
#else
        (void)initial;
        (void)final_state;
        (void)params;
        (void)pot;
        (void)inplace;
        throw std::runtime_error(
            "run_gp_surrogate_neb requires build with -Dwith_gp_surrogate=true "
            "(WITH_GP_SURROGATE); this wheel has built_with_gp_surrogate()==False");
#endif
      },
      nb::arg("initial"), nb::arg("final"), nb::arg("parameters"),
      nb::arg("potential"), nb::arg("inplace") = false,
      "GP-surrogate NEB via GPSurrogateJob (compile-gated WITH_GP_SURROGATE). "
      "Returns NudgedElasticBand. Raises RuntimeError if feature off.");
}

} // namespace eonc::pybind
