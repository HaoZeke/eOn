"""ClientEON pipeline as explicit Python steps.

Each function is one stage of the standalone ``eonclient`` main loop. Compose
them in notebooks/drivers instead of a single opaque entry point.

Typical job path (any JobType including NEB)::

    params = load_parameters("config.ini")
    t0 = steady_clock_now()
    files = run_job(params)          # make_job + Job.run + drop Job
    write_potcall_summary()
    append_timing("results.dat", t0)

Matter path (minimization without Job wrapper)::

    params = load_parameters("config.ini")
    pot = make_potential(params.potential, params)
    m = Matter(pot, params)
    m.con2matter("pos.con")
    ok = m.relax(write_movie=True, prefix_movie="minimization",
                 prefix_checkpoint="pos")
    m.matter2con("min.con")
    write_minimization_results(params, pot, m, converged=ok)
    del m, pot
    write_potcall_summary()
    append_timing("results.dat", t0)
"""

from __future__ import annotations

import os
from contextlib import chdir
from pathlib import Path
from typing import Any, Sequence

from pyeonclient._core import (
    JobType,
    Matter,
    Parameters,
    Potential,
    RunStatus,
    append_results_timing,
    get_process_times,
    make_job,
    make_potential,
    steady_clock_now,
    write_potcall_summary,
)


def load_parameters(path: str | Path = "config.ini") -> Parameters:
    """Step: load ``config.ini`` (or JSON via Parameters.load_json)."""
    p = Parameters()
    p.load(str(path))
    return p


def run_job(params: Parameters) -> list[str]:
    """Step: ``make_job`` → ``Job.run`` → destroy Job (Potential released).

    Does not write potcall summary or timing; call those next.
    """
    job = make_job(params)
    try:
        return list(job.run())
    finally:
        # Explicit drop so PotRegistry sees destruction before summary.
        del job


def append_timing(
    results_path: str | Path = "results.dat",
    t0: float | None = None,
    *,
    elapsed: float | None = None,
) -> None:
    """Step: append ``time_seconds`` / ``user_time`` / ``system_time``.

    Pass *t0* from :func:`steady_clock_now` before the work, or *elapsed*
    wall seconds directly.
    """
    if elapsed is None:
        if t0 is None:
            raise ValueError("append_timing requires t0= or elapsed=")
        elapsed = steady_clock_now() - t0
    _real, user, system = get_process_times()
    append_results_timing(str(results_path), float(elapsed), float(user), float(system))


def write_minimization_results(
    params: Parameters,
    pot: Potential,
    matter: Matter,
    *,
    converged: bool,
    path: str | Path = "results.dat",
) -> None:
    """Step: write MinimizationJob-style results.dat body (no timing footer)."""
    status = RunStatus.GOOD if converged else RunStatus.FAIL_MAX_ITERATIONS
    # magic_enum-style names (same strings as MinimizationJob)
    pot_name = str(params.potential).rsplit(".", 1)[-1]
    status_name = str(status).rsplit(".", 1)[-1]
    status_i = 0 if status == RunStatus.GOOD else 1
    lines = [
        f"{status_i} termination_reason",
        f"{status_name} termination_reason_text",
        "minimization job_type",
        f"{pot_name} potential_type",
        f"{int(pot.force_call_counter)} total_force_calls",
    ]
    if status != RunStatus.FAIL_POTENTIAL_FAILED:
        lines.append(f"{matter.potential_energy} potential_energy")
    Path(path).write_text("\n".join(lines) + "\n")


def minimize_workdir(
    workdir: str | Path,
    *,
    pos: str = "pos.con",
    min_con: str = "min.con",
    movie_prefix: str = "minimization",
    write_movie: bool = True,
) -> list[str]:
    """Compose Matter steps for a minimization workdir (cookbook min_reactant).

    Steps: load_parameters → make_potential → Matter → con2matter → relax →
    matter2con → write_minimization_results → del handles → potcalls → timing.
    """
    work = Path(workdir).resolve()
    files: list[str] = []
    with chdir(work):
        t0 = steady_clock_now()
        params = load_parameters("config.ini")
        pot = make_potential(params.potential, params)
        matter = Matter(pot, params)
        st = matter.con2matter(pos)
        # IoStatus may be enum; accept truthy / Ok
        if hasattr(st, "name") and st.name not in ("Ok", "OK"):
            from pyeonclient._core import io_ok

            if not io_ok(st):
                raise RuntimeError(f"con2matter failed for {pos}: {st}")
        quiet = bool(getattr(params, "quiet", False))
        ckpt = bool(getattr(params, "checkpoint", False))
        converged = bool(
            matter.relax(
                quiet=quiet,
                write_movie=write_movie,
                checkpoint=ckpt,
                prefix_movie=movie_prefix,
                prefix_checkpoint="pos",
            )
        )
        matter.matter2con(min_con)
        files.append(min_con)
        write_minimization_results(params, pot, matter, converged=converged)
        files.append("results.dat")
        del matter
        del pot
        write_potcall_summary("_potcalls.json")
        files.append("_potcalls.json")
        append_timing("results.dat", t0)
    return files


def run_job_in_directory(
    workdir: str | Path,
    params: Parameters | None = None,
) -> list[str]:
    """Compose ClientEON job steps in *workdir* (any job type, including NEB).

    Steps: (optional load config.ini) → run_job → write_potcall_summary →
    append_timing. Prefer :func:`minimize_workdir` when you want Matter.relax
    control rather than MinimizationJob.
    """
    work = Path(workdir).resolve()
    with chdir(work):
        t0 = steady_clock_now()
        if params is None:
            params = Parameters()
        ini = work / "config.ini"
        if ini.is_file():
            params.load(str(ini))
        files = run_job(params)
        write_potcall_summary("_potcalls.json")
        if "_potcalls.json" not in files:
            files.append("_potcalls.json")
        append_timing("results.dat", t0)
        return files


def run_eon_cwd() -> list[str]:
    """Cookbook helper: run whatever job config.ini requests in cwd."""
    return run_job_in_directory(Path("."), Parameters())
