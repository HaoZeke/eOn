"""pyeonclient — nanobind surface for the eOn C++ client.

Core (``_core``)
----------------
Parameters, Potential, Matter, Job, enums, and **per-step** ClientEON helpers:
``make_job``, ``write_potcall_summary``, ``get_process_times``,
``append_results_timing``, ``steady_clock_now``.

Steps (``pyeonclient.steps``)
-----------------------------
Compose the client pipeline in Python (load → job/matter → potcalls → timing).
Use these from notebooks instead of a single opaque binary wrapper.

Bridges: :mod:`pyeonclient.bridge`, :mod:`pyeonclient.ase_bridge`.
"""

from __future__ import annotations

try:
    from pyeonclient._core import (  # type: ignore F401
        IoStatus,
        Job,
        JobType,
        Matter,
        NEBInit,
        OptType,
        Parameters,
        PbcConvention,
        Potential,
        PotType,
        RunStatus,
        __version__,
        append_results_timing,
        get_process_times,
        io_ok,
        io_status_name,
        job_type_from_name,
        job_type_name,
        make_job,
        make_potential,
        pot_type_from_name,
        pot_type_name,
        steady_clock_now,
        write_potcall_summary,
        NudgedElasticBand,
        NEBStatus,
        neb_read_file_paths,
        neb_load_path_from_files,
        neb_linear_path,
        neb_write_results,
        pot_registry_total_force_calls,
    )
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyeonclient._core is not built. Configure with "
        "-Dwith_pyeonclient=true (nanobind; stable ABI / free-threaded)."
    ) from e

from pyeonclient.bridge import (
    from_structure,
    matter_to_structure,
    structure_to_matter,
    to_structure,
)
from pyeonclient.ase_bridge import (  # noqa: E402
    ase_to_matter,
    ase_to_structure,
    conframe_to_matter,
    matter_to_ase,
    matter_to_conframe,
    structure_to_ase,
)
from pyeonclient.steps import (  # noqa: E402
    append_timing,
    load_parameters,
    minimize_workdir,
    run_eon_cwd,
    run_job,
    run_job_in_directory,
    write_minimization_results,
    neb_workdir,
)

to_ase = matter_to_ase
from_ase = ase_to_matter

__all__ = [
    "IoStatus",
    "Job",
    "JobType",
    "Matter",
    "NEBInit",
    "NEBStatus",
    "NudgedElasticBand",
    "OptType",
    "Parameters",
    "PbcConvention",
    "Potential",
    "PotType",
    "RunStatus",
    "append_results_timing",
    "append_timing",
    "ase_to_matter",
    "ase_to_structure",
    "conframe_to_matter",
    "from_ase",
    "from_structure",
    "get_process_times",
    "io_ok",
    "io_status_name",
    "job_type_from_name",
    "job_type_name",
    "load_parameters",
    "make_job",
    "make_potential",
    "matter_to_ase",
    "matter_to_conframe",
    "matter_to_structure",
    "minimize_workdir",
    "neb_linear_path",
    "neb_load_path_from_files",
    "neb_read_file_paths",
    "neb_workdir",
    "neb_write_results",
    "pot_registry_total_force_calls",
    "pot_type_from_name",
    "pot_type_name",
    "run_eon_cwd",
    "run_job",
    "run_job_in_directory",
    "steady_clock_now",
    "structure_to_ase",
    "structure_to_matter",
    "to_ase",
    "to_structure",
    "write_minimization_results",
    "write_potcall_summary",
    "__version__",
]
