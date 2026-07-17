"""pyeonclient — in-process eOn client (nanobind).

First-class algorithm surface (ASE-shaped geometries as :class:`Matter`)::

    # Minimization
    matter = from_ase(atoms, pot, params); matter.relax()

    # Min-mode
    dimer = ImprovedDimer(matter, params, pot)
    dimer.compute(matter, direction)          # direction (n, 3)
    ev, mode = dimer.eigenvalue, dimer.eigenvector

    # Single-ended saddle
    ss = MinModeSaddleSearch(matter, mode, E_react, params, pot)
    status = ss.run()                         # mutates matter → saddle

    # NEB band
    path = [from_ase(img, pot, params) for img in images]
    neb = NudgedElasticBand(path, params, pot); neb.compute()

    # Hessian / HTST prefactors
    H = Hessian(params, matter)
    freqs = H.get_freqs(matter, all_free_atoms(matter))
    pref1, pref2 = get_prefactors(params, min1, saddle, min2)

Modules: ``_core`` (compiled), :mod:`pyeonclient.ase_bridge`,
:mod:`pyeonclient.bridge`, :mod:`pyeonclient.steps` (optional workdir helpers).

Docs: https://eondocs.org/user_guide/pyeonclient.html
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
        # NEB
        NudgedElasticBand,
        NEBStatus,
        neb_read_file_paths,
        neb_load_path_from_files,
        neb_linear_path,
        neb_write_results,
        pot_registry_total_force_calls,
        # Min-mode
        Dimer,
        ImprovedDimer,
        Lanczos,
        Davidson,
        # Saddle
        MinModeSaddleSearch,
        SaddleStatus,
        saddle_status_message,
        # Analysis
        Hessian,
        get_prefactors,
        moved_atoms,
        all_free_atoms,
        # Build probes
        built_with_metatomic,
        built_with_rgpot,
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
    rgpot_metatomic_workdir,
    rgpot_metatomic_neb_workdir,
)

to_ase = matter_to_ase
from_ase = ase_to_matter

__all__ = [
    # Core types
    "IoStatus",
    "Job",
    "JobType",
    "Matter",
    "Parameters",
    "PbcConvention",
    "Potential",
    "PotType",
    "OptType",
    "RunStatus",
    # NEB
    "NEBInit",
    "NEBStatus",
    "NudgedElasticBand",
    "neb_linear_path",
    "neb_load_path_from_files",
    "neb_read_file_paths",
    "neb_write_results",
    # Min-mode
    "Dimer",
    "ImprovedDimer",
    "Lanczos",
    "Davidson",
    # Saddle
    "MinModeSaddleSearch",
    "SaddleStatus",
    "saddle_status_message",
    # Analysis
    "Hessian",
    "get_prefactors",
    "moved_atoms",
    "all_free_atoms",
    # Factories / steps
    "make_job",
    "make_potential",
    "load_parameters",
    "run_job",
    "run_job_in_directory",
    "run_eon_cwd",
    "minimize_workdir",
    "neb_workdir",
    "rgpot_metatomic_workdir",
    "rgpot_metatomic_neb_workdir",
    "write_minimization_results",
    "write_potcall_summary",
    "append_results_timing",
    "append_timing",
    "get_process_times",
    "steady_clock_now",
    "pot_registry_total_force_calls",
    # Enum helpers
    "io_ok",
    "io_status_name",
    "job_type_from_name",
    "job_type_name",
    "pot_type_from_name",
    "pot_type_name",
    # Bridges
    "from_ase",
    "to_ase",
    "ase_to_matter",
    "matter_to_ase",
    "ase_to_structure",
    "structure_to_ase",
    "conframe_to_matter",
    "matter_to_conframe",
    "from_structure",
    "to_structure",
    "structure_to_matter",
    "matter_to_structure",
    # Probes
    "built_with_metatomic",
    "built_with_rgpot",
    "__version__",
]
