"""pyeonclient — complete nanobind surface for the eOn C++ client.

Core (``_core`` extension)
--------------------------
Parameters, Potential, Matter, Job, make_job, run_job_in_directory, enums.

Bridges
-------
* :mod:`pyeonclient.bridge` — Structure ↔ Matter (no ASE)
* :mod:`pyeonclient.ase_bridge` — ASE ↔ Matter / Structure (optional ASE)

Build: ``-Dwith_pyeonclient=true`` (nanobind only; abi3 / free-threaded).
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
        io_ok,
        io_status_name,
        job_type_from_name,
        job_type_name,
        make_job,
        make_potential,
        pot_type_from_name,
        pot_type_name,
        run_job_in_directory,
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

# ASE helpers: importable names; functions raise if ASE missing
from pyeonclient.ase_bridge import (  # noqa: E402
    ase_to_matter,
    ase_to_structure,
    conframe_to_matter,
    matter_to_ase,
    matter_to_conframe,
    structure_to_ase,
)

# Friendly aliases
to_ase = matter_to_ase
from_ase = ase_to_matter

__all__ = [
    "IoStatus",
    "Job",
    "JobType",
    "Matter",
    "NEBInit",
    "OptType",
    "Parameters",
    "PbcConvention",
    "Potential",
    "PotType",
    "RunStatus",
    "ase_to_matter",
    "ase_to_structure",
    "conframe_to_matter",
    "from_ase",
    "from_structure",
    "io_ok",
    "io_status_name",
    "job_type_from_name",
    "job_type_name",
    "make_job",
    "make_potential",
    "matter_to_ase",
    "matter_to_conframe",
    "matter_to_structure",
    "pot_type_from_name",
    "pot_type_name",
    "run_job_in_directory",
    "structure_to_ase",
    "structure_to_matter",
    "to_ase",
    "to_structure",
    "__version__",
]
