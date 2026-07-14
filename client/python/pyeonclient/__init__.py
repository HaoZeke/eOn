"""pyeonclient — nanobind bindings to the eOn C++ client core.

Public surface (from ``_core`` extension)
-----------------------------------------
* ``Parameters`` — config.ini / JSON parameter bag
* ``Potential`` / ``make_potential`` — forcefield factory
* ``Matter`` — structure + energy/forces + ``con2matter`` / ``matter2con``
* ``IoStatus``, ``PotType``, ``JobType``, ``PbcConvention``, ``RunStatus``

Server bridge
-------------
* ``structure_to_matter`` / ``matter_to_structure`` in :mod:`pyeonclient.bridge`
  convert eOn server :class:`~eon.structure.Structure` objects without ASE.

Build with ``-Dwith_pyeonclient=true`` (nanobind only; stable ABI abi3 on CPython ≥ 3.12, free-threaded when ``Py_GIL_DISABLED``). No pybind11.
"""

from __future__ import annotations

try:
    from pyeonclient._core import (  # type: ignore F401
        IoStatus,
        JobType,
        Matter,
        Parameters,
        PbcConvention,
        Potential,
        PotType,
        RunStatus,
        __version__,
        io_ok,
        io_status_name,
        make_potential,
        pot_type_from_name,
        pot_type_name,
    )
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyeonclient._core is not built. Configure with "
        "-Dwith_pyeonclient=true (and install nanobind), then rebuild."
    ) from e

from pyeonclient.bridge import matter_to_structure, structure_to_matter

__all__ = [
    "IoStatus",
    "JobType",
    "Matter",
    "Parameters",
    "PbcConvention",
    "Potential",
    "PotType",
    "RunStatus",
    "io_ok",
    "io_status_name",
    "make_potential",
    "matter_to_structure",
    "pot_type_from_name",
    "pot_type_name",
    "structure_to_matter",
    "__version__",
]
