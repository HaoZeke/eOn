"""Enums shared by API specs (and eventually config models).

String values match pyeonclient / C++ routing (not necessarily INI keys).
"""

from __future__ import annotations

from enum import Enum


class MinModeMethod(str, Enum):
    improved = "improved"
    classic = "classic"
    lanczos = "lanczos"
    davidson = "davidson"


class Accelerant(str, Enum):
    """Optional accelerant on dimer/NEB stacks.

    ``gp`` is Gaussian-process only (not Prefactor, not Parallel Replica).
    For dimer, ``gp`` requires :attr:`MinModeMethod.improved`.
    """

    none = "none"
    gp = "gp"


class PathInit(str, Enum):
    """Maps to client ``NEBInit``."""

    linear = "linear"
    idpp = "idpp"
    idpp_collective = "idpp_collective"
    sidpp = "sidpp"
    sidpp_zbl = "sidpp_zbl"
    file = "file"  # path list via Parameters.neb_initial_path
