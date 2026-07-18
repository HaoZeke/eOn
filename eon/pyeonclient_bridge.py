"""Optional bridge: eOn server Structure <-> pyeonclient.Matter."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from eon.structure import Structure


def available() -> bool:
    try:
        import pyeonclient  # noqa: F401

        return True
    except ImportError:
        return False


def to_matter(structure: "Structure", potential: Any, parameters: Any):
    from pyeonclient.bridge import structure_to_matter

    return structure_to_matter(structure, potential, parameters)


def from_matter(matter: Any) -> "Structure":
    from pyeonclient.bridge import matter_to_structure

    return matter_to_structure(matter)
