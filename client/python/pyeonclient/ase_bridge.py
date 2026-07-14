"""ASE ↔ Matter / Structure / readcon helpers (optional ASE dependency).

Prefer these over ad-hoc conversion. When ASE is not installed, functions
raise ``ImportError`` with a clear message (server path does not require ASE).
"""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms as AseAtoms


def _require_ase():
    try:
        import ase
        from ase import Atoms
        from ase.constraints import FixAtoms
    except ImportError as e:
        raise ImportError(
            "ASE is required for to_ase/from_ase helpers "
            "(optional: pip/pixi install ase)"
        ) from e
    return ase, Atoms, FixAtoms


def matter_to_ase(matter: Any, *, pbc: Optional[bool] = None) -> "AseAtoms":
    """Convert a live :class:`pyeonclient.Matter` to ``ase.Atoms``."""
    _, Atoms, FixAtoms = _require_ase()

    pos = np.asarray(matter.positions, dtype=float)
    cell = np.asarray(matter.cell, dtype=float)
    numbers = np.asarray(matter.atomic_numbers, dtype=int)
    masses = np.asarray(matter.masses, dtype=float)
    fixed = np.asarray(matter.fixed, dtype=int).reshape(-1)

    use_pbc = bool(matter.periodic) if pbc is None else bool(pbc)
    atoms = Atoms(
        numbers=numbers,
        positions=pos,
        cell=cell,
        pbc=use_pbc,
        masses=masses,
    )
    mask = fixed != 0
    if np.any(mask):
        atoms.set_constraint(FixAtoms(mask=mask.tolist()))
    return atoms


def ase_to_matter(
    atoms: "AseAtoms",
    potential: Any,
    parameters: Any,
) -> Any:
    """Build :class:`pyeonclient.Matter` from ``ase.Atoms``."""
    from pyeonclient import Matter

    _, _, FixAtoms = _require_ase()

    n = len(atoms)
    m = Matter(potential, parameters)
    m.resize(n)
    m.positions = np.ascontiguousarray(atoms.get_positions(), dtype=np.float64)
    m.cell = np.ascontiguousarray(np.asarray(atoms.get_cell()), dtype=np.float64)
    m.masses = np.ascontiguousarray(atoms.get_masses(), dtype=np.float64)
    m.atomic_numbers = np.ascontiguousarray(
        atoms.get_atomic_numbers(), dtype=np.int64
    )

    fixed = np.zeros(n, dtype=np.int64)
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            idx = np.asarray(c.get_indices(), dtype=int).reshape(-1)
            fixed[idx] = 1
    m.fixed = fixed
    m.periodic = bool(np.any(atoms.pbc))
    return m


def structure_to_ase(structure: Any, *, pbc: bool = True) -> "AseAtoms":
    """Convert eOn server :class:`~eon.structure.Structure` to ``ase.Atoms``."""
    _, Atoms, FixAtoms = _require_ase()

    from pyeonclient.bridge import _symbol_to_z

    n = len(structure)
    numbers = np.array(
        [_symbol_to_z(name) for name in structure.names], dtype=int
    )
    free = np.asarray(structure.free, dtype=float).reshape(-1)
    atoms = Atoms(
        numbers=numbers,
        positions=np.asarray(structure.r, dtype=float),
        cell=np.asarray(structure.box, dtype=float),
        pbc=pbc,
        masses=np.asarray(structure.mass, dtype=float),
    )
    mask = free < 0.5
    if np.any(mask):
        atoms.set_constraint(FixAtoms(mask=mask.tolist()))
    return atoms


def ase_to_structure(atoms: "AseAtoms") -> Any:
    """Convert ``ase.Atoms`` to eOn server :class:`~eon.structure.Structure`."""
    from eon.structure import Structure
    from pyeonclient.bridge import _z_to_symbol

    _, _, FixAtoms = _require_ase()

    n = len(atoms)
    s = Structure(n)
    s.r = np.asarray(atoms.get_positions(), dtype=float).copy()
    s.box = np.asarray(atoms.get_cell(), dtype=float).copy()
    s.mass = np.asarray(atoms.get_masses(), dtype=float).copy()
    z = np.asarray(atoms.get_atomic_numbers(), dtype=int)
    s.names = [_z_to_symbol(int(zi)) for zi in z]

    free = np.ones(n, dtype=float)
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            free[np.asarray(c.get_indices(), dtype=int)] = 0.0
    s.free = free
    return s


def matter_to_conframe(matter: Any):
    """Matter → readcon.ConFrame (via temporary .con round-trip or arrays).

    Uses in-memory write_con_string path through Structure for fidelity.
    """
    import readcon
    from pyeonclient.bridge import matter_to_structure

    s = matter_to_structure(matter)
    return s.to_conframe()


def conframe_to_matter(frame, potential: Any, parameters: Any):
    """readcon.ConFrame → Matter."""
    from eon.structure import Structure
    from pyeonclient.bridge import structure_to_matter

    s = Structure.from_conframe(frame)
    return structure_to_matter(s, potential, parameters)
