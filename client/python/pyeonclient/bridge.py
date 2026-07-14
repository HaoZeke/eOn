"""Structure (server) <-> Matter (client) without ASE."""

from __future__ import annotations

from typing import Any

import numpy as np

# Minimal Z table for common eOn fixtures; extended via elements int keys when available.
_SYMBOL_Z = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Ti": 22,
    "Fe": 26,
    "Cu": 29,
    "Zn": 30,
    "Pt": 78,
    "Au": 79,
}


def _symbol_to_z(sym: str) -> int:
    s = (sym or "").strip()
    if not s:
        raise KeyError("empty chemical symbol")
    if s in _SYMBOL_Z:
        return _SYMBOL_Z[s]
    t = s[:1].upper() + s[1:].lower()
    if t in _SYMBOL_Z:
        return _SYMBOL_Z[t]
    # Fall back to eon.atoms.elements int keys (symbol -> entry -> reverse)
    try:
        from eon.atoms import elements

        info = elements.get(s) or elements.get(t)
        if info is not None:
            # elements is also keyed by Z; find matching symbol
            for k, v in elements.items():
                if isinstance(k, int) and v is info:
                    return int(k)
    except Exception:
        pass
    raise KeyError(
        f"no atomic number for symbol {sym!r}; extend pyeonclient.bridge._SYMBOL_Z"
    )


def _z_to_symbol(z: int) -> str:
    try:
        from eon.atoms import elements

        info = elements.get(int(z))
        if info is not None:
            return str(info["symbol"])
    except Exception:
        pass
    for sym, zi in _SYMBOL_Z.items():
        if zi == int(z):
            return sym
    return f"Z{int(z)}"


def structure_to_matter(structure: Any, potential: Any, parameters: Any) -> Any:
    """Build a live :class:`Matter` from a server Structure/Atoms-like object."""
    from pyeonclient import Matter

    n = len(structure)
    m = Matter(potential, parameters)
    m.resize(n)
    m.positions = np.ascontiguousarray(structure.r, dtype=np.float64)
    m.cell = np.ascontiguousarray(structure.box, dtype=np.float64)
    m.masses = np.ascontiguousarray(structure.mass, dtype=np.float64)
    free = np.asarray(structure.free, dtype=np.float64).reshape(-1)
    m.fixed = (free < 0.5).astype(np.int64)
    z = np.empty(n, dtype=np.int64)
    for i, name in enumerate(structure.names):
        z[i] = _symbol_to_z(name)
    m.atomic_numbers = z
    m.periodic = True
    return m


def matter_to_structure(matter: Any) -> Any:
    """Convert a live :class:`Matter` to :class:`eon.structure.Structure`."""
    from eon.structure import Structure

    n = int(matter.n_atoms)
    s = Structure(n)
    s.r = np.asarray(matter.positions, dtype=float).copy()
    s.box = np.asarray(matter.cell, dtype=float).copy()
    s.mass = np.asarray(matter.masses, dtype=float).copy()
    fixed = np.asarray(matter.fixed, dtype=np.int64).reshape(-1)
    s.free = np.where(fixed != 0, 0.0, 1.0)
    z = np.asarray(matter.atomic_numbers, dtype=np.int64).reshape(-1)
    s.names = [_z_to_symbol(int(z[i])) for i in range(n)]
    return s
