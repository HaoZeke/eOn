"""ASE ↔ Matter / Structure / readcon helpers (optional ASE dependency).

Hot path is C++ nanobind (``matter_from_ase`` / ``matter_to_ase`` /
``make_potential_from_ase``): bulk ndarray buffers, zero-copy Matter views.

Seamless calculator path::

    import pyeonclient as pyec
    from ase.calculators.emt import EMT

    atoms.calc = EMT()
    matter = pyec.from_ase(atoms)              # C++ bulk fill + ASE pot wrap
    # or
    pot = pyec.potential_from_ase(atoms.calc)
    matter = pyec.from_ase(atoms, pot, params)

    # Zero-copy geometry views into Matter (mutate + mark_geometry_dirty)
    matter.positions[:] = atoms.get_positions()
    matter.mark_geometry_dirty()
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


def potential_from_ase(calculator: Any, parameters: Any = None) -> Any:
    """Wrap a live ASE Calculator as an eOn :class:`Potential`.

    Does **not** require a file-based ASE_POT script or ``-Dwith_ase``.
    Holds a reference to ``calculator``; not thread-safe across images.
    """
    if calculator is None:
        raise ValueError("potential_from_ase: calculator is None")
    from pyeonclient._core import make_potential_from_ase

    return make_potential_from_ase(calculator)


def matter_to_ase(matter: Any, *, pbc: Optional[bool] = None) -> "AseAtoms":
    """Convert a live :class:`pyeonclient.Matter` to ``ase.Atoms`` (C++ path)."""
    _require_ase()
    from pyeonclient._core import matter_to_ase as _matter_to_ase

    return _matter_to_ase(matter, pbc)


def ase_to_matter(
    atoms: "AseAtoms",
    potential: Any = None,
    parameters: Any = None,
) -> Any:
    """Build :class:`pyeonclient.Matter` from ``ase.Atoms`` (C++ bulk path).

    If ``potential`` is omitted, uses ``atoms.calc`` via
    :func:`potential_from_ase`. If ``parameters`` is omitted, a default
    :class:`Parameters` is created.
    """
    from pyeonclient import Parameters
    from pyeonclient._core import matter_from_ase

    _require_ase()

    if parameters is None:
        parameters = Parameters()
    if potential is None:
        calc = getattr(atoms, "calc", None)
        if calc is None:
            raise ValueError(
                "from_ase/ase_to_matter: pass potential= or attach atoms.calc "
                "(ASE Calculator) so eOn can evaluate forces"
            )
        potential = potential_from_ase(calc, parameters)

    return matter_from_ase(atoms, potential, parameters)


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
    """Matter → readcon.ConFrame (via Structure arrays)."""
    from pyeonclient.bridge import matter_to_structure

    s = matter_to_structure(matter)
    return s.to_conframe()


def conframe_to_matter(frame, potential: Any, parameters: Any):
    """readcon.ConFrame → Matter."""
    from eon.structure import Structure
    from pyeonclient.bridge import structure_to_matter

    s = Structure.from_conframe(frame)
    return structure_to_matter(s, potential, parameters)
