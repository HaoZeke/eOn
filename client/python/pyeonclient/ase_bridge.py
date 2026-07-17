"""ASE calculators as Matter potentials (tight pyeonclient integration).

Design
------
eOn ``Matter`` owns geometry. An ASE ``Calculator`` is the force engine via
``AsePotential`` / ``attach_ase_calculator`` / ``from_ase``. MD, relax, NEB,
and dimer then call ASE **through** ``Matter.potential_energy`` / forces —
you do not drive a parallel ``ase.Atoms`` world.

Efficiency (C++ ``AseCalcPotential``)
-------------------------------------
* One peer ``ase.Atoms`` per calculator (not rebuilt each force call).
* ``bind_matter(m)`` wires shared positions when ASE allows and uses
  ``Matter.geometry_generation`` for ASE ``system_changes`` (only
  ``positions`` / ``cell`` when possible, not full ``all_changes``).
* Foreign buffers (``Potential.get_ef``) still use bulk NumPy views.

Preferred usage::

    import pyeonclient as pyec
    from ase.calculators.emt import EMT

    # 1) From ASE structure + calc (binds automatically)
    atoms.calc = EMT()
    m = pyec.from_ase(atoms)
    m.relax()                       # ASE forces under the hood

    # 2) Explicit Matter-first
    handle = pyec.ase_potential(EMT())
    m = pyec.Matter(handle.potential, pyec.Parameters())
    m.resize(n)
    m.positions = ...
    handle.bind_matter(m)
    e = m.potential_energy

    # 3) Attach calc to existing Matter
    pyec.attach_ase_calculator(m, EMT())
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


def ase_potential(calculator: Any) -> Any:
    """Wrap an ASE Calculator as a first-class pyeonclient potential handle.

    Returns an ``AsePotential`` with ``.potential`` (for ``Matter``) and
    ``.bind_matter(m)`` for shared geometry + system_changes caching.
    """
    if calculator is None:
        raise ValueError("ase_potential: calculator is None")
    from pyeonclient._core import ase_potential as _ase_potential

    return _ase_potential(calculator)


def potential_from_ase(calculator: Any, parameters: Any = None) -> Any:
    """Wrap a live ASE Calculator as an eOn :class:`Potential`.

    Prefer :func:`ase_potential` + ``bind_matter`` for Matter workflows.
    Does **not** require ``-Dwith_ase`` at compile time.
    """
    if calculator is None:
        raise ValueError("potential_from_ase: calculator is None")
    from pyeonclient._core import make_potential_from_ase

    return make_potential_from_ase(calculator)


def attach_ase_calculator(matter: Any, calculator: Any) -> None:
    """Set ``matter``'s potential to ``calculator`` and bind shared geometry."""
    if calculator is None:
        raise ValueError("attach_ase_calculator: calculator is None")
    from pyeonclient._core import attach_ase_calculator as _attach

    _attach(matter, calculator)


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
    """Build :class:`pyeonclient.Matter` from ``ase.Atoms``.

    If ``potential`` is omitted, uses ``atoms.calc``. After construction the
    ASE pot is **bound** to the Matter so subsequent force calls use shared
    geometry + system_changes caching.
    """
    from pyeonclient import Parameters
    from pyeonclient._core import bind_ase_matter, matter_from_ase

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

    matter = matter_from_ase(atoms, potential, parameters)
    bind_ase_matter(potential, matter)
    return matter


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
