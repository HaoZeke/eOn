"""Seamless ASE Calculator → Matter Potential."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")
pytest.importorskip("ase")
from ase import Atoms
from ase.calculators.lj import LennardJones


def test_potential_from_ase_lj_energy_forces():
    calc = LennardJones(epsilon=1.0, sigma=1.0, rc=10.0, ro=0.0, smooth=False)
    pot = pyec.potential_from_ase(calc)
    assert pot is not None
    # two atoms
    pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
    z = np.array([13, 13], dtype=np.int64)  # Al Z
    cell = np.eye(3) * 10.0
    E, F = pot.get_ef(pos, z, cell)
    assert np.isfinite(E)
    assert F.shape == (2, 3)
    # opposite forces along x
    assert F[0, 0] > 0 and F[1, 0] < 0


def test_from_ase_uses_atoms_calc():
    atoms = Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        cell=np.eye(3) * 12.0,
        pbc=True,
    )
    atoms.calc = LennardJones(epsilon=1.0, sigma=1.0, rc=5.0, smooth=False)
    # seamless: no explicit pot/params
    m = pyec.from_ase(atoms)
    assert m.n_atoms == 2
    e = float(m.potential_energy)
    assert np.isfinite(e)
    # force path via relax step or max_force
    f = float(m.max_force)
    assert np.isfinite(f)


def test_from_ase_without_calc_requires_pot():
    atoms = Atoms("H", positions=[[0, 0, 0]], cell=np.eye(3) * 5, pbc=False)
    with pytest.raises(ValueError, match="atoms.calc|potential"):
        pyec.from_ase(atoms)
