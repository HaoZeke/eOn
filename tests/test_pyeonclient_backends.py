"""Key→factory backends registry + ASE pot smoke."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")
from pyeonclient.backends import list_backends, make_backend


def test_list_backends_has_keys():
    keys = list_backends()
    assert "ase" in keys
    assert "lj" in keys
    assert "metatomic" in keys
    assert "rgpot_metatomic" in keys
    assert "ase_metatomic" in keys


def test_make_backend_lj_on_matter():
    pot = make_backend("lj")
    params = pyec.Parameters()
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 15.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    e = float(m.potential_energy)
    assert np.isfinite(e)
    E, F = pot.get_ef(
        np.asarray(m.positions),
        np.asarray(m.atomic_numbers, dtype=np.int64),
        np.asarray(m.cell),
    )
    assert np.isfinite(E)
    assert F.shape == (2, 3)


def test_make_backend_ase_lj():
    pytest.importorskip("ase")
    from ase.calculators.lj import LennardJones

    pot = make_backend(
        "ase",
        calculator=LennardJones(epsilon=1.0, sigma=1.0, rc=5.0, smooth=False),
    )
    pos = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]], dtype=np.float64)
    z = np.array([1, 1], dtype=np.int64)
    cell = np.eye(3) * 12.0
    E, F = pot.get_ef(pos, z, cell)
    assert np.isfinite(E)
    assert F.shape == (2, 3)
    assert F[0, 0] * F[1, 0] < 0  # opposite pair force


def test_unknown_backend_raises():
    with pytest.raises(KeyError, match="unknown backend"):
        make_backend("not-a-real-backend")


def test_make_backend_lj_params_kwarg_sets_pottype():
    params = pyec.Parameters()
    pot = make_backend("lj", params=params)
    assert params.potential == pyec.PotType.LJ
    assert pot is not None
    e, f = pot.get_ef(
        np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]),
        np.array([1, 1], dtype=np.int64),
        np.eye(3) * 12.0,
    )
    assert np.isfinite(e)
    assert f.shape == (2, 3)


def test_matter_forces_free_and_max_force_lj():
    """GIL-safe free-force path: eval under release, numpy with GIL held."""
    params = pyec.Parameters()
    pot = make_backend("lj", params=params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 15.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    m.fixed = np.zeros(2, dtype=np.int64)
    ff = np.asarray(m.forces_free)
    assert ff.shape == (2, 3)
    assert np.isfinite(ff).all()
    assert np.isfinite(float(m.max_force))
