"""Dimer/NEB validation: light rules always; DimerSpec needs pydantic optional."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _lj():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    m.cell = np.eye(3) * 20
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    m.fixed = np.zeros(2, dtype=np.int64)
    return m, pot, params


def test_light_gp_requires_improved_without_pydantic_spec():
    m, pot, params = _lj()
    # ValueError (light) or ValidationError (if pydantic installed)
    exc: tuple = (ValueError,)
    try:
        from pydantic import ValidationError

        exc = (ValueError, ValidationError)
    except ImportError:
        pass
    with pytest.raises(exc, match="improved"):
        pyec.Dimer(m, params, pot, method="classic", accelerant="gp")
    with pytest.raises(exc, match="improved"):
        pyec.Dimer(m, params, pot, method="lanczos", accelerant="gp")


def test_dimer_default_works_without_spec():
    m, pot, params = _lj()
    d = pyec.Dimer(m, params, pot)
    d.compute(m, np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]))
    assert np.isfinite(d.eigenvalue)


@pytest.mark.skipif(
    pyec.DimerSpec is None,
    reason="optional: pip/uv install \'pyeonclient[models]\' (pydantic)",
)
def test_dimer_spec_gp_with_classic_fails():
    from pydantic import ValidationError

    with pytest.raises(ValidationError, match="improved"):
        pyec.DimerSpec(method="classic", accelerant="gp")


@pytest.mark.skipif(pyec.DimerSpec is None, reason="optional pydantic models")
def test_dimer_spec_defaults_and_gp_improved():
    s = pyec.DimerSpec()
    assert s.method.value == "improved"
    s2 = pyec.DimerSpec(method="improved", accelerant="gp")
    assert s2.core_kwargs()["accelerant"] == "gp"


@pytest.mark.skipif(pyec.NebSpec is None, reason="optional pydantic models")
def test_neb_spec_path_init():
    s = pyec.NebSpec(path_init="idpp", n_images=5)
    params = pyec.Parameters()
    s.apply_to_parameters(params)
    assert params.neb_images == 5
    assert params.neb_init_method == pyec.NEBInit.IDPP


@pytest.mark.skipif(pyec.NebSpec is None, reason="optional pydantic models")
def test_neb_spec_energy_weighted_ci_mmf_and_file():
    """NebSpec maps cookbook knobs onto Parameters (shipped apply_to_parameters)."""
    s = pyec.NebSpec(
        n_images=8,
        path_init="file",
        path_list="idppPath.dat",
        energy_weighted=True,
        ew_ksp_min=0.9,
        ew_ksp_max=9.0,
        ci_mmf=True,
        ci_mmf_after=0.1,
        max_iterations=50,
        force_tolerance=0.02,
        max_move=0.05,
        write_movies=True,
        random_seed=42,
    )
    params = pyec.Parameters()
    out = s.apply_to_parameters(params)
    assert out is params
    assert params.neb_images == 8
    assert params.neb_init_method == pyec.NEBInit.FILE
    assert params.neb_initial_path == "idppPath.dat"
    assert params.neb_energy_weighted is True
    assert abs(params.neb_ew_ksp_min - 0.9) < 1e-12
    assert abs(params.neb_ew_ksp_max - 9.0) < 1e-12
    assert params.neb_ci_mmf is True
    assert params.neb_max_iterations == 50
    assert abs(params.neb_force_tolerance - 0.02) < 1e-12
    assert abs(params.opt_max_move - 0.05) < 1e-12
    assert params.write_movies is True
    assert params.random_seed == 42


def test_make_backend_lj_shares_params_potential_type():
    """make_backend('lj') with shared Parameters keeps pot type consistent."""
    from pyeonclient.backends import make_backend

    params = pyec.Parameters()
    pot = make_backend("lj", params=params)
    assert params.potential == pyec.PotType.LJ
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.25, 0.0, 0.0]])
    m.cell = np.eye(3) * 18.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    e = float(m.potential_energy)
    assert np.isfinite(e)
