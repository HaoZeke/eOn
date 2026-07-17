"""Class-first jobs + inplace=False default (no workdir)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _lj_dimer(seed: int = 0):
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.random_seed = seed
    params.temperature = 300.0
    params.opt_max_iterations = 200
    params.opt_converged_force = 0.1
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    m.fixed = np.zeros(2, dtype=np.int64)
    m.periodic = False
    return m, pot, params


def _short_md(params):
    params.dynamics_steps = 4
    params.dynamics_time_step = 1.0
    return params


def test_relax_default_does_not_mutate():
    m, pot, params = _lj_dimer()
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] = 1.6
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    out, ok = m.relax(inplace=False)
    assert np.allclose(before, np.asarray(m.positions))


def test_relax_inplace_mutates():
    m, pot, params = _lj_dimer()
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] = 1.6
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    out, ok = m.relax(inplace=True)
    assert not np.allclose(before, np.asarray(m.positions))


def test_molecular_dynamics_non_mutating():
    m, pot, params = _lj_dimer(1)
    _short_md(params)
    before = np.asarray(m.positions).copy()
    out = pyec.MolecularDynamics(m, params, pot).run(inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert np.asarray(out.positions).shape == before.shape


def test_monte_carlo_non_mutating():
    m, pot, params = _lj_dimer(2)
    params.monte_carlo_steps = 5
    params.monte_carlo_step_size = 0.01
    before = np.asarray(m.positions).copy()
    out = pyec.MonteCarlo(m, params, pot).run(inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert np.asarray(out.positions).shape == before.shape


def test_basin_hopping_non_mutating():
    m, pot, params = _lj_dimer(3)
    params.basin_hopping_steps = 2
    params.basin_hopping_displacement = 0.05
    before = np.asarray(m.positions).copy()
    out = pyec.BasinHopping(m, params, pot).run(inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert np.asarray(out.positions).shape == before.shape


def test_process_search_non_mutating():
    m, pot, params = _lj_dimer(4)
    params.process_search_minimize_first = False
    params.saddle_max_iterations = 10
    params.saddle_displace_magnitude = 0.1
    params.dimer_rotations_max = 2
    mode = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    before = np.asarray(m.positions).copy()
    reactant, saddle, status = pyec.ProcessSearch(
        m, mode, params, pot
    ).run(inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert isinstance(status, int)


@pytest.mark.parametrize(
    "Cls",
    [
        "TAD",
        "ParallelReplica",
        "SafeHyperdynamics",
        "ReplicaExchange",
    ],
)
def test_long_timescale_class_run_non_mutating(Cls, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    m, pot, params = _lj_dimer(10)
    _short_md(params)
    before = np.asarray(m.positions).copy()
    out = getattr(pyec, Cls)(m, params, pot).run(inplace=False)
    assert np.allclose(before, np.asarray(m.positions)), Cls
    assert np.asarray(out.positions).shape == before.shape
    assert np.isfinite(float(out.potential_energy))


@pytest.mark.parametrize("Cls", ["TAD", "ParallelReplica"])
def test_long_timescale_inplace(Cls, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    m, pot, params = _lj_dimer(11)
    _short_md(params)
    out = getattr(pyec, Cls)(m, params, pot).run(inplace=True)
    assert np.asarray(m.positions).shape == (2, 3)
    assert np.isfinite(float(m.potential_energy))


def test_structure_comparison():
    m1, pot, params = _lj_dimer(6)
    m2 = pyec.Matter(m1)
    assert pyec.structures_equal(m1, m2)
    assert pyec.structure_distance(m1, m2) == pytest.approx(0.0, abs=1e-9)


def test_neb_accelerant_gp_gated():
    m1, pot, params = _lj_dimer(7)
    m2 = pyec.Matter(m1)
    pos = np.asarray(m2.positions).copy()
    pos[1, 0] += 0.3
    m2.positions = pos
    if not pyec.built_with_gp_surrogate():
        with pytest.raises(RuntimeError, match="gp_surrogate|WITH_GP|accelerant"):
            pyec.NEB(m1, m2, params, pot, accelerant="gp").compute()
    # plain NEB construct still works
    params.neb_images = 2
    params.neb_init_method = pyec.NEBInit.LINEAR
    params.neb_minimize_endpoints = False
    params.opt_max_iterations = 5
    params.opt_converged_force = 1.0
    neb = pyec.NEB(m1, m2, params, pot)
    # compute may not converge; construct is enough if force eval works
    st = neb.compute()
    assert neb.n_path >= 2


def test_no_run_star_free_functions():
    for name in (
        "run_tad",
        "run_dynamics",
        "run_monte_carlo",
        "run_basin_hopping",
        "run_gp_surrogate_neb",
        "run_parallel_replica",
    ):
        assert not hasattr(pyec, name), f"{name} must not be public"


def test_public_class_symbols():
    for name in (
        "Dimer",
        "NEB",
        "TAD",
        "MolecularDynamics",
        "MonteCarlo",
        "BasinHopping",
        "ProcessSearch",
        "NudgedElasticBand",
        "Hessian",
        "built_with_gp_surrogate",
        "built_with_gprd",
        "neb_idpp_path",
    ):
        assert hasattr(pyec, name), name
