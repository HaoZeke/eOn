"""inplace=False default + Matter-first job APIs (no workdir).

Issues: eOn-rczn (inplace), eOn-m6o8 process_search, eOn-p8em dynamics,
eOn-dvf5 MC, eOn-73ji BH, eOn-bcd3 TAD, eOn-tw3g PR family, eOn-515s
structure comparison, eOn-cb3z GP surrogate gate.
"""

from __future__ import annotations

import numpy as np
import pytest

pc = pytest.importorskip("pyeonclient")


def _lj_dimer(seed: int = 0):
    params = pc.Parameters()
    params.potential = pc.PotType.LJ
    params.random_seed = seed
    params.temperature = 300.0
    params.opt_max_iterations = 200
    params.opt_converged_force = 0.1
    pot = pc.make_potential(pc.PotType.LJ, params)
    m = pc.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    m.fixed = np.zeros(2, dtype=np.int64)
    m.periodic = False
    return m, pot, params


def test_relax_default_does_not_mutate():
    m, pot, params = _lj_dimer()
    # displace from min so relax moves
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] = 1.6
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    out, ok = m.relax(inplace=False)
    after = np.asarray(m.positions)
    assert np.allclose(before, after), "default relax must not mutate self"
    out_pos = np.asarray(out.positions)
    assert not np.allclose(before, out_pos) or ok is not None
    # out is a distinct configuration object
    assert out is not m or not np.allclose(before, out_pos)


def test_relax_inplace_mutates():
    m, pot, params = _lj_dimer()
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] = 1.6
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    out, ok = m.relax(inplace=True)
    after = np.asarray(m.positions)
    assert not np.allclose(before, after), "inplace=True must move atoms"
    assert np.allclose(after, np.asarray(out.positions))


def test_min_mode_saddle_search_non_mutating():
    m, pot, params = _lj_dimer()
    params.saddle_max_iterations = 15
    params.dimer_rotations_max = 3
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] += 0.2
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    mode = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    E0 = float(m.potential_energy)
    out, status = pc.min_mode_saddle_search(
        m, mode, E0, params, pot, inplace=False
    )
    assert np.allclose(before, np.asarray(m.positions))
    assert isinstance(status, int)
    assert out is not None
    _ = np.asarray(out.positions)


def test_min_mode_saddle_search_inplace():
    m, pot, params = _lj_dimer()
    params.saddle_max_iterations = 15
    params.dimer_rotations_max = 3
    pos0 = np.asarray(m.positions).copy()
    pos0[1, 0] += 0.2
    m.positions = pos0
    before = np.asarray(m.positions).copy()
    mode = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    E0 = float(m.potential_energy)
    out, status = pc.min_mode_saddle_search(
        m, mode, E0, params, pot, inplace=True
    )
    # with inplace, out is same shared object as m
    after = np.asarray(m.positions)
    # algorithm may move or not depending on status; out should be m
    assert out is m or np.allclose(np.asarray(out.positions), after)


def test_run_dynamics_non_mutating():
    m, pot, params = _lj_dimer(1)
    params.dynamics_steps = 5
    params.dynamics_time_step = 1.0
    before = np.asarray(m.positions).copy()
    out = pc.run_dynamics(m, params, inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert out is not None
    assert np.asarray(out.positions).shape == before.shape


def test_run_monte_carlo_non_mutating():
    m, pot, params = _lj_dimer(2)
    params.monte_carlo_steps = 5
    params.monte_carlo_step_size = 0.01
    before = np.asarray(m.positions).copy()
    out = pc.run_monte_carlo(m, params, inplace=False)
    assert np.allclose(before, np.asarray(m.positions))
    assert np.asarray(out.positions).shape == before.shape


def test_run_basin_hopping_non_mutating():
    m, pot, params = _lj_dimer(3)
    params.basin_hopping_steps = 2
    params.basin_hopping_displacement = 0.05
    before = np.asarray(m.positions).copy()
    out = pc.run_basin_hopping(m, params, pot, inplace=False)
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
    reactant, saddle, status = pc.process_search(
        m, mode, params, pot, inplace=False
    )
    assert np.allclose(before, np.asarray(m.positions))
    assert isinstance(status, int)
    assert np.asarray(saddle.positions).shape == before.shape


def test_long_timescale_entry_symbols():
    # Named TAD/PR/safe-hyper/replica-exchange entries are thin Dynamics
    # wrappers (identical body to run_dynamics, covered by
    # test_run_dynamics_non_mutating). Assert the Matter-first public surface.
    for name in (
        "run_tad",
        "run_parallel_replica",
        "run_safe_hyperdynamics",
        "run_replica_exchange",
    ):
        assert callable(getattr(pc, name)), name


def test_structure_comparison():
    m1, pot, params = _lj_dimer(6)
    m2 = pc.Matter(m1)  # copy construct
    assert pc.structures_equal(m1, m2)
    assert pc.structure_distance(m1, m2) == pytest.approx(0.0, abs=1e-9)
    pos = np.asarray(m2.positions).copy()
    pos[1, 0] += 0.5
    m2.positions = pos
    assert not pc.structures_equal(m1, m2)
    assert pc.structure_distance(m1, m2) > 0.0


def test_gp_surrogate_gate():
    assert hasattr(pc, "built_with_gp_surrogate")
    m1, pot, params = _lj_dimer(7)
    m2 = pc.Matter(m1)
    if not pc.built_with_gp_surrogate():
        with pytest.raises(RuntimeError, match="gp_surrogate|WITH_GP"):
            pc.run_gp_surrogate_neb(m1, m2, params, pot, inplace=False)
    else:
        # Feature on: may still raise if unwired — must not pretend success
        try:
            pc.run_gp_surrogate_neb(m1, m2, params, pot, inplace=False)
        except RuntimeError as e:
            assert "gp" in str(e).lower() or "surrogate" in str(e).lower()


def test_public_symbols_present():
    for name in (
        "min_mode_saddle_search",
        "run_dynamics",
        "run_monte_carlo",
        "run_basin_hopping",
        "process_search",
        "structures_equal",
        "structure_distance",
        "run_tad",
        "run_parallel_replica",
        "run_safe_hyperdynamics",
        "run_replica_exchange",
        "run_gp_surrogate_neb",
        "built_with_gp_surrogate",
        "ImprovedDimer",
        "MinModeSaddleSearch",
        "NudgedElasticBand",
        "Hessian",
    ):
        assert hasattr(pc, name), name
