"""Job factory / run_job_in_directory (skip without pyeonclient)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pc = pytest.importorskip("pyeonclient")


def test_job_type_names():
    assert pc.job_type_from_name("minimization") == pc.JobType.Minimization
    assert "Minimization" in pc.job_type_name(pc.JobType.Minimization)


def test_make_job_minimization():
    params = pc.Parameters()
    params.job = pc.JobType.Minimization
    params.potential = pc.PotType.LJ
    params.quiet = True
    job = pc.make_job(params)
    assert job.get_type() == pc.JobType.Minimization


def test_run_job_in_directory_minimize(tmp_path, pot_params_files):
    """Full client MinimizationJob via cwd (same as eonclient)."""
    work = pot_params_files
    params = pc.Parameters()
    params.job = pc.JobType.Minimization
    params.potential = pc.PotType.LJ
    params.quiet = True
    params.write_log = False
    files = pc.run_job_in_directory(str(work), params)
    assert any("min.con" in f or f.endswith("min.con") for f in files) or (
        work / "min.con"
    ).is_file()
    assert (work / "min.con").is_file() or any("min" in f for f in files)


@pytest.fixture
def pot_params_files(tmp_path):
    """Write config.ini + pos.con for a 2-atom LJ minimize."""
    from eon.structure import Structure
    from eon import fileio as fio

    s = Structure(2)
    s.r = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    s.box = np.eye(3) * 20.0
    s.mass = np.ones(2)
    s.free = np.ones(2)
    s.names = ["H", "H"]
    fio.savecon(str(tmp_path / "pos.con"), s)
    (tmp_path / "config.ini").write_text(
        "[Main]\n"
        "job = minimization\n"
        "quiet = true\n"
        "write_log = false\n"
        "[Potential]\n"
        "potential = lj\n"
        "[Optimizer]\n"
        "max_iterations = 200\n"
        "converged_force = 0.1\n"
    )
    return tmp_path
