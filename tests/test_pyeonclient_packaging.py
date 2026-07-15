"""Packaging metadata for PyPI-able pyeonclient (torch-style, not conda-forge-gated)."""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject-pyeonclient.toml"


def _load_toml():
    try:
        import tomllib
    except ImportError:  # pragma: no cover
        import tomli as tomllib  # type: ignore
    return tomllib.loads(PYPROJECT.read_text())


def test_pyproject_identity_and_extras():
    assert PYPROJECT.is_file()
    data = _load_toml()
    proj = data["project"]
    assert proj["name"] == "pyeonclient"
    # version lockstep with module (when extension present)
    assert proj["version"] >= "0.3.0"
    deps = proj["dependencies"]
    assert any("numpy" in d for d in deps)
    assert any("readcon" in d for d in deps)
    opt = proj["optional-dependencies"]
    assert "metatomic" in opt
    assert any("torch" in d for d in opt["metatomic"])
    assert any("metatomic" in d for d in opt["metatomic"])
    assert "ase" in opt
    # default meson args: RGPOT on, metatomic off (base wheel)
    setup = data["tool"]["meson-python"]["args"]["setup"]
    assert "-Dwith_pyeonclient=true" in setup
    assert "-Dwith_rgpot=true" in setup
    assert "-Dwith_metatomic=false" in setup
    assert "-Dinstall_eon_server=false" in setup


def test_build_wheel_script_exists():
    script = ROOT / "scripts" / "pyeonclient_build_wheel.sh"
    assert script.is_file()
    text = script.read_text()
    assert "PYEONCLIENT_VARIANT" in text
    assert "metatomic" in text
    assert "with_metatomic=true" in text


def test_runtime_feature_flags_when_built():
    pc = pytest.importorskip("pyeonclient")
    assert hasattr(pc, "built_with_metatomic")
    assert hasattr(pc, "built_with_rgpot")
    assert isinstance(pc.built_with_metatomic(), bool)
    assert isinstance(pc.built_with_rgpot(), bool)
    # version lockstep with pyproject
    data = _load_toml()
    assert pc.__version__ == data["project"]["version"]
