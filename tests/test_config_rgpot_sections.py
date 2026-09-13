"""Server config.yaml must accept documented rgpot D3/D4/Expr/MOPAC INI."""

from __future__ import annotations

from pathlib import Path

import yaml

from eon.config import ConfigClass

REPO = Path(__file__).resolve().parents[1]
YAML_PATH = REPO / "eon" / "config.yaml"

EXPR_INI = """\
[Main]
job = minimization
temperature = 300

[Potential]
potential = expr

[ExprPot]
expression = 0.5*lj + d3
terms = lj:lj, d3:d3

[D3Pot]
functional = pbe0
atm = false
damping = zero

[D4Pot]
functional = pbe
atm = true
charge = 1.0

[MOPACPot]
charge = 0
spin = 0
model = 4
engine_path =
"""


def _yaml() -> dict:
    return yaml.load(YAML_PATH.read_text(), Loader=yaml.BaseLoader)


def test_potential_values_include_dftd_expr_mopac():
    values = _yaml()["Potential"]["options"]["potential"]["values"]
    for name in ("dftd3", "dftd4", "expr", "mopac"):
        assert name in values, f"Potential.potential values missing {name}"


def test_d3_d4_expr_mopac_sections_match_client_keys():
    y = _yaml()
    assert set(y["D3Pot"]["options"]) == {"functional", "atm", "damping"}
    assert set(y["D4Pot"]["options"]) == {"functional", "atm", "charge"}
    assert set(y["ExprPot"]["options"]) == {"expression", "terms"}
    assert set(y["MOPACPot"]["options"]) == {
        "charge",
        "spin",
        "model",
        "engine_path",
    }
    assert y["D3Pot"]["options"]["functional"]["default"] == "pbe"
    assert y["D3Pot"]["options"]["damping"]["default"] == "bj"
    assert y["D4Pot"]["options"]["charge"]["default"] == "0.0"
    assert y["MOPACPot"]["options"]["model"]["default"] == "4"


def test_configclass_init_accepts_expr_and_child_sections(tmp_path):
    ini = tmp_path / "config.ini"
    ini.write_text(EXPR_INI)
    cfg = ConfigClass()
    cfg.init(str(ini))
    assert cfg.init_done


def test_configclass_init_accepts_standalone_dftd_and_mopac(tmp_path):
    for pot in ("dftd3", "dftd4", "mopac"):
        ini = tmp_path / f"{pot}.ini"
        ini.write_text(
            "[Main]\njob = minimization\ntemperature = 300\n"
            f"[Potential]\npotential = {pot}\n"
        )
        cfg = ConfigClass()
        cfg.init(str(ini))
        assert cfg.init_done
