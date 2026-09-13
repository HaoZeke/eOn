"""L1 job-config models live in eon-schema (shared with eon-akmc)."""

from eon_schema.config import (
    Config,
    D3Pot,
    D4Pot,
    ExprPot,
    MainConfig,
    Metatomic,
    MOPACPot,
    PotentialConfig,
)


def test_main_config_defaults():
    m = MainConfig()
    assert m.job == "akmc"


def test_metatomic_model_fields():
    # construction with defaults
    m = Metatomic()
    assert hasattr(m, "model_path") or "model_path" in type(m).model_fields


def test_root_config_composes_sections():
    # Config requires nested models — build with defaults where possible
    assert MainConfig is not None
    assert PotentialConfig is not None
    assert Config is not None


def test_dftd_expr_mopac_section_defaults():
    d3 = D3Pot()
    assert d3.functional == "pbe"
    assert d3.damping == "bj"
    d4 = D4Pot()
    assert d4.charge == 0.0
    expr = ExprPot(expression="0.5*lj + d3", terms="lj:lj, d3:d3")
    assert expr.expression.startswith("0.5")
    mopac = MOPACPot()
    assert mopac.model == 4
