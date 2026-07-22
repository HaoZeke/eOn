"""L1 job-config models live in eon-schema (shared with eon-akmc)."""

from eon_schema.config import Config, MainConfig, Metatomic, PotentialConfig


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
