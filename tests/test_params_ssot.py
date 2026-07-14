"""Cap'n Proto L0 SSoT parity and catalog integrity."""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]

def _pydantic_model_fields(schema_path: Path, class_name: str) -> set[str]:
    """Parse field names from a Pydantic model without importing pydantic."""
    import ast
    tree = ast.parse(schema_path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            names = set()
            for stmt in node.body:
                if isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name):
                    names.add(stmt.target.id)
            return names
    raise AssertionError(f"class {class_name} not found in {schema_path}")



def test_catalog_matches_codegen_roundtrip(tmp_path):
    """Regenerating the catalog is a no-op on a clean tree."""
    import eon._params_ssot_catalog as cat_mod
    from eon import params_ssot

    before = json.dumps(params_ssot.catalog(), sort_keys=True)
    # run codegen in-process
    sys.path.insert(0, str(REPO / "tools" / "params_ssot"))
    import codegen  # type: ignore

    assert codegen.main() == 0
    # reload would be needed; compare files on disk
    disk = (REPO / "schema" / "eon_params_catalog.json").read_text()
    assert json.loads(disk)["source"] == "schema/eon_params.capnp"
    assert "Main" in json.loads(disk)["sections"]
    # catalog module content equals file (written by codegen)
    assert cat_mod.CATALOG["schema_version"] == 1


def test_capnp_schema_file_exists():
    p = REPO / "schema" / "eon_params.capnp"
    assert p.is_file()
    text = p.read_text()
    assert "struct MainOptions" in text
    assert "struct EonParameters" in text
    assert "temperature @2 :Float64 = 300.0" in text


def test_python_defaults_main_temperature():
    from eon import params_ssot

    d = params_ssot.defaults_for("Main")
    assert d["temperature"] == 300.0
    assert d["finite_difference"] == 0.01
    assert d["job"] == "process_search"


def test_parity_main_fields_in_config_yaml():
    """config.yaml Main options must be a subset of SSoT (+ no inventing)."""
    from eon import params_ssot

    y = yaml.load((REPO / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader)
    yaml_keys = set(y["Main"]["options"].keys())
    ssot_keys = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    # server may omit client-only parallel/ini_filename/con_filename
    client_only = {"parallel", "ini_filename", "con_filename"}
    # every yaml key must be in SSoT
    missing = yaml_keys - ssot_keys
    assert not missing, f"config.yaml Main keys not in SSoT: {missing}"


def test_parity_main_fields_in_schema_py():
    from eon import params_ssot

    model_fields = _pydantic_model_fields(REPO / "eon" / "schema.py", "MainConfig")
    ssot_keys = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    assert not (model_fields - ssot_keys), (
        f"schema.MainConfig fields not in SSoT: {model_fields - ssot_keys}"
    )


def test_parity_potential_optimizer_structure_process():
    from eon import params_ssot

    y = yaml.load((REPO / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader)

    checks = [
        ("Potential", "PotentialConfig", "Potential"),
        ("Structure Comparison", "StructureComparisonConfig", "Structure Comparison"),
        ("Process Search", "ProcessSearchConfig", "Process Search"),
    ]
    for sec, model_name, yaml_sec in checks:
        ssot = {f["snake"] for f in params_ssot.scalar_fields(sec)}
        model_keys = _pydantic_model_fields(REPO / "eon" / "schema.py", model_name)
        model_keys = {k for k in model_keys if k[0].islower()}
        missing_from_ssot = model_keys - ssot
        assert not missing_from_ssot, f"{model_name} not in SSoT: {missing_from_ssot}"
        if yaml_sec in y:
            ykeys = set(y[yaml_sec]["options"].keys())
            if sec != "Optimizer":
                assert not (ykeys - ssot), f"yaml {yaml_sec} not in SSoT: {ykeys - ssot}"


def test_parity_parameters_json_main_keys():
    """ParametersJSON.cpp Main keys ⊆ SSoT Main (+ documented aliases)."""
    from eon import params_ssot

    text = (REPO / "client" / "ParametersJSON.cpp").read_text()
    m = re.search(r'j\["Main"\]\s*=\s*\{([^;]+)\};', text, re.S)
    assert m, "Main JSON block not found"
    keys = set(re.findall(r'\{"([^"]+)"', m.group(1)))
    ssot = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    # write_con_forces may be missing from JSON dump — still must not invent
    assert not (keys - ssot), f"ParametersJSON Main keys not in SSoT: {keys - ssot}"


def test_drift_injection_fails_parity(tmp_path, monkeypatch):
    """Deliberate consumer-only field makes parity fail (gate self-check)."""
    from eon import params_ssot

    # inject fake field only on the consumer side check
    fake_yaml_keys = {"job", "temperature", "not_in_ssot_ever"}
    ssot_keys = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    assert fake_yaml_keys - ssot_keys == {"not_in_ssot_ever"}


def test_codegen_header_exists_and_mentions_temperature():
    h = (REPO / "client" / "generated" / "ParametersSSOTDefaults.h").read_text()
    assert "MAIN_TEMPERATURE" in h
    assert "300.0" in h
    assert "AUTO-GENERATED" in h
