"""Cap'n Proto L0 SSoT parity, defaults, regen no-op, and real drift gate."""
from __future__ import annotations

import ast
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
GENERATED = [
    REPO / "schema" / "eon_params_catalog.json",
    REPO / "eon" / "_params_ssot_catalog.py",
    REPO / "client" / "generated" / "ParametersSSOTDefaults.h",
    REPO / "client" / "generated" / "ParametersSSOTFieldIndex.inc",
]


def _file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pydantic_model_fields(schema_path: Path, class_name: str) -> set[str]:
    tree = ast.parse(schema_path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            names = set()
            for stmt in node.body:
                if isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name):
                    names.add(stmt.target.id)
            return names
    raise AssertionError(f"class {class_name} not found in {schema_path}")


def _load_yaml_options(section: str) -> dict:
    y = yaml.load((REPO / "eon" / "config.yaml").read_text(), Loader=yaml.BaseLoader)
    return y[section]["options"]


def test_capnp_schema_file_exists():
    p = REPO / "schema" / "eon_params.capnp"
    assert p.is_file()
    text = p.read_text()
    assert "struct MainOptions" in text
    assert "struct EonParameters" in text
    assert "temperature @2 :Float64 = 300.0" in text


def test_codegen_is_noop_on_committed_artifacts(tmp_path):
    """Regenerating from .capnp must not change committed consumer files."""
    before = {p: _file_sha(p) for p in GENERATED if p.is_file()}
    assert len(before) == len(GENERATED), f"missing generated files: {set(GENERATED) - set(before)}"

    env = dict(**__import__("os").environ)
    env["PYTHONPATH"] = str(REPO) + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    proc = subprocess.run(
        [sys.executable, str(REPO / "tools" / "params_ssot" / "codegen.py")],
        cwd=str(REPO),
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    after = {p: _file_sha(p) for p in GENERATED}
    for p in GENERATED:
        assert before[p] == after[p], f"codegen dirty: {p.relative_to(REPO)}"


def test_python_defaults_main_temperature():
    from eon import params_ssot

    d = params_ssot.defaults_for("Main")
    assert d["temperature"] == 300.0
    assert d["finite_difference"] == 0.01
    assert d["job"] == "process_search"


@pytest.mark.parametrize(
    "section",
    ["Main", "Potential", "Optimizer", "Structure Comparison", "Process Search"],
)
def test_parity_yaml_keys_subset_of_ssot(section):
    from eon import params_ssot

    opts = _load_yaml_options(section)
    allowed = params_ssot.allowed_yaml_keys(section)
    # Main yaml may omit client-only keys that are still in SSoT
    extra = set(opts.keys()) - allowed
    assert not extra, f"config.yaml [{section}] keys not in SSoT/aliases: {extra}"


@pytest.mark.parametrize(
    "section",
    ["Main", "Potential", "Optimizer", "Structure Comparison", "Process Search"],
)
def test_parity_yaml_defaults_match_ssot(section):
    from eon import params_ssot

    opts = _load_yaml_options(section)
    expected = params_ssot.yaml_default_map(section)
    for key, exp in expected.items():
        if key not in opts:
            # client-only keys may be absent from server config.yaml
            if section == "Main" and key in {
                "parallel",
                "ini_filename",
                "con_filename",
            }:
                continue
            if section == "Potential" and key in {"potentials_path"}:
                continue
            # Optimizer nested scalars without flat alias only live nested
            if section == "Optimizer" and key in {
                "memory",
                "inverse_curvature",
                "max_inverse_curvature",
                "auto_scale",
                "angle_reset",
                "distance_reset",
                "no_overshooting",
                "knock_out_max_move",
                "line_search",
                "line_converged",
                "line_search_max_iter",
                "max_iter_before_reset",
                "steepest_descent",
                "alpha",
                "two_point",
            }:
                continue
            pytest.fail(f"yaml [{section}] missing SSoT key {key}")
        if key not in opts:
            continue
        got = params_ssot.normalize_yaml_default(opts[key])
        assert got == exp or (isinstance(exp, float) and float(got) == float(exp)), (
            f"yaml [{section}].{key} default {got!r} != SSoT {exp!r}"
        )


def test_parity_main_fields_in_schema_py():
    from eon import params_ssot

    model_fields = _pydantic_model_fields(REPO / "eon" / "schema.py", "MainConfig")
    ssot_keys = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    assert not (model_fields - ssot_keys), (
        f"schema.MainConfig fields not in SSoT: {model_fields - ssot_keys}"
    )


def test_parity_schema_defaults_main_potential_optimizer():
    """Pydantic defaults for key covered fields match SSoT."""
    # Import only if pydantic present; else parse defaults from AST/source
    from eon import params_ssot

    schema_text = (REPO / "eon" / "schema.py").read_text()
    # lammps_logging default=False
    m = re.search(
        r"lammps_logging:\s*bool\s*=\s*Field\(\s*default=(True|False)",
        schema_text,
    )
    assert m, "lammps_logging Field not found"
    assert m.group(1) == "False", f"schema lammps_logging default {m.group(1)} != SSoT False"

    # OptimizerConfig time_step default=1.0
    m = re.search(
        r"class OptimizerConfig\(BaseModel\):.*?time_step:\s*float\s*=\s*Field\(\s*default=([0-9.]+)",
        schema_text,
        re.S,
    )
    assert m, "OptimizerConfig.time_step not found"
    assert float(m.group(1)) == params_ssot.defaults_for("Optimizer")["time_step"]


def test_parity_potential_structure_process_models():
    from eon import params_ssot

    checks = [
        ("Potential", "PotentialConfig"),
        ("Structure Comparison", "StructureComparisonConfig"),
        ("Process Search", "ProcessSearchConfig"),
        ("Optimizer", "OptimizerConfig"),
    ]
    for sec, model_name in checks:
        ssot = {f["snake"] for f in params_ssot.scalar_fields(sec)}
        model_keys = _pydantic_model_fields(REPO / "eon" / "schema.py", model_name)
        model_keys = {k for k in model_keys if k[0].islower()}
        missing = model_keys - ssot
        assert not missing, f"{model_name} not in SSoT: {missing}"


def test_parity_parameters_json_main_keys():
    from eon import params_ssot

    text = (REPO / "client" / "ParametersJSON.cpp").read_text()
    m = re.search(r'j\["Main"\]\s*=\s*\{([^;]+)\};', text, re.S)
    assert m, "Main JSON block not found"
    keys = set(re.findall(r'\{"([^"]+)"', m.group(1)))
    ssot = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    assert not (keys - ssot), f"ParametersJSON Main keys not in SSoT: {keys - ssot}"


def test_parity_parameters_h_main_members():
    """Parameters.h main_options_t members map into SSoT Main snakes."""
    from eon import params_ssot

    text = (REPO / "client" / "Parameters.h").read_text()
    m = re.search(r"struct main_options_t \{(.*?)\}\s*main_options;", text, re.S)
    assert m, "main_options_t not found"
    body = m.group(1)
    # C++ members: Type name{default}; or Type name;
    members = re.findall(r"\b(?:bool|long|double|int|size_t|std::string|JobType)\s+(\w+)\b", body)
    # Map camelCase/C++ names to snake
    cxx_to_snake = {
        "job": "job",
        "randomSeed": "random_seed",
        "temperature": "temperature",
        "quiet": "quiet",
        "writeLog": "write_log",
        "checkpoint": "checkpoint",
        "iniFilename": "ini_filename",
        "conFilename": "con_filename",
        "finiteDifference": "finite_difference",
        "maxForceCalls": "max_force_calls",
        "removeNetForce": "remove_net_force",
        "parallel": "parallel",
        "writeConForces": "write_con_forces",
    }
    ssot = {f["snake"] for f in params_ssot.scalar_fields("Main")}
    for mem in members:
        snake = cxx_to_snake.get(mem)
        assert snake is not None, f"unmapped main_options member {mem}"
        assert snake in ssot, f"Parameters.h main_options.{mem} not in SSoT"


def test_field_index_inc_is_codegen_output():
    inc = (REPO / "client" / "generated" / "ParametersSSOTFieldIndex.inc").read_text()
    assert "AUTO-GENERATED" in inc
    assert '"Main.temperature"' in inc
    assert '"Optimizer.lbfgs_memory"' in inc  # flat alias


def test_real_drift_injection_fails_then_restores(tmp_path):
    """Mutate config.yaml consumer, assert parity fails, restore file."""
    from eon import params_ssot

    yaml_path = REPO / "eon" / "config.yaml"
    original = yaml_path.read_text()
    backup = tmp_path / "config.yaml.bak"
    backup.write_text(original)

    try:
        # inject illegal key under Main
        poisoned = original.replace(
            "        remove_net_force:\n            kind: boolean\n            default: True\n",
            "        remove_net_force:\n            kind: boolean\n            default: True\n\n"
            "        not_in_ssot_ever:\n            kind: boolean\n            default: false\n",
        )
        if poisoned == original:
            # try lowercase true
            poisoned = original.replace(
                "        remove_net_force:\n            kind: boolean\n            default: true\n",
                "        remove_net_force:\n            kind: boolean\n            default: true\n\n"
                "        not_in_ssot_ever:\n            kind: boolean\n            default: false\n",
            )
        assert poisoned != original, "failed to inject drift into config.yaml"
        yaml_path.write_text(poisoned)

        opts = _load_yaml_options("Main")
        allowed = params_ssot.allowed_yaml_keys("Main")
        extra = set(opts.keys()) - allowed
        assert "not_in_ssot_ever" in extra, f"drift not detected: extra={extra}"
    finally:
        yaml_path.write_text(backup.read_text())

    # restored
    opts2 = _load_yaml_options("Main")
    assert "not_in_ssot_ever" not in opts2


def test_codegen_header_exists_and_mentions_temperature():
    h = (REPO / "client" / "generated" / "ParametersSSOTDefaults.h").read_text()
    assert "MAIN_TEMPERATURE" in h
    assert "300.0" in h
    assert "AUTO-GENERATED" in h
