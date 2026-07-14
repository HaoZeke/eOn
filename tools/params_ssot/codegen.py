#!/usr/bin/env python3
"""Generate consumer projections from schema/eon_params.capnp (L0 SSoT).

Outputs (relative to repo root):
  schema/eon_params_catalog.json
  eon/_params_ssot_catalog.py
  client/generated/ParametersSSOTDefaults.h

Does not require pycapnp: parses the .capnp text for field ordinals, types,
and defaults. Nested structs become dotted paths (e.g. optimizer.lbfgs.memory).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CAPNP = REPO / "schema" / "eon_params.capnp"

# Map capnp struct name -> catalog section key (INI-ish)
SECTION_MAP = {
    "MainOptions": "Main",
    "PotentialOptions": "Potential",
    "StructureComparisonOptions": "Structure Comparison",
    "ProcessSearchOptions": "Process Search",
    "OptimizerOptions": "Optimizer",
    "OptimizerLbfgsOptions": "Optimizer.LBFGS",
    "OptimizerCgOptions": "Optimizer.CG",
    "OptimizerQuickminOptions": "Optimizer.Quickmin",
    "OptimizerSdOptions": "Optimizer.SD",
}

# camelCase / lower first field names -> snake_case used by INI/JSON/Python
FIELD_SNAKE = {
    "randomSeed": "random_seed",
    "writeLog": "write_log",
    "iniFilename": "ini_filename",
    "conFilename": "con_filename",
    "finiteDifference": "finite_difference",
    "maxForceCalls": "max_force_calls",
    "removeNetForce": "remove_net_force",
    "writeConForces": "write_con_forces",
    "mpiPollPeriod": "mpi_poll_period",
    "lammpsLogging": "lammps_logging",
    "lammpsThreads": "lammps_threads",
    "emtRasmussen": "emt_rasmussen",
    "logPotential": "log_potential",
    "extPotPath": "ext_pot_path",
    "potentialsPath": "potentials_path",
    "distanceDifference": "distance_difference",
    "neighborCutoff": "neighbor_cutoff",
    "checkRotation": "check_rotation",
    "indistinguishableAtoms": "indistinguishable_atoms",
    "energyDifference": "energy_difference",
    "removeTranslation": "remove_translation",
    "useCovalent": "use_covalent",
    "covalentScale": "covalent_scale",
    "bruteNeighbors": "brute_neighbors",
    "minimizeFirst": "minimize_first",
    "minimizationOffset": "minimization_offset",
    "optMethod": "opt_method",
    "convergenceMetric": "convergence_metric",
    "maxIterations": "max_iterations",
    "maxMove": "max_move",
    "convergedForce": "converged_force",
    "timeStep": "time_step",
    "maxTimeStep": "max_time_step",
    "inverseCurvature": "inverse_curvature",
    "maxInverseCurvature": "max_inverse_curvature",
    "autoScale": "auto_scale",
    "angleReset": "angle_reset",
    "distanceReset": "distance_reset",
    "noOvershooting": "no_overshooting",
    "knockOutMaxMove": "knock_out_max_move",
    "lineSearch": "line_search",
    "lineConverged": "line_converged",
    "lineSearchMaxIter": "line_search_max_iter",
    "maxIterBeforeReset": "max_iter_before_reset",
    "steepestDescent": "steepest_descent",
    "twoPoint": "two_point",
    "schemaVersion": "schema_version",
}


def to_snake(name: str) -> str:
    if name in FIELD_SNAKE:
        return FIELD_SNAKE[name]
    # generic camelCase
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def parse_default(raw: str | None, type_name: str):
    if raw is None:
        return None
    raw = raw.strip()
    if type_name == "Bool":
        return raw.lower() == "true"
    if type_name in ("Float64", "Float32"):
        return float(raw)
    if type_name in ("Int64", "Int32", "UInt64", "UInt32", "Int16", "UInt16"):
        return int(raw)
    if type_name == "Text":
        if raw.startswith('"') and raw.endswith('"'):
            return raw[1:-1]
        return raw
    return raw


def parse_capnp(text: str) -> dict:
    """Return {struct: [{name, ordinal, type, default, snake}]}"""
    structs: dict[str, list] = {}
    # strip comments
    lines = []
    for line in text.splitlines():
        if "#" in line:
            # keep string hashes? rare in this file
            if '"' not in line.split("#")[0]:
                line = line.split("#")[0]
        lines.append(line)
    body = "\n".join(lines)

    struct_re = re.compile(
        r"struct\s+(\w+)\s*\{([^}]*)\}", re.S
    )
    field_re = re.compile(
        r"(\w+)\s+@(\d+)\s*:\s*([A-Za-z0-9_.]+)"
        r"(?:\s*=\s*([^;]+))?\s*;",
        re.S,
    )
    for m in struct_re.finditer(body):
        sname = m.group(1)
        if sname == "EonParameters":
            # root only stores nested refs; skip nested field types without defaults
            fields = []
            for fm in field_re.finditer(m.group(2)):
                fname, ord_, typ, default = fm.groups()
                fields.append(
                    {
                        "name": fname,
                        "ordinal": int(ord_),
                        "type": typ,
                        "default": parse_default(default, typ) if default else None,
                        "snake": to_snake(fname),
                    }
                )
            structs[sname] = fields
            continue
        fields = []
        for fm in field_re.finditer(m.group(2)):
            fname, ord_, typ, default = fm.groups()
            # skip nested struct members without scalar default (lbfgs @7 :OptimizerLbfgsOptions)
            if typ in SECTION_MAP or typ.endswith("Options"):
                if default is None and typ[0].isupper():
                    fields.append(
                        {
                            "name": fname,
                            "ordinal": int(ord_),
                            "type": typ,
                            "default": None,
                            "snake": to_snake(fname),
                            "nested": True,
                        }
                    )
                    continue
            fields.append(
                {
                    "name": fname,
                    "ordinal": int(ord_),
                    "type": typ,
                    "default": parse_default(default, typ) if default else None,
                    "snake": to_snake(fname),
                }
            )
        structs[sname] = fields
    return structs


def build_catalog(structs: dict) -> dict:
    sections = {}
    for sname, fields in structs.items():
        if sname == "EonParameters":
            continue
        sec = SECTION_MAP.get(sname)
        if not sec:
            continue
        sections[sec] = {
            "struct": sname,
            "fields": [
                {
                    "snake": f["snake"],
                    "capnp": f["name"],
                    "ordinal": f["ordinal"],
                    "type": f["type"],
                    "default": f["default"],
                    **({"nested": True} if f.get("nested") else {}),
                }
                for f in fields
            ],
        }
    return {
        "source": "schema/eon_params.capnp",
        "schema_version": 1,
        "sections": sections,
    }


def emit_python(catalog: dict) -> str:
    body = json.dumps(catalog, indent=2, sort_keys=True)
    py_body = (
        body.replace(": true", ": True")
        .replace(": false", ": False")
        .replace(": null", ": None")
    )
    return (
        '"""AUTO-GENERATED from schema/eon_params.capnp — do not edit.\n\n'
        "Regenerate: python tools/params_ssot/codegen.py\n"
        '"""\n'
        "from __future__ import annotations\n\n"
        f"CATALOG = {py_body}\n"
    )


def emit_cpp_header(catalog: dict) -> str:
    lines = [
        "// AUTO-GENERATED from schema/eon_params.capnp — do not edit.",
        "// Regenerate: python tools/params_ssot/codegen.py",
        "#pragma once",
        "#include <cstdint>",
        "#include <string_view>",
        "namespace eonc::params_ssot {",
        "struct GeneratedDefaults {",
    ]

    def cpp_val(v, typ: str) -> str:
        if v is None:
            return "/* nested */"
        if typ == "Bool":
            return "true" if v else "false"
        if typ == "Text":
            return f'std::string_view{{"{v}"}}'
        if typ in ("Float64", "Float32"):
            return f"{float(v)}"
        return str(int(v))

    # Flat constants for Main/Potential etc
    for sec, data in catalog["sections"].items():
        safe = re.sub(r"[^A-Za-z0-9]+", "_", sec)
        for f in data["fields"]:
            if f.get("nested") or f["default"] is None:
                continue
            cname = f"{safe}_{f['snake']}".upper()
            lines.append(
                f"  static constexpr auto {cname} = {cpp_val(f['default'], f['type'])};"
            )
    lines.append("};")
    lines.append("} // namespace eonc::params_ssot")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    text = CAPNP.read_text(encoding="utf-8")
    structs = parse_capnp(text)
    if "MainOptions" not in structs:
        print("error: failed to parse MainOptions from", CAPNP, file=sys.stderr)
        return 1
    catalog = build_catalog(structs)

    out_json = REPO / "schema" / "eon_params_catalog.json"
    out_py = REPO / "eon" / "_params_ssot_catalog.py"
    out_h = REPO / "client" / "generated" / "ParametersSSOTDefaults.h"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_h.parent.mkdir(parents=True, exist_ok=True)

    out_json.write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")
    out_py.write_text(emit_python(catalog))
    out_h.write_text(emit_cpp_header(catalog))
    print("wrote", out_json.relative_to(REPO))
    print("wrote", out_py.relative_to(REPO))
    print("wrote", out_h.relative_to(REPO))
    nfields = sum(
        len([f for f in s["fields"] if not f.get("nested")])
        for s in catalog["sections"].values()
    )
    print(f"sections={len(catalog['sections'])} scalar_fields={nfields}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
