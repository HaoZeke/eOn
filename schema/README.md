# eOn Cap'n Proto parameter SSoT

`eon_params.capnp` is the **sole authoring surface** for covered simulation
parameter names, types, defaults, and sections (L0 field graph).

## Covered groups
Main, Potential, Optimizer (+ LBFGS/CG/Quickmin/SD), Structure Comparison,
Process Search.

## Regenerate consumers
```bash
python tools/params_ssot/codegen.py
```
Writes:
- `schema/eon_params_catalog.json`
- `eon/_params_ssot_catalog.py`
- `client/generated/ParametersSSOTDefaults.h`

## Adapters
- Client: `Parameters::load` / JSON (INI/JSON files)
- Python: `eon.config` / Pydantic (`eon/schema.py`) — field names parity-gated

## Do not
- Add keys to covered sections of `config.yaml` or `schema.py` without updating
  this Cap'n Proto file and regenerating.
