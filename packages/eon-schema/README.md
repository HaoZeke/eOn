# eon-schema

Shared **eOn** Cap’n Proto L0 SSoT (**vendored** from monorepo `schema/`),
enums, and optional pydantic **API** models for `eon-akmc` and `pyeonclient`.

## Fat tree vs this package

| | Fat monorepo tarball | This PyPI package |
|--|----------------------|-------------------|
| Purpose | conda-forge `eon`, full C++/server builds | Shared Python models + SSoT copy |
| Cap’n Proto | Full monorepo including `schema/` | Vendored under `src/eon_schema/ssot/` |
| Feedstock | **Uses fat tarball only** | Not required for 0.1.x |

**Authoring** is monorepo `schema/eon_params.capnp`. After edits:

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh
```

Release contract: always produce the **fat** `eon-v*.tar.xz` for
conda-forge; publish **splits** (`eon-schema`, `pyeonclient`, `eon-akmc`, …)
independently when their APIs change. Tree layout may be cleaned up as long as
the fat archive remains buildable.

See **[PUBLISHING.md](PUBLISHING.md)** for release trains and PyPI steps.

## Install

```bash
pip install eon-schema
pip install 'eon-schema[pydantic]'
# uv: uv pip install 'eon-schema[pydantic]'
```

Optional: `RGPYCRUMBS_AUTO_DEPS=1` with `rgpycrumbs` for on-demand pydantic install.

## Quick use

```python
from eon_schema.ssot import capnp_path, load_catalog
from eon_schema.fields import MinModeMethod, Accelerant
from eon_schema.api import DimerSpec  # needs pydantic

print(capnp_path())
print(list(load_catalog()["sections"])[:3])
print(DimerSpec(method=MinModeMethod.improved, accelerant=Accelerant.gp).core_kwargs())
```

## Layout

```text
src/eon_schema/
  ssot/          # vendored Cap'n Proto L0 + catalog (from monorepo schema/)
  fields/        # enums (no pydantic)
  api/           # DimerSpec, NebSpec (pydantic optional)
  config/        # placeholder: L1 still in eon-akmc eon.schema
  _deps.py       # ensure_import (rgpycrumbs if present)
```
