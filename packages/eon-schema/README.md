# eon-schema

Shared **eOn** Cap’n Proto L0 SSoT (**authoring home**), enums, and optional
pydantic **API** models for `eon-akmc` and `pyeonclient`.

## Fat tree vs this package

| | Fat monorepo tarball | This PyPI package |
|--|----------------------|-------------------|
| Purpose | conda-forge `eon`, full C++/server builds | Shared Python models + SSoT copy |
| Cap’n Proto | `schema/` **mirror** in the tarball | Authoring under `src/eon_schema/ssot/` |
| Feedstock | **Uses fat tarball only** | Not required for 0.1.x |

After editing the SSoT here, sync into the monorepo mirror so fat releases stay complete:

```bash
./packages/eon-schema/scripts/sync_ssot_to_tree.sh
python tools/params_ssot/codegen.py
```

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
  ssot/          # Cap'n Proto L0 authoring + catalog
  fields/        # enums (no pydantic)
  api/           # DimerSpec, NebSpec (pydantic optional)
  config/        # placeholder: L1 still in eon-akmc eon.schema
  _deps.py       # ensure_import (rgpycrumbs if present)
```
