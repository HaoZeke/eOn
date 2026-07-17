# Publishing: fat full-tree tarball vs split PyPI packages

eOn is a **monorepo**. One git tag can feed several distribution shapes:

| Shape | Artifact | Consumer |
|-------|----------|----------|
| **Fat tree** | `eon-vX.Y.Z.tar.xz` (`git archive` of the monorepo) | **conda-forge `eon-feedstock`**, full source builds, EasyBuild |
| **Split** | PyPI `eon-akmc`, `pyeonclient`, `eon-schema`, … | pip/uv pure-Python or focused wheels |

conda-forge does **not** need a separate `eon-schema` feedstock for 0.1.x: the
fat tarball still contains `schema/`, `eon/`, `client/`, and `packages/`.

## Fat tarball (conda-forge / full eOn cut)

Unchanged process — see monorepo [`docs/source/devdocs/release.md`](../../docs/source/devdocs/release.md):

1. Conventional commits + towncrier + `cog bump` → version lockstep.
2. Tag `vX.Y.Z` (prefer signed).
3. `release.yml` runs roughly:

   ```bash
   git archive --format=tar "vX.Y.Z" | xz -9 > "eon-vX.Y.Z.tar.xz"
   ```

4. GitHub Release attaches that asset; feedstock uses its URL + sha256.
5. PyPI job publishes **`eon-akmc`** (server import `eon`) from the same tag
   when configured — **not** the only Python project in the repo.

**Fat-tree requirements for SSoT**

- Authoring: `packages/eon-schema/src/eon_schema/ssot/eon_params.capnp`
- Mirror into tree for historical paths + tarball layout:

  ```bash
  ./packages/eon-schema/scripts/sync_ssot_to_tree.sh
  python tools/params_ssot/codegen.py
  ```

- Commit **both** package SSoT and `schema/` mirror so `git archive` is complete.

Feedstock `recipe.yaml` continues to point at:

```yaml
url: https://github.com/TheochemUI/eOn/releases/download/v${{ version }}/eon-v${{ version }}.tar.xz
```

No new conda package is required when you only add PyPI splits.

## Split: PyPI `eon-schema` (independent semver)

Independent of a `v2.y.z` cut. Bump `packages/eon-schema/pyproject.toml`
`version` when the public API or vendored SSoT changes.

### Build + check

```bash
cd packages/eon-schema
./scripts/sync_ssot_to_tree.sh   # keep fat mirror in sync if you edited package SSoT
python tools/params_ssot/codegen.py   # from repo root, if field graph changed
cd packages/eon-schema
python -m build
twine check dist/*
# smoke
python -m pip install dist/eon_schema-*.whl 'pydantic>=2'
python -c "from eon_schema.ssot import capnp_path; from eon_schema.api import DimerSpec; print(capnp_path()); print(DimerSpec())"
```

### Upload

```bash
twine upload dist/*
# or Trusted Publishing / uv publish
```

Optional monorepo tag for archaeology only: `eon-schema-v0.1.0`  
(do **not** use this tag as the feedstock source).

### PyPI project setup (first time)

1. Create project **`eon-schema`** on pypi.org (2FA).
2. Preferred: Trusted Publisher for GitHub `TheochemUI/eOn` (or your release
   fork), dedicated workflow or `workflow_dispatch` for this package.
3. Fallback: project-scoped API token → `TWINE_USERNAME=__token__`.

## Split: other PyPI projects (same monorepo)

| Project | How it is cut | Notes |
|---------|----------------|-------|
| `eon-akmc` | Full-tree tag + meson-python / release.yml | Import `eon`; **not** name `eon` on PyPI |
| `pyeonclient` | `pyproject-pyeonclient.toml` + wheel CI | Optional `[models]` → `eon-schema[pydantic]` |
| `eon-schema` | This directory | Zero hard deps; optional pydantic |

## What must stay in the fat tarball

Everything feedstock and in-tree builds need, including after cleanup:

- `schema/` mirror of Cap’n Proto + catalog
- `packages/eon-schema/` (authoring + package sources)
- `eon/` (server + `eon/schema.py` Config)
- `client/` (C++ + generated headers)
- `tools/params_ssot/`
- meson/CMake, pixi, etc.

Do not `export-ignore` `packages/` or `schema/` in `.gitattributes`.

## Checklist: full eOn release (fat)

- [ ] SSoT edited under `packages/eon-schema/.../ssot/` if changed
- [ ] `sync_ssot_to_tree.sh` + `codegen.py` committed
- [ ] towncrier / cog / `release_assert.py`
- [ ] Tag `vX.Y.Z` → fat `eon-v*.tar.xz` → feedstock sha256
- [ ] Optional same-day: bump/publish split packages if their public API changed

## Checklist: eon-schema-only release

- [ ] Version bump in `packages/eon-schema/pyproject.toml`
- [ ] `CHANGELOG.md` entry
- [ ] SSoT sync + build + twine
- [ ] Upload to PyPI
- [ ] **No** feedstock PR unless you intentionally add `eon-schema` conda later

## Later (docs monorepo polish)

Deeper monorepo docs (package map, CI matrix, dependency graph) can live under
`docs/source/devdocs/` without changing this dual fat/split release contract.
