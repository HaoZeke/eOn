# Publishing: fat monorepo tarball vs split PyPI packages

eOn is a **monorepo**. Releases can ship **both**:

| Shape | Artifact | Consumer |
|-------|----------|----------|
| **Fat tree** | `eon-vX.Y.Z.tar.xz` (`git archive` of the monorepo) | **conda-forge `eon-feedstock`**, full source builds, EasyBuild |
| **Split** | PyPI `eon-akmc`, `pyeonclient`, `eon-schema`, … | pip/uv pure-Python or focused wheels |

Layout cleanup (where Cap’n Proto lives, package paths, etc.) is allowed.
The release contract is: **one complete fat archive for the conda recipe**,
even when you also publish split PyPI packages from the same monorepo.

conda-forge does **not** need a separate `eon-schema` feedstock for 0.1.x:
the fat tarball already contains whatever the monorepo needs to build `eon`.

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

Feedstock `recipe.yaml` continues to point at:

```yaml
url: https://github.com/TheochemUI/eOn/releases/download/v${{ version }}/eon-v${{ version }}.tar.xz
```

No new conda package is required when you only add PyPI splits.

### Cap’n Proto SSoT and the fat archive

- **Authoring:** `schema/eon_params.capnp`
- After edits:

  ```bash
  python tools/params_ssot/codegen.py
  ./packages/eon-schema/scripts/sync_ssot_into_package.sh
  ```

- Commit generated catalogs / headers and the package vendored copy so
  `git archive` is self-contained for feedstock **and** for split builds.

Do not `export-ignore` paths the fat build needs (e.g. `schema/`, `client/`,
`eon/`, `packages/` if present) in `.gitattributes`.

## Split: PyPI `eon-schema` (independent semver)

Independent of a `v2.y.z` cut. Bump `packages/eon-schema/pyproject.toml`
`version` when the public API or vendored SSoT changes.

### Build + check

```bash
# from monorepo root if SSoT changed:
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh

cd packages/eon-schema
python -m build
twine check dist/*
python -m pip install dist/eon_schema-*.whl 'pydantic>=2'
python -c "from eon_schema.ssot import capnp_path; from eon_schema.api import DimerSpec; print(capnp_path()); print(DimerSpec())"
```

### Upload

```bash
twine upload dist/*
# or Trusted Publishing / uv publish
```

Optional monorepo tag for archaeology only: `eon-schema-v0.1.0`  
(do **not** use this tag as the feedstock source — feedstock wants the **fat**
`eon-v*` archive).

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

## Checklist: full eOn release (fat + optional splits)

- [ ] SSoT edits under `schema/` if needed; codegen + package sync committed
- [ ] towncrier / cog / `release_assert.py`
- [ ] Tag `vX.Y.Z` → fat `eon-v*.tar.xz` → feedstock sha256
- [ ] Optional same-day: bump/publish split packages if their public API changed

## Checklist: eon-schema-only release

- [ ] Version bump in `packages/eon-schema/pyproject.toml`
- [ ] `CHANGELOG.md` entry
- [ ] Vendored SSoT in sync with monorepo `schema/`
- [ ] Build + twine + upload to PyPI
- [ ] **No** feedstock PR unless you intentionally add `eon-schema` conda later
