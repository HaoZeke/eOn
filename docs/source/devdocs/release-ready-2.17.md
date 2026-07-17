---
myst:
  html_meta:
    description: "Pre-authorization checklist for eOn 2.17.0 / eon-schema 0.2.1 / pyeonclient 0.3.2."
---

# Release-ready checklist (2.17 train)

Versions are **synchronized in-tree**. Publish only after explicit authorization.

## Version lockstep (done in tree)

| Surface | Version |
|---------|---------|
| monorepo `pyproject.toml` / `pixi.toml` / towncrier | **2.17.0** |
| `CHANGELOG.md` section | **2.17.0** |
| eon-schema | **0.2.1** |
| pyeonclient pyproject + `module.cpp` | **0.3.2** |
| pyeonclient `[models]` pin | `eon-schema>=0.2.1` |

## Pre-push gates

```bash
python3 scripts/release_assert.py 2.17.0 --require-changelog
cd packages/eon-schema && python -m build && twine check dist/*
# pyeonclient: manylinux wheel on rg.terra / CI (not laptop)
```

## Tag order (after authorize + OIDC pending publishers)

```bash
# 1) eon-schema
git tag -s eon-schema-v0.2.1 -m "eon-schema 0.2.1"
git push origin eon-schema-v0.2.1   # → eon-schema-publish.yml OIDC

# 2) pyeonclient (after schema is on PyPI)
git tag -s pyeonclient-v0.3.2 -m "pyeonclient 0.3.2"
git push origin pyeonclient-v0.3.2  # → pyeonclient-wheels.yml OIDC

# 3) monorepo fat + eon-akmc
git tag -s v2.17.0 -m "eOn 2.17.0"
git push origin v2.17.0             # → release.yml tarball + eon-akmc sdist OIDC
```

## OIDC GitHub Environments to create

- `release` (eon-akmc)
- `pypi-pyeonclient`
- `pypi-eon-schema`

Configure matching **Trusted publishers** on pypi.org — see {doc}`pypi-oidc`.

## Not done until you authorize

- Push branch / tags to GitHub  
- PyPI uploads  
- conda-forge feedstock PR for `eon` 2.17.0  
