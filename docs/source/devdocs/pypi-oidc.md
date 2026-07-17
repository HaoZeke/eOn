---
myst:
  html_meta:
    description: "PyPI Trusted Publishing (OIDC) setup for eon-akmc, pyeonclient, and eon-schema."
    keywords: "eOn PyPI, OIDC, Trusted Publishing, eon-akmc, pyeonclient, eon-schema"
---

# PyPI Trusted Publishing (OIDC)

All three distribution names use **OIDC Trusted Publishing** via
`pypa/gh-action-pypi-publish` and GitHub Environments. Prefer this over
long-lived `PYPI_API_TOKEN` secrets.

## Projects and workflows

| PyPI project | GitHub Environment | Workflow | Trigger |
|--------------|--------------------|----------|---------|
| `eon-akmc` | `release` | `release.yml` | tag `v*` (stable only) |
| `pyeonclient` | `pypi-pyeonclient` | `pyeonclient-wheels.yml` | tag `pyeonclient-v*` or dispatch `publish=true` |
| `eon-schema` | `pypi-eon-schema` | `eon-schema-publish.yml` | tag `eon-schema-v*` or dispatch `publish=true` |

Do **not** create a PyPI project named `eon` (occupied by Epidemics on Networks).

## One-time setup on pypi.org (per project)

For each project above, under **Publishing → Trusted publishers → Add a new
pending publisher**:

| Field | Value |
|-------|--------|
| Owner | `TheochemUI` (or your fork owner) |
| Repository | `eOn` |
| Workflow name | exact file name from table (e.g. `eon-schema-publish.yml`) |
| Environment name | exact environment from table |

Save. The first successful workflow run activates the publisher.

## One-time setup on GitHub

Create Environments (Settings → Environments):

- `release` — optional required reviewers for fat/eon-akmc cut  
- `pypi-pyeonclient`  
- `pypi-eon-schema`  

Workflow permissions: `id-token: write` is set on the publish jobs.

Optional: keep `PYPI_API_TOKEN` only as emergency fallback; OIDC jobs do not
need it when Trusted Publishing is active.

## Publish order (2.17 train)

1. `eon-schema-v0.2.1` → pure Python  
2. `pyeonclient-v0.3.2` → abi3 base wheels  
3. `v2.17.0` → fat tarball + eon-akmc sdist (then conda-forge feedstock)

## Verify without publish

- PR builds for package path changes  
- `workflow_dispatch` with `publish=false`  
- Tag push is production once Trusted Publishers are configured
