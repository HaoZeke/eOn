# pyeonclient on PyPI (standalone package)

`pyeonclient` is a **separate** PyPI project from `eon-akmc`.

| Project | Import | Contents |
|---------|--------|----------|
| `eon-akmc` | `eon` | Server (AKMC); no client binary |
| `pyeonclient` | `pyeonclient` | Nanobind Matter/Parameters/Potential/Jobs |

## Wheel matrix (exactly three manylinux wheels)

| Artifact | Tag | Notes |
|----------|-----|--------|
| Stable ABI | `cp312-abi3-manylinux_*` | Built on 3.12; works on GIL CPython ≥ 3.12 |
| Free-threaded 3.13 | `cp313t-manylinux_*` | `Py_GIL_DISABLED` |
| Free-threaded 3.14 | `cp314t-manylinux_*` | `Py_GIL_DISABLED` |

Meson selects abi3 vs `NB_FREE_THREADED` from the build interpreter
(`sysconfig.Py_GIL_DISABLED`). There is no free-threaded stable ABI yet.

## Source layout

| File | Role |
|------|------|
| `pyproject-pyeonclient.toml` | Package identity for `pyeonclient` |
| `pyproject.toml` | Unchanged `eon-akmc` identity |
| `scripts/pyeonclient_prepare_pyproject.sh` | Swap root pyproject for builds |
| `.github/workflows/pyeonclient-wheels.yml` | Build 3 wheels + sdist; publish |

## Release

1. Bump `version` in `pyproject-pyeonclient.toml` and
   `client/python/bind/module.cpp` `__version__` (keep lockstep).
2. Tag: `git tag -s pyeonclient-v0.2.0 -m "pyeonclient 0.2.0"`
3. Push tag → workflow builds wheels and publishes via Trusted Publishing.

### One-time PyPI Trusted Publisher setup

On [pypi.org](https://pypi.org):

1. Create project **`pyeonclient`** (or first upload creates it).
2. Publishing → Trusted publishers → GitHub:
   - Owner: `TheochemUI` (or fork)
   - Repository: `eOn`
   - Workflow: `pyeonclient-wheels.yml`
   - Environment: `pypi-pyeonclient`

Also create the GitHub Environment `pypi-pyeonclient` (protection rules as desired).

## Local build

```bash
./scripts/pyeonclient_prepare_pyproject.sh
python -m pip install -U build nanobind numpy meson ninja
python -m build   # sdist + wheel for *this* interpreter
./scripts/pyeonclient_prepare_pyproject.sh restore
```

## Install

```bash
pip install pyeonclient           # abi3 or matching free-threaded wheel
pip install 'pyeonclient[ase]'    # optional ASE converters
```
