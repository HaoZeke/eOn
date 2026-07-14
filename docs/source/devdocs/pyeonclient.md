# pyeonclient (nanobind client bindings)

`pyeonclient` exposes the eOn **C++ client core** to Python with
[nanobind](https://nanobind.readthedocs.io/) only (no pybind11).

| Type | Role |
|------|------|
| `Parameters` | `config.ini` / JSON parameters |
| `Potential` / `make_potential` | forcefield factory |
| `Matter` | structure, energy/forces, `.con` I/O |

## ABI

- **Stable ABI (abi3)** for normal CPython ≥ 3.12.
- **Free-threaded** (`NB_FREE_THREADED`) when the build interpreter has
  `Py_GIL_DISABLED` (CPython ≥ 3.13t). There is no stable ABI for free-threaded
  CPython yet (nanobind / PEP 803); the build picks free-threaded and skips
  abi3 in that case.

## Build

```bash
meson setup build -Dwith_pyeonclient=true
meson compile -C build
```

Requires `nanobind>=2.2` on the build Python (`pixi.toml` pins it).

## Tests

```bash
pytest tests/test_pyeonclient_matter.py -v
```

## Server bridge

```python
from eon import pyeonclient_bridge as br
# br.to_matter(structure, pot, params)
# br.from_matter(matter)
```
