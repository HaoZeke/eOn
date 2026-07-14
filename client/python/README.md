## Install

```bash
pip install pyeonclient
pip install 'pyeonclient[ase]'  # optional
```

# pyeonclient

Nanobind-only bindings for the eOn C++ client library (`libeonclib`).

**No pybind11.** This module supersedes the old unused `with_pybind11` option.

## ABI policy

| Interpreter | Extension build |
|-------------|-----------------|
| CPython ≥ 3.12 (GIL) | **Stable ABI (abi3)** — `Py_LIMITED_API=0x030C0000` |
| CPython ≥ 3.13 free-threaded (`Py_GIL_DISABLED`) | **Free-threaded** (`NB_FREE_THREADED`); stable ABI not available yet |

Meson detects free-threaded vs GIL builds from `sysconfig.Py_GIL_DISABLED` and
sets exactly one of the two modes (same rule as nanobind’s CMake
`STABLE_ABI` + `FREE_THREADED` pair).

## Package layout

```
client/python/
  meson.build              # -Dwith_pyeonclient=true
  pyeonclient/
    __init__.py            # re-exports
    bridge.py              # Structure <-> Matter (server)
  bind/
    module.cpp             # NB_MODULE(_core, ...)
    bind_*.cpp
    eigen_numpy.hpp
```

Artifact: `pyeonclient/_core*.so` (or `_core*.abi3.so`).

## Build

```bash
# needs nanobind>=2.2 on the build Python
meson setup build -Dwith_pyeonclient=true
meson compile -C build
pytest tests/test_pyeonclient_matter.py -v
```

## Design

Bottom-up surface: **Matter**, **Parameters**, **Potential**, enums/I/O —
not a thin job wrapper. Jobs and communicators layer on Matter later.
`.con` I/O uses ConFileIO / readcon-core (same path as `eonclient`).
