## Install

```bash
# Base wheel (PyPI): RGPOT pot linked; no torch in DT_NEEDED
pip install pyeonclient

# Runtime deps for Metatomic *models* (torch from PyPI — same as torch itself)
pip install 'pyeonclient[metatomic]'

# Optional ASE converters
pip install 'pyeonclient[ase]'
```

Feature probes on a built extension::

    import pyeonclient as pc
    pc.built_with_rgpot()       # True on default PyPI wheels
    pc.built_with_metatomic()   # True only on metatomic wheels / from-source fat builds

### Metatomic-linked wheel (from source, torch-style)

Prebuilt manylinux *base* wheels omit linked Metatomic (size/CUDA matrix).
To link Metatomic (fat pot + engine) against **your** PyPI torch::

```bash
pip install torch metatomic-torch metatensor-torch vesin
pip install -U build nanobind numpy meson ninja meson-python
# from a git checkout of eOn:
PYEONCLIENT_VARIANT=metatomic ./scripts/pyeonclient_build_wheel.sh
pip install dist/pyeonclient-*.whl
```

RGPOT metatomic (dlopen ``libmetatomic_engine.so``) works on **base** wheels when
the engine is built separately and pointed at with ``[RgpotPot] engine_path`` /
``RGPOT_METATOMIC_ENGINE``.


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


## Cookbook / workdir jobs (pip only — no conda-forge eOn)

Base PyPI wheels ship **RGPOT** linked. Metatomic forces use **dlopen** of
``libmetatomic_engine.so`` from the separate ``rgpot`` package (torch **2.7+**
multi-ABI engines). You do **not** need the conda-forge ``eon`` package or the
``eonclient`` binary.

```bash
pip install pyeonclient 'rgpot>=2.4.1' 'pyeonclient[metatomic]'
# or: pip install torch metatomic-torch metatensor-torch metatensor-core vesin
```

Same workdir layout as ``eonclient`` (``config.ini`` + ``pos.con`` / NEB files)::

```python
import pyeonclient as pc
import rgpot

# Optional: pin engine explicitly (defaults via RGPOT_METATOMIC_ENGINE / picker)
eng = rgpot.default_metatomic_engine_path()
pc.run_job_in_directory(".")  # reads config.ini in cwd
# helpers: pc.rgpot_metatomic_workdir(...), pc.rgpot_metatomic_neb_workdir(...)
```

``config.ini`` Metatomic via RGPOT::

```ini
[Potential]
potential = RGPOT
[RgpotPot]
backend = metatomic
model_path = /path/to/model.pt
# engine_path = ...  # optional; otherwise rgpot multi-ABI picker / env
```

CLI helper::

    python scripts/run_eon_workdir_pyeonclient.py /path/to/workdir

A fat ``-Dwith_metatomic=true`` wheel (torch linked into the pot) remains an
optional local/CI artifact; the pip product path is **base wheel + rgpot engines**.

## ClientEON steps (compose in Python)

Do not treat the client as a single black box. Bound steps::

    params = pc.load_parameters("config.ini")   # or Parameters().load(...)
    t0 = pc.steady_clock_now()
    files = pc.run_job(params)                 # make_job + Job.run + drop
    pc.write_potcall_summary()
    pc.append_timing("results.dat", t0)

Minimization via Matter (no Job wrapper)::

    pc.minimize_workdir("min_reactant")  # con2matter → relax → results → potcalls → timing

NEB still uses ``Job.run`` (NudgedElasticBandJob) composed with the same potcall
and timing steps via ``run_job_in_directory`` / ``run_eon_cwd``.
