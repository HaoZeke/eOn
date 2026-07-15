---
myst:
  html_meta:
    "description": "RGPOT metatomic backend: optional libmetatomic_engine.so vs fat potential = Metatomic."
---

# RGPOT Metatomic (optional engine)

Two ways to run Metatomic models in eOn. **Both stay supported.** rgpot is not
on conda-forge, so the fat native Metatomic path is the packaging default;
the RGPOT engine is an optional thin-host path that loads the same science via
`libmetatomic_engine.so`.

| Path | Config | Torch on host `eonclient`? | Use when |
|------|--------|----------------------------|----------|
| **Fat / native** | `potential = Metatomic` + `[Metatomic]` | Yes (`-Dwith_metatomic`) | Default, fastest, conda-forge |
| **RGPOT / engine** | `potential = RGPOT`, `backend = metatomic` | No on thin hosts (`with_metatomic=false`); engine `.so` holds torch | Optional plugin; same energies |

```{code-block} ini
[Potential]
potential = RGPOT

[RgpotPot]
backend = metatomic
engine_path = /path/to/libmetatomic_engine.so
model_path = /path/to/model.pt
device = cpu
```

Or set `RGPOT_METATOMIC_ENGINE` / `METATOMIC_ENGINE` and put the model under
`[Metatomic] model_path` (also accepted when backend is metatomic).

## Build

```{code-block} bash
# Fat host: native Metatomic pot + engine for optional RGPOT consumers
meson setup build-mta -Dwith_metatomic=true -Dwith_rgpot=true
meson compile -C build-mta
# -> client/libmetatomic_pot.so  (potential = Metatomic)
# -> client/libmetatomic_engine.so  (RGPOT backend=metatomic)

# Thin host: RGPOT only; no torch at link time
meson setup build-thin -Dwith_metatomic=false -Dwith_rgpot=true
meson compile -C build-thin
# engine still comes from a fat build or a packager that ships libmetatomic_engine.so
```

Do **not** remove native Metatomic from eOn: it is the conda-forge path.

## Runtime notes

- Prefer **build-tree pot plugins** over any env that also installs
  `librgpot_pot.so` (e.g. cookbook `.nox/.../lib`). Putting a fat install prefix
  first on `LD_LIBRARY_PATH` can interpose a pre-metatomic `librgpot_pot` and
  yield `unknown backend 'metatomic'`.
- The engine needs torch **when dlopened** (`LD_LIBRARY_PATH` / rpath of the
  engine and its metatensor stack). Thin host DT_NEEDED stays free of torch.
- The engine wraps eOn `MetatomicPotential` through the C ABI
  (`rgpot_mta_create` / `rgpot_mta_force`); energies match fat
  `potential = Metatomic` on the same model (verified to 1e-3 eV on PET-MAD).
