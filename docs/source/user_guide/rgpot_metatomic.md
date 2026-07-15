---
myst:
  html_meta:
    "description": "RGPOT metatomic backend: optional libmetatomic_engine.so vs fat potential = Metatomic."
---

# RGPOT Metatomic (optional engine)

Two ways to run Metatomic models in eOn:

| Path | Config | Torch linked into `eonclient`? | Use when |
|------|--------|--------------------------------|----------|
| **Fat / native** | `potential = Metatomic` + `[Metatomic]` | Yes (`-Dwith_metatomic`) | Default, fastest, conda-forge |
| **RGPOT / engine** | `potential = RGPOT`, `backend = metatomic` | No for thin hosts; engine `.so` holds torch | Optional plugin; same science |

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

Build the engine with `-Dwith_metatomic=true` (produces `libmetatomic_engine.so`
next to the fat pot). Do **not** remove native Metatomic: rgpot is not on
conda-forge, so fat Metatomic remains the packaging path.
