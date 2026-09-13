---
myst:
  html_meta:
    "description": "eOn v3.2.1: IDPP path init honors frozen atoms; pyeonclient can set optimizer and refine methods on NEB and dimer jobs."
    "keywords": "eOn v3.2.1, IDPP constraints, pyeonclient opt_method, NEB refine"
---

## [v3.2.1] - 2026-09-13

Patch on `v3.2.0`. IDPP, collective IDPP, and SIDPP no longer move
whole-atom frozen sites when a mover sits next to them
(TheochemUI/eOn#410). pyeonclient `Parameters` can set `opt_method`,
`neb_opt_method`, `refine_opt_method`, and `refine_threshold`
(TheochemUI/eOn#406).

The version chore (`cog bump`) still has to run on `main` after this
branch merges. These notes exist so that cut is not blocked on writing.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
