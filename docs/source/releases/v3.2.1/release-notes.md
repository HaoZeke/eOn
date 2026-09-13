# eOn 3.2.1

A patch for two holes that sat on `main` since 3.2.0.

See monorepo `CHANGELOG.md` section 3.2.1 for the fragment list.

## IDPP path init keeps frozen atoms fixed

`IDPPObjectiveFunction` and `CollectiveIDPPObjectiveFunction` optimized
the full 3N coordinate vector and wrote every atom. Pair forces from a
mover next to a frozen site dragged that site. The objectives now use
the same free-atom DOF contract as `MatterObjectiveFunction`. Linear
and SIDPP interpolation are masked with `getFree()`.

`pytest tests/test_pyeonclient_neb_paths.py` and
`tests/test_pyeonclient_opt_params.py` (7 passed), including IDPP,
collective IDPP, and SIDPP on a three-atom slab with two frozen atoms
and a mover at 0.9 A from one of them.

## pyeonclient can pick the job optimizer

INI already had `[Optimizer]` and Refine. The nanobind `Parameters`
object did not. Dimers and NEBs were stuck on compiled defaults. The
bindings now expose `opt_method`, `neb_opt_method`,
`refine_opt_method`, and `refine_threshold` as `OptType` values.

## Dependency floors

The meson wrap and Python floors move with the current releases:
rgpot **v3.1.2** (ForceBatch is on that tag; do not track a deleted
branch), readcon-core **v0.14.10** / PyPI readcon **0.14.9**, vesin
**0.6.1**, Highway **1.4.0**, potentials-schema **v1.15.1**, xtb
**v6.7.1**.

rgpot 3.1 ships in-process D3 and D4 kernels. eOn now constructs them:
`potential = dftd3` or `dftd4`, with `[D3Pot]` / `[D4Pot]` for functional,
ATM, D3 damping, and D4 charge. `potential = expr` builds `rgpot::ExprPot`
from `[ExprPot] expression` plus a named `terms` list (`lj:lj, d3:d3`).
