---
myst:
  html_meta:
    "description": "pyeonclient Python API: Matter, min-mode (Dimer/Lanczos/Davidson), MinModeSaddleSearch, NEB, Hessian, Prefactor."
    "keywords": "pyeonclient, Dimer, saddle search, NEB, Hessian, Matter, ASE"
---

# pyeonclient (Python client API)

```{versionadded} 2.16
In-process client: first-class algorithm objects on :class:`~pyeonclient.Matter`,
not a working directory and not the `eonclient` binary alone.
```

[pyeonclient](https://pypi.org/project/pyeonclient/) exposes the eOn **C++
client algorithms** to Python with nanobind. Geometries are
{py:class}`~pyeonclient.Matter`. Each method is a live class you construct and
call — the same engines the binary uses.

| Algorithm | Python class / function | Binary analogue |
|-----------|-------------------------|-----------------|
| Local min | `Matter.relax()` | `job = minimization` |
| Min-mode | `Dimer`, `ImprovedDimer`, `Lanczos`, `Davidson` | dimer / lanczos / davidson |
| Saddle search | `MinModeSaddleSearch` | `job = saddle_search` |
| NEB | `NudgedElasticBand` | `job = nudged_elastic_band` |
| Hessian | `Hessian` | `job = hessian` |
| Prefactor | `get_prefactors` | `job = prefactor` |
| Opaque batch | `make_job` + `Job.run` | any `JobType` via factory |

`make_job(params)` still constructs every `JobType` the C++ factory supports
(process search, dynamics, basin hopping, TAD, …) for workdir-oriented runs.
Prefer the first-class classes when you want Matter in / Matter out.

## Install

```{code-block} bash
pip install pyeonclient
pip install 'pyeonclient[ase]'
pip install 'rgpot>=2.4.2'   # multi-ABI engines for Metatomic via RGPOT
```

## Shared setup

```{code-block} python
import numpy as np
import pyeonclient as pc

params = pc.Parameters()
params.potential = pc.PotType.LJ          # or RGPOT, METATOMIC, …
params.opt_converged_force = 0.01
params.opt_max_iterations = 1000
pot = pc.make_potential(params.potential, params)

matter = pc.Matter(pot, params)
matter.resize(n)
matter.positions = ...                    # (n, 3) float64
matter.cell = np.eye(3) * L
matter.atomic_numbers = ...
matter.masses = ...
# or: matter = pc.from_ase(atoms, pot, params)
```

## Minimization

```{code-block} python
ok = matter.relax()                       # like LBFGS(atoms).run(fmax=...)
print(ok, matter.potential_energy, matter.max_force)
```

## Min-mode (dimer / Lanczos / Davidson)

```{code-block} python
params.dimer_improved = True              # ImprovedDimer (default)
params.dimer_rotations_max = 10
params.dimer_opt_method = "cg"            # cg | lbfgs | sd

direction = np.random.default_rng(0).normal(size=matter.positions.shape)
dimer = pc.ImprovedDimer(matter, params, pot)
dimer.compute(matter, direction)
print(dimer.eigenvalue)                   # curvature
print(dimer.eigenvector)                  # (n, 3) mode
print(dimer.stats_rotations, dimer.stats_torque)

# Alternatives with the same API:
# pc.Dimer(matter, params, pot)
# pc.Lanczos(matter, params, pot)
# pc.Davidson(matter, params, pot)
```

Parameters: `dimer_improved`, `dimer_rotation_angle`, `dimer_converged_angle`,
`dimer_max_iterations`, `dimer_opt_method`, `dimer_rotations_max` /
`_min`, `dimer_torque_max` / `_min`, `dimer_remove_rotation`.

## Single-ended saddle search

```{code-block} python
params.saddle_minmode_method = "dimer"    # dimer | lanczos | davidson
params.saddle_max_iterations = 1000
params.saddle_max_energy = 20.0

# Displace from the minimum, then search
mode = dimer.eigenvector                  # or a random (n, 3) guess
E_react = matter.potential_energy
ss = pc.MinModeSaddleSearch(matter, mode, E_react, params, pot)
status = ss.run()                         # mutates matter → saddle geometry
print(pc.saddle_status_message(status))
print(ss.eigenvalue, ss.iteration, ss.force_calls)
saddle = matter                           # now at the saddle
```

`SaddleStatus` enumerates the same codes as the C++ client (`GOOD`,
`BAD_MAX_ITERATIONS`, `DIMER_LOST_MODE`, …).

Parameters: `saddle_method`, `saddle_minmode_method`, `saddle_max_iterations`,
`saddle_max_energy`, `saddle_displace_magnitude`, `saddle_displace_type`,
`saddle_converged_force`, `saddle_remove_rotation`.

## NEB

```{code-block} python
path = [pc.from_ase(img, pot, params) for img in images]  # list[Matter]
neb = pc.NudgedElasticBand(path, params, pot)
status = neb.compute()
path = list(neb.path_images())
energies = [neb.image_energy(i) for i in range(neb.n_path)]
```

See {doc}`neb` for energy-weighted springs and OCINEB (`neb_ci_mmf`, …).

## Hessian and prefactors

```{code-block} python
atoms = pc.all_free_atoms(matter)         # int64 indices
hess = pc.Hessian(params, matter)
H = hess.get_hessian(matter, atoms)       # (3n, 3n)
freqs = hess.get_freqs(matter, atoms)

pref1, pref2 = pc.get_prefactors(params, min1, saddle, min2)
```

Parameters: `hessian_atom_list`, `hessian_zero_freq_value`, `prefactor_rate`
(`htst` / `qqhtst`), `prefactor_filter_scheme`, `prefactor_within_radius`, …

## Job factory (full JobType surface)

Every type in `JobType` that the C++ `makeJob` factory implements is available
as a workdir-oriented job:

```{code-block} python
params.job = pc.JobType.Saddle_Search     # or Process_Search, Dynamics, …
job = pc.make_job(params)                 # reads cwd files like eonclient
files = job.run()
pc.write_potcall_summary()
```

Use this for process search, dynamics, basin hopping, TAD, parallel replica,
Monte Carlo, GP surrogate, structure comparison, etc., until those gain
dedicated Matter-first classes. The algorithms above are already first-class.

## ASE bridge

```{code-block} python
matter = pc.from_ase(atoms, pot, params)
atoms = pc.to_ase(matter)
```

Requires `pip install 'pyeonclient[ase]'`.

## Related

- {doc}`dimer`, {doc}`neb`, {doc}`hessian`, {doc}`prefactor`, {doc}`saddle_search`
- {doc}`/devdocs/pyeonclient-pypi` — wheels
- [Cookbook: eOn + PET-MAD NEB](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
