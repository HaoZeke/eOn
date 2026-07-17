# pyeonclient

In-process eOn client for Python. First-class algorithm objects on **Matter**
(not a workdir wrapper).

| Class / function | Role |
|------------------|------|
| `Matter` | Geometry + PEF (ASE `Atoms` analogue) |
| `ImprovedDimer` / `Dimer` / `Lanczos` / `Davidson` | Min-mode |
| `MinModeSaddleSearch` | Single-ended saddle search |
| `NudgedElasticBand` | NEB band (`list[Matter]`) |
| `Hessian` / `get_prefactors` | Vibrational analysis / HTST |
| `make_job` + `Job.run` | Full `JobType` factory (workdir jobs) |

Docs: [Python API](https://eondocs.org/user_guide/pyeonclient.html)

## Install

```bash
pip install pyeonclient 'pyeonclient[ase]' 'rgpot>=2.4.2'
```

## Example (dimer → saddle)

```python
import numpy as np
import pyeonclient as pc

params = pc.Parameters()
params.potential = pc.PotType.LJ
pot = pc.make_potential(params.potential, params)
matter = pc.Matter(pot, params)
# ... set positions / cell / numbers ...

matter.relax()
mode0 = np.random.default_rng(0).normal(size=matter.positions.shape)
dimer = pc.ImprovedDimer(matter, params, pot)
dimer.compute(matter, mode0)

ss = pc.MinModeSaddleSearch(
    matter, dimer.eigenvector, matter.potential_energy, params, pot
)
status = ss.run()
print(pc.saddle_status_message(status), ss.eigenvalue)
```

## NEB

```python
path = [pc.from_ase(img, pot, params) for img in images]
neb = pc.NudgedElasticBand(path, params, pot)
neb.compute()
path = list(neb.path_images())
```

## Build

```bash
meson setup build -Dwith_pyeonclient=true
meson compile -C build
pytest tests/test_pyeonclient_eigenmode.py tests/test_pyeonclient_neb.py -v
```
