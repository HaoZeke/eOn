# EasyBuild packaging for eOn

GROMACS / LAMMPS / Siesta-style easyconfigs for contributing eOn to
[easybuild-easyconfigs](https://github.com/easybuilders/easybuild-easyconfigs).

## Primary product (conda-forge parity)

| File | Role |
|------|------|
| `easyconfigs/e/eOn/eOn-2.16.0-foss-2024a.eb` | **Full default**: metatomic + xTB + serve + rgpot + Fortran/CuH2 |
| `easyconfigs/m/metatensor/metatensor-0.1.17-GCCcore-13.3.0.eb` | Companion C library |
| `easyconfigs/m/metatensor-torch/metatensor-torch-0.10.0-foss-2024a.eb` | Companion TorchScript bindings |
| `easyconfigs/m/metatomic-torch/metatomic-torch-0.1.15-foss-2024a.eb` | Companion metatomic C++ API |
| `easyconfigs/q/quill/quill-11.1.0-GCCcore-13.3.0.eb` | C++ logging |
| `easyconfigs/m/Meson/Meson-1.8.2-GCCcore-13.3.0.eb` | Meson ≥1.8 (2024a only ships 1.4.0) |

Optional lite leftover (not for submission): `eOn-2.16.0-gfbf-2024a-core.eb`.

**Toolchain:** `foss/2024a` (matches PyTorch + GROMACS/LAMMPS). Feature flags match
conda-forge `eon-feedstock` `build.sh`. Multi-source stages `rgpot` and
`readcon-core` into `subprojects/` (Siesta/QE style).

CUDA torch is a **versionsuffix** sibling when needed (like GROMACS `-CUDA-…`),
not a dozen option-off stubs.

## Prep with eb-stack (required workflow)

```bash
eb-stack check-recipe \
  --recipe easyconfigs/e/eOn/eOn-2.16.0-foss-2024a.eb \
  --easyconfigs /path/to/easybuild-easyconfigs \
  --easyconfigs easyconfigs \
  --require-configopt=-Dwith_metatomic=true \
  --require-configopt=-Dwith_xtb=true \
  --require-configopt=-Dwith_serve=true \
  --require-configopt=-Dwith_rgpot=true \
  --require-configopt=-Dpip_metatomic=false

# If deps are missing, scaffold drafts into the overlay and re-run:
eb-stack check-recipe \
  --recipe easyconfigs/e/eOn/eOn-2.16.0-foss-2024a.eb \
  --easyconfigs /path/to/easybuild-easyconfigs \
  --easyconfigs easyconfigs \
  --scaffold-missing easyconfigs
```

## Build

```bash
export EASYBUILD_ROBOT_PATHS=$PWD/easyconfigs:$EASYBUILD_ROBOT_PATHS
eb metatensor-0.1.17-GCCcore-13.3.0.eb --robot
eb metatensor-torch-0.10.0-foss-2024a.eb --robot
eb metatomic-torch-0.1.15-foss-2024a.eb --robot
eb quill-11.1.0-GCCcore-13.3.0.eb --robot
eb Meson-1.8.2-GCCcore-13.3.0.eb --robot
# LTO + fat deps: cap parallel on tight nodes
eb eOn-2.16.0-foss-2024a.eb --robot --parallel=4
```

Sanity:

```bash
eonclient --help
eonclient --version
ldd $(which eonclient) | grep -Ei 'libxtb|metatomic|capnp|torch'
python -c 'import eon'
```
