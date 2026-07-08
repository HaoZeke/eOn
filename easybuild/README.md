# EasyBuild packaging for eOn

Draft easyconfigs for contributing eOn to
[easybuild-easyconfigs](https://github.com/easybuilders/easybuild-easyconfigs).

| File | Role |
|------|------|
| `easyconfigs/q/quill/quill-11.1.0-GCCcore-13.3.0.eb` | C++ logging dependency (missing from upstream EB) |
| `easyconfigs/m/Meson/Meson-1.8.2-GCCcore-13.3.0.eb` | Meson ≥1.8 for eOn `meson_version` floor (2024a only has 1.4.0) |
| `easyconfigs/m/Meson/Meson-1.8.2_reenable-binutils-workaround.patch` | Upstream EB patch for Meson 1.8.2 |
| `easyconfigs/e/eOn/eOn-2.16.0-gfbf-2024a.eb` | Core optimized eOn client + Python AKMC package |

**Toolchain:** `gfbf/2024a` (GCC 13.3 + FlexiBLAS + FFTW; no MPI — core AKMC does not need OpenMPI).

**Core profile:** release build, Fortran pots, CuH2 on; optional engines
(xTB, ARTn, metatomic, ASE-embedded, VASP/AMS, serve/RPC, parallel NEB)
off. Full flag table and offline/EESSI notes live in the internal packaging
doc (`Software/eOn/eon-easybuild-v2.16.0-packaging.org`).

## Pre-flight with eb-stack (recommended)

Before `eb --robot`, validate parse + robot coverage + packaging gates:

```bash
eb-stack check-recipe \
  --recipe easyconfigs/e/eOn/eOn-2.16.0-gfbf-2024a.eb \
  --easyconfigs /path/to/easybuild-easyconfigs \
  --easyconfigs easyconfigs \
  --require-configopt=-Dwith_fortran=true \
  --require-configopt=-Dwith_tests=false
```

Expect all runtime/build deps found (including `quill` on `GCCcore-13.3.0`)
and packaging gate OK (easyblock, moduleclass, checksums, configopts).

## Build (developer)

```bash
export EASYBUILD_ROBOT_PATHS=$PWD/easyconfigs:$EASYBUILD_ROBOT_PATHS
eb quill-11.1.0-GCCcore-13.3.0.eb --robot
eb Meson-1.8.2-GCCcore-13.3.0.eb --robot
# Release builds enable -O3 -flto=auto in client/meson.build; cap parallel
# on memory-constrained builders:
eb eOn-2.16.0-gfbf-2024a.eb --robot --parallel=4
```

Sanity (after module load):

```bash
eonclient --help
python -c 'import eon'
```

## Notes for reviewers

- Companion **quill** recipe is required; not yet in upstream easyconfigs.
- `configopts` is a single plain string (not `' '.join`) so tooling can gate flags.
- `nlohmann_json` and `readcon-core` still come from Meson wraps (`Rust` builddep).
- Portable: `native_arch=false`, `fast_math=false`, `use_mkl=false`.