# RGPot + CPMD BLYP (cpmdc) point energy example

Small water molecule point calculation through eOn's **RGPot** potential, using
rgpot's **CPMDPot** frontend and the **libcpmdc** embedded OpenCPMD engine.

## Requirements

- eOn built with `-Dwith_rgpot=true` (Meson) or `-DWITH_RGPOT=ON` (CMake)
- `libcpmdc.so` on the loader path, or set one of:
  - `CPMDC_LIBRARY=/path/to/libcpmdc.so`
  - `RGPOT_CPMDC_ENGINE=/path/to/libcpmdc.so`
  - `[RGPot] engine_path` in `config.ini`

## Run

```bash
export CPMDC_LIBRARY=/path/to/libcpmdc.so
# eonclient from the build directory, for example:
/path/to/bbdir/client/eonclient
```

Expect finite energy in `results.dat` / client logs and a successful Point job.
