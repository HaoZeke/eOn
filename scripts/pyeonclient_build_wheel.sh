#!/usr/bin/env bash
# Build a pyeonclient wheel for the *current* interpreter (PyPI-style).
#
# Variants (env):
#   PYEONCLIENT_VARIANT=base        # default: with_rgpot, no torch link
#   PYEONCLIENT_VARIANT=metatomic   # with_metatomic + with_rgpot (needs torch
#                                   # + metatomic-torch on the build env)
#
# Examples (same idea as building against a local torch install):
#   pip install -U build nanobind numpy meson ninja meson-python
#   pip install 'torch' 'metatomic-torch' 'metatensor-torch' 'vesin'  # CPU ok
#   PYEONCLIENT_VARIANT=metatomic ./scripts/pyeonclient_build_wheel.sh
#
# Override any meson flag:
#   PYEONCLIENT_SETUP_ARGS='-Dwith_metatomic=true -Dtorch_path=$CONDA_PREFIX' \
#     ./scripts/pyeonclient_build_wheel.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
VARIANT="${PYEONCLIENT_VARIANT:-base}"

bash scripts/pyeonclient_prepare_pyproject.sh

cleanup() { bash scripts/pyeonclient_prepare_pyproject.sh restore || true; }
trap cleanup EXIT

python -m pip install -U pip build nanobind 'numpy>=1.26.4' meson ninja meson-python

SETUP_ARGS=(
  "-Dwith_pyeonclient=true"
  "-Dinstall_eon_server=false"
  "-Dwith_tests=false"
  "-Dwith_ase=false"
  "-Dwith_xtb=false"
  "-Dwith_serve=false"
  "-Dwith_fortran=true"
  "-Dwith_cuh2=true"
  "-Dwith_rgpot=true"
  "-Dwrap_mode=forcefallback"
)

case "$VARIANT" in
  base)
    SETUP_ARGS+=("-Dwith_metatomic=false")
    ;;
  metatomic)
    SETUP_ARGS+=("-Dwith_metatomic=true" "-Dpip_metatomic=true")
    # Prefer an existing torch install for headers/libs
    if [[ -n "${TORCH_PATH:-}" ]]; then
      SETUP_ARGS+=("-Dtorch_path=${TORCH_PATH}")
      SETUP_ARGS+=("-Dpip_metatomic=false")
    fi
    python -c "import torch; print('torch', torch.__version__)" || {
      echo "metatomic variant needs torch on PYTHONPATH (pip install torch)" >&2
      exit 2
    }
    ;;
  *)
    echo "unknown PYEONCLIENT_VARIANT=$VARIANT (base|metatomic)" >&2
    exit 2
    ;;
esac

if [[ -n "${PYEONCLIENT_SETUP_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  SETUP_ARGS+=($PYEONCLIENT_SETUP_ARGS)
fi

CFG=()
for a in "${SETUP_ARGS[@]}"; do
  CFG+=(--config-setting="setup-args=${a}")
done

echo "building pyeonclient variant=${VARIANT} args=${SETUP_ARGS[*]}"
python -m build --wheel "${CFG[@]}"
ls -la dist/
python - <<'PY'
import glob, zipfile, sys
whls = sorted(glob.glob("dist/*.whl"))
if not whls:
    sys.exit("no wheel produced")
print("wheel", whls[-1])
with zipfile.ZipFile(whls[-1]) as z:
    names = [n for n in z.namelist() if "pyeonclient" in n][:20]
    print("sample members", names)
PY
