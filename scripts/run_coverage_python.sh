#!/usr/bin/env bash
# pytest-cov over the pure-Python eon AKMC package (import eon).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT_XML="${1:-python_coverage.xml}"

python3 -m pip install -U pip -q 2>/dev/null || true
python3 -m pip install -q pytest pytest-cov coverage numpy pyyaml 'readcon>=0.8' 2>/dev/null || true

# Prefer installed package; fall back to source tree on PYTHONPATH.
export PYTHONPATH="${ROOT}${PYTHONPATH:+:$PYTHONPATH}"

echo "==> pytest eon package coverage"
# Lightweight unit suite (no full eonclient AKMC runs).
# test_one_pt / test_akmc_pt need the client binary and are optional.
python3 -m pytest \
  tests/test_config_metadata.py \
  tests/test_displacement_atom_list.py \
  -q --tb=short \
  --cov=eon \
  --cov-report=xml:"$ROOT/$OUT_XML" \
  --cov-report=term-missing \
  --cov-branch

test -s "$ROOT/$OUT_XML"
python3 - "$ROOT/$OUT_XML" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
text = p.read_text(encoding="utf-8", errors="replace")
print(f"python coverage XML size={p.stat().st_size}")
if "filename=" not in text:
    raise SystemExit("empty coverage XML")
# rough line summary if present
print("OK python coverage artifact")
PY
echo "OK wrote $ROOT/$OUT_XML"
