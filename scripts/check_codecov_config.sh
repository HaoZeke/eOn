#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
YML="$ROOT/codecov.yml"
WF="$ROOT/.github/workflows/coverage.yml"
fail=0
die() { echo "ERROR: $*" >&2; fail=1; }
ok() { echo "OK: $*"; }
[[ -f "$YML" ]] || die "missing codecov.yml"
[[ -f "$WF" ]] || die "missing coverage.yml"
for flag in python cpp fortran; do
  grep -qE "name:[[:space:]]*${flag}" "$YML" && ok "codecov.yml $flag" || die "missing flag $flag"
  if grep -E "flags:[[:space:]]*${flag}" "$WF" | grep -vq '^\s*#'; then
    ok "coverage.yml flags $flag"
  else
    die "coverage.yml missing flags $flag"
  fi
done
grep -q 'use_oidc: true' "$WF" || die "missing use_oidc"
grep -q 'id-token: write' "$WF" || die "missing id-token"
grep -q 'fail_ci_if_error: false' "$WF" || die "missing soft-fail"
grep -q 'app.codecov.io' "$WF" || die "missing dashboard note"
for s in run_coverage_python.sh run_coverage_cpp.sh run_coverage_fortran.sh; do
  [[ -x "$ROOT/scripts/$s" ]] || die "missing $s"
  ok "scripts/$s"
done
[[ "$fail" -eq 0 ]] || exit 1
echo "check_codecov_config: all checks passed"
