#!/usr/bin/env bash
# Swap root pyproject.toml for the pyeonclient distribution identity.
# Usage: scripts/pyeonclient_prepare_pyproject.sh [restore]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
if [[ "${1:-}" == "restore" ]]; then
  if [[ -f pyproject.eon-akmc.toml.bak ]]; then
    mv -f pyproject.eon-akmc.toml.bak pyproject.toml
    echo "restored eon-akmc pyproject.toml"
  fi
  exit 0
fi
if [[ ! -f pyproject-pyeonclient.toml ]]; then
  echo "missing pyproject-pyeonclient.toml" >&2
  exit 1
fi
if [[ ! -f pyproject.eon-akmc.toml.bak ]]; then
  cp -a pyproject.toml pyproject.eon-akmc.toml.bak
fi
cp -a pyproject-pyeonclient.toml pyproject.toml
echo "active pyproject: pyeonclient ($(python3 -c 'import tomllib; print(tomllib.load(open("pyproject.toml","rb"))["project"]["name"], tomllib.load(open("pyproject.toml","rb"))["project"]["version"])'))"
