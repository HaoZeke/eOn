#!/usr/bin/env bash
# DEPRECATED direction: tree was once canonical. Prefer editing package SSoT
# then sync_ssot_to_tree.sh. This script still pulls tree → package if needed
# for one-time recovery.
set -euo pipefail
PKG="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$(cd "$PKG/../.." && pwd)"
cp -a "$ROOT/schema/eon_params.capnp" "$PKG/src/eon_schema/ssot/"
cp -a "$ROOT/schema/eon_params_catalog.json" "$PKG/src/eon_schema/ssot/"
echo "synced tree schema/ → package ssot (recovery only)"
