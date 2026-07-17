#!/usr/bin/env bash
# Alias for sync_ssot_into_package.sh (schema/ → package vendored copy).
exec "$(cd "$(dirname "$0")" && pwd)/sync_ssot_into_package.sh"
