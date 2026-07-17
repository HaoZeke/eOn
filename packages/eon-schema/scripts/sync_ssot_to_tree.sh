#!/usr/bin/env bash
# Deprecated name. Authoring is monorepo schema/; package is vendored.
# This used to push package → schema/. Use sync_ssot_into_package.sh instead.
set -euo pipefail
echo "error: authoring is schema/eon_params.capnp (monorepo), not the package." >&2
echo "  1. edit schema/eon_params.capnp" >&2
echo "  2. python tools/params_ssot/codegen.py" >&2
echo "  3. ./packages/eon-schema/scripts/sync_ssot_into_package.sh" >&2
exit 1
