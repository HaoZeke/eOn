#!/usr/bin/env bash
# Copy monorepo SSoT (schema/) into the eon-schema split package so the
# PyPI wheel ships a self-contained Cap'n Proto + catalog copy.
# Authoring is always schema/eon_params.capnp — edit there, then codegen.
set -euo pipefail
PKG="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$(cd "$PKG/../.." && pwd)"
SRC="$ROOT/schema"
DEST="$PKG/src/eon_schema/ssot"
cp -a "$SRC/eon_params.capnp" "$DEST/"
cp -a "$SRC/eon_params_catalog.json" "$DEST/"
cat > "$DEST/README.md" << 'NOTE'
# Cap'n Proto L0 SSoT (vendored for this package)

**Authoring home:** monorepo `schema/eon_params.capnp`

This directory is a **vendored copy** shipped in the PyPI `eon-schema` wheel
so installs work without the full monorepo. Do not edit field graphs here.

After editing monorepo `schema/`:

```bash
python tools/params_ssot/codegen.py
./packages/eon-schema/scripts/sync_ssot_into_package.sh
```

Fat conda-forge releases use a full monorepo `git archive` tarball (not this package alone).
NOTE
echo "synced $SRC → $DEST (vendored into eon-schema split)"
