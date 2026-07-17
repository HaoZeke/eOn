#!/usr/bin/env bash
# Publish SSoT from packages/eon-schema (authoring) into monorepo schema/
# so fat release tarballs and tools/params_ssot/codegen.py keep working.
set -euo pipefail
PKG="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$(cd "$PKG/../.." && pwd)"
SRC="$PKG/src/eon_schema/ssot"
DEST="$ROOT/schema"
cp -a "$SRC/eon_params.capnp" "$DEST/"
cp -a "$SRC/eon_params_catalog.json" "$DEST/"
# Keep a pointer README in schema/
cat > "$DEST/README.md" << 'NOTE'
# eOn Cap'n Proto parameter SSoT (tree mirror)

**Authoring home:** `packages/eon-schema/src/eon_schema/ssot/eon_params.capnp`

This directory is a **mirror** for:
- fat release tarballs (`eon-v*.tar.xz`) used by conda-forge `eon-feedstock`
- `tools/params_ssot/codegen.py` historical paths

After editing the package SSoT, run:

```bash
./packages/eon-schema/scripts/sync_ssot_to_tree.sh
python tools/params_ssot/codegen.py
```

Do not add fields only here without updating the package authoring file.
NOTE
echo "synced $SRC → $DEST (fat-tree mirror)"
