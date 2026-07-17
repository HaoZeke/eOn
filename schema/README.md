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
