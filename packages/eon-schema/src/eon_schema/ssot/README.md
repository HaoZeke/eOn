# Cap’n Proto L0 SSoT (authoring)

Edit **`eon_params.capnp`** here first.

Then from the monorepo root:

```bash
./packages/eon-schema/scripts/sync_ssot_to_tree.sh
python tools/params_ssot/codegen.py
```

That updates the fat-tree mirror `schema/` (conda-forge tarball) and generated
client/Python catalogs.
