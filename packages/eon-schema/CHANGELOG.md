# Changelog

## [Unreleased]

### Added

- `NebSpec` covers CI / energy-weighted springs / OCI-MMF / optimizer knobs and
  `PathInit.file` + `path_list` for class-first cookbook NEB without hand-set
  `Parameters.neb_*` fields.

## [0.2.0] — 2026-07-17

### Added

- INI helpers in ``eon_schema.config.ini``: ``write_ini``, ``read_ini``,
  ``defaults_from_catalog``, ``hydrate_ini``, ``unknown_ini_keys``,
  ``model_to_ini_section`` / ``models_to_ini`` / ``write_models_ini``
  (shared config write path for rgpycrumbs without eon-akmc).

### Changed

- **L1 job-config models live here** (`eon_schema.config`), not only in eon-akmc.
  `eon.schema` and `pyeonclient.models` are re-exports.
- `pydantic>=2` is a **hard** dependency (L1 + L2). Extra `[pydantic]` is a no-op
  kept for older install strings.

### Added

- Full former `eon/schema.py` surface under `eon_schema.config.models`
  (`MainConfig`, `Metatomic`, root `Config`, …).

## [0.1.0] — 2026-07-17

### Added

- Initial PyPI package `eon-schema`.
- Vendored Cap’n Proto L0 copy (`eon_params.capnp` + catalog JSON) from monorepo
  `schema/` (authoring SSoT). Fat `eon-v*` tarball for conda-forge remains the
  full monorepo archive; this package is a PyPI split only.
- Enums: `MinModeMethod`, `Accelerant`, `PathInit`.
- Optional pydantic L2 API: `DimerSpec`, `NebSpec` (extra `[pydantic]`).
- Soft optional dep loading via rgpycrumbs `ensure_import` when available.
- Publishing notes: fat monorepo tarball (conda-forge) + split PyPI packages.
