# Changelog

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
