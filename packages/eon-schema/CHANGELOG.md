# Changelog

## [0.1.0] — 2026-07-17

### Added

- Initial PyPI package `eon-schema`.
- Vendored Cap’n Proto L0 copy (`eon_params.capnp` + catalog JSON) from monorepo
  `schema/` (authoring remains in the full eOn tree).
- Enums: `MinModeMethod`, `Accelerant`, `PathInit`.
- Optional pydantic L2 API: `DimerSpec`, `NebSpec` (extra `[pydantic]`).
- Soft optional dep loading via rgpycrumbs `ensure_import` when available.
- Publishing / release notes for feedstock-safe dual shipping.
