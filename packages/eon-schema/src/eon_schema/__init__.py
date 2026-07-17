"""eon-schema — shared eOn parameter SSoT copy + API models.

Layers
------
* **L0** — Cap'n Proto field graph (vendored copy of monorepo ``schema/``).
* **L2** — in-process API composition (``DimerSpec``, ``NebSpec``).
* **L1** — full job config pydantic models remain in ``eon.schema`` (eon-akmc)
  until migrated; this package does not replace the full-tree eOn release.

Import::

    import eon_schema
    from eon_schema.api import DimerSpec, NebSpec  # needs pydantic
    from eon_schema.ssot import capnp_path, catalog_path
"""

from __future__ import annotations

__version__ = "0.1.0"

from eon_schema.ssot import capnp_path, catalog_path, load_catalog

__all__ = [
    "__version__",
    "capnp_path",
    "catalog_path",
    "load_catalog",
]
