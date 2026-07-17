"""eon-schema — shared eOn parameter SSoT + job config + API models.

Layers
------
* **L0** — Cap'n Proto field graph (vendored copy of monorepo ``schema/``).
* **L1** — full job config pydantic models (``eon_schema.config``); also
  re-exported as ``eon.schema`` from eon-akmc for Sphinx / legacy imports.
* **L2** — in-process API composition (``DimerSpec``, ``NebSpec``).

Import::

    import eon_schema
    from eon_schema.ssot import capnp_path, load_catalog
    from eon_schema.config import MainConfig, Config  # needs pydantic
    from eon_schema.api import DimerSpec, NebSpec
"""

from __future__ import annotations

__version__ = "0.2.0"

from eon_schema.ssot import capnp_path, catalog_path, load_catalog

__all__ = [
    "__version__",
    "capnp_path",
    "catalog_path",
    "load_catalog",
]
