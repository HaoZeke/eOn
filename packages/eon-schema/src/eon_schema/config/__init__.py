"""L1 full job-config models.

Still live in ``eon.schema`` (eon-akmc full tree) so conda-forge ``eon`` and
existing consumers are unaffected. Future migration lands here; for now this
module documents the boundary.
"""

from __future__ import annotations


def eon_akmc_schema():
    """Import the full pydantic job config from eon-akmc when installed."""
    try:
        import eon.schema as schema  # type: ignore
    except ImportError as e:
        raise ImportError(
            "Full L1 job config still ships in eon-akmc (import eon.schema). "
            "Install eon-akmc / conda-forge eon for Config models. "
            "eon-schema L2 API specs do not require eon-akmc."
        ) from e
    return schema
