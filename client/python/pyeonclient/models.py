"""Chemist-facing specs — re-exported from shared ``eon-schema``.

Install::

    pip install 'pyeonclient[models]'   # pulls eon-schema
    # or: pip install eon-schema

Both pyeonclient and eon-akmc use the same models. Do not re-define specs here.
"""

from __future__ import annotations

try:
    from eon_schema.api import (
        DimerSpec,
        NebSpec,
        parse_dimer_spec,
        parse_neb_spec,
    )
    from eon_schema.fields import Accelerant, MinModeMethod, PathInit
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyeonclient.models requires the shared eon-schema package.\n"
        "  pip install 'pyeonclient[models]'\n"
        "  # or: pip install eon-schema\n"
        "  # monorepo: pip install -e packages/eon-schema\n"
        "Core Dimer/NEB work without models; C++ enforces "
        'accelerant="gp" only with method="improved".'
    ) from e

__all__ = [
    "Accelerant",
    "DimerSpec",
    "MinModeMethod",
    "NebSpec",
    "PathInit",
    "parse_dimer_spec",
    "parse_neb_spec",
]
