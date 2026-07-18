"""Pydantic job-config models for eOn (docs / validation).

**Authoritative implementation:** :mod:`eon_schema.config` (package
``eon-schema``). Both **eon-akmc** and **pyeonclient** consume that package;
this module only re-exports so Sphinx ``eon.schema.*`` and legacy imports stay
stable.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Monorepo: import without a prior editable install of eon-schema.
_pkg_src = Path(__file__).resolve().parents[1] / "packages" / "eon-schema" / "src"
if _pkg_src.is_dir() and str(_pkg_src) not in sys.path:
    sys.path.insert(0, str(_pkg_src))

from eon_schema.config import *  # noqa: F401,F403
from eon_schema.config import __all__ as __all__  # noqa: E402
