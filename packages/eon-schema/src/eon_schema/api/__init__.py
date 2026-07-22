"""L2 API composition models (in-process constructors)."""

from __future__ import annotations

from typing import Any

from eon_schema.api.dimer import parse_dimer_spec
from eon_schema.api.neb import parse_neb_spec


def __getattr__(name: str) -> Any:
    if name == "DimerSpec":
        from eon_schema.api.dimer import _dimer_spec_cls

        return _dimer_spec_cls()
    if name == "NebSpec":
        from eon_schema.api.neb import _neb_spec_cls

        return _neb_spec_cls()
    raise AttributeError(name)


__all__ = [
    "DimerSpec",
    "NebSpec",
    "parse_dimer_spec",
    "parse_neb_spec",
]
