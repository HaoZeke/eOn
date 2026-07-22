"""Optional dependency loading (prefer rgpycrumbs; no hard dep)."""

from __future__ import annotations

import importlib
import os
from typing import Any, Optional


def ensure_import(
    module: str,
    *,
    install_spec: Optional[str] = None,
    extra_hint: Optional[str] = None,
) -> Any:
    """Import *module*, optionally auto-installing via rgpycrumbs.

    Resolution order:
    1. ``rgpycrumbs.api.ensure_import`` / ``rgpycrumbs._aux.ensure_import``
       when rgpycrumbs is installed (honors ``RGPYCRUMBS_AUTO_DEPS=1`` /
       ``RGPKGS_AUTO_DEPS``).
    2. Plain ``importlib.import_module``.
    3. ``ImportError`` with pip/uv install instructions.

    Does **not** hard-require ``uv`` or ``rgpycrumbs`` on the base wheel.
    """
    spec = install_spec or module
    try:
        try:
            from rgpycrumbs.api import ensure_import as _ei  # type: ignore
        except ImportError:
            from rgpycrumbs._aux import ensure_import as _ei  # type: ignore
        # Best-effort call — rgpycrumbs signatures vary slightly by version.
        try:
            return _ei(module, install_spec=spec)
        except TypeError:
            try:
                return _ei(module, spec)
            except TypeError:
                return _ei(module)
    except ImportError:
        pass

    try:
        return importlib.import_module(module)
    except ImportError as e:
        auto = os.environ.get("RGPYCRUMBS_AUTO_DEPS") or os.environ.get(
            "RGPKGS_AUTO_DEPS"
        )
        hints = [
            f"Missing optional dependency {module!r}.",
            f"  pip install {spec!r}",
            f"  uv pip install {spec!r}",
            "  # or: pip install 'eon-schema[pydantic]'",
        ]
        if not auto:
            hints.append(
                "Opt-in auto-install (requires rgpycrumbs + uv on PATH): "
                "export RGPYCRUMBS_AUTO_DEPS=1"
            )
        if extra_hint:
            hints.append(extra_hint)
        raise ImportError("\n".join(hints)) from e
