"""Vendored Cap'n Proto L0 SSoT (copy of monorepo schema/)."""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path
from typing import Any


def _ssot_dir() -> Path:
    # importlib.resources for installed wheel; __file__ for editable
    try:
        root = resources.files("eon_schema.ssot")
        # Traversable → Path when possible
        p = Path(str(root))
        if (p / "eon_params.capnp").is_file():
            return p
    except Exception:
        pass
    return Path(__file__).resolve().parent


def capnp_path() -> Path:
    """Path to vendored ``eon_params.capnp``."""
    return _ssot_dir() / "eon_params.capnp"


def catalog_path() -> Path:
    """Path to vendored ``eon_params_catalog.json``."""
    return _ssot_dir() / "eon_params_catalog.json"


def load_catalog() -> dict[str, Any]:
    """Load the generated L0 field catalog (JSON)."""
    with catalog_path().open(encoding="utf-8") as f:
        return json.load(f)


def covered_sections() -> list[str]:
    cat = load_catalog()
    return list(cat.get("sections", {}).keys())
