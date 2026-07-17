#!/usr/bin/env python3
"""Ensure pyeonclient version is lockstep across pyproject + C++ module."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib  # type: ignore

    meta = tomllib.loads((ROOT / "pyproject-pyeonclient.toml").read_text())
    v_py = meta["project"]["version"]
    mod = (ROOT / "client/python/bind/module.cpp").read_text()
    m = re.search(r'__version__"\)\s*=\s*"([^"]+)"', mod)
    if not m:
        print("ERROR: no __version__ in module.cpp", file=sys.stderr)
        return 1
    v_cpp = m.group(1)
    if v_py != v_cpp:
        print(f"ERROR: version mismatch pyproject={v_py} module.cpp={v_cpp}", file=sys.stderr)
        return 1
    print(f"pyeonclient version OK: {v_py}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
