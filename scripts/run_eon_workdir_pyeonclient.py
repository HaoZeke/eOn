#!/usr/bin/env python3
"""Run an eOn workdir by composing pyeonclient ClientEON steps.

Uses :func:`pyeonclient.run_job_in_directory` (load → job → potcalls → timing)
or ``--matter-min`` for the Matter.relax minimization path.

Usage::

    python scripts/run_eon_workdir_pyeonclient.py /path/to/min_reactant
    python scripts/run_eon_workdir_pyeonclient.py /path/to/neb_dir
    python scripts/run_eon_workdir_pyeonclient.py min_reactant --matter-min
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("workdir", type=Path, nargs="?", default=Path("."))
    p.add_argument(
        "--matter-min",
        action="store_true",
        help="Use Matter.relax steps instead of MinimizationJob",
    )
    args = p.parse_args(argv)
    work = args.workdir.resolve()
    if not work.is_dir():
        print(f"error: not a directory: {work}", file=sys.stderr)
        return 2
    if not (work / "config.ini").is_file():
        print(f"error: missing config.ini in {work}", file=sys.stderr)
        return 2
    try:
        import pyeonclient as pc
    except ImportError as e:
        print("error: pyeonclient not importable:", e, file=sys.stderr)
        return 3
    if args.matter_min:
        files = pc.minimize_workdir(work)
    else:
        files = pc.run_job_in_directory(work)
    for f in files:
        print(f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
