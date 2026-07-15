#!/usr/bin/env python3
"""Run an eOn workdir (config.ini + .con files) via pyeonclient.

Same I/O contract as ``eonclient``: chdir to *workdir*, load config.ini,
run the configured job, write results.dat and job-specific outputs.

Usage::

    python scripts/run_eon_workdir_pyeonclient.py /path/to/min_reactant
    python scripts/run_eon_workdir_pyeonclient.py .   # cwd

Requires a built ``pyeonclient`` with the potentials used in config.ini
(e.g. ``-Dwith_metatomic=true`` for Metatomic jobs).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "workdir",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Directory with config.ini (and pos.con / reactant.con / …)",
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
        print(
            "error: pyeonclient not importable (build with "
            "-Dwith_pyeonclient=true):",
            e,
            file=sys.stderr,
        )
        return 3
    files = pc.run_job_in_directory(str(work), pc.Parameters())
    for f in files:
        print(f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
