# Copyright (c) eOn Development Team
# SPDX-License-Identifier: BSD-3-Clause
"""Compare a filtered prefactor with the all-free-atom prefactor.

The client under test is ``--eonclient``. A geometry directory supplies
``reactant.con``, ``saddle.con``, ``product.con``, and the potential
files. The unfixed client ignores ``all_free_atoms``, so the two
prefactors match and this command exits 1. The fixed client multiplies
every free atom when the flag is true, so the prefactors differ and
this command exits 0.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path


def clock() -> str:
    zone = timezone(timedelta(hours=-6))
    return datetime.now(zone).strftime("%Y-%m-%dT%H:%M:%S%z")


def write_config(dest: Path, all_free: bool) -> None:
    flag = "true" if all_free else "false"
    dest.joinpath("config.ini").write_text(
        "\n".join(
            [
                "[Main]",
                "job = prefactor",
                "temperature = 160.0",
                "finite_difference = 0.01",
                "",
                "[Potential]",
                "potential = lammps",
                "lammps_threads = 1",
                "",
                "[Communicator]",
                "type = local",
                "num_jobs = 1",
                "number_of_CPUs = 1",
                "",
                "[Prefactor]",
                "default_value = 0.0",
                f"all_free_atoms = {flag}",
                "filter_scheme = fraction",
                "filter_fraction = 0.90",
                "",
            ]
        )
    )


def copy_geometry(src: Path, dest: Path) -> None:
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(
        src,
        dest,
        ignore=shutil.ignore_patterns(
            "client.log",
            "client_quill.log",
            "results.dat",
            "freq.dat",
            "freqs.dat",
            "hessian.dat",
            "config.ini",
        ),
    )


def prefactor(text: str) -> tuple[float, float]:
    forward = reverse = None
    for line in text.splitlines():
        if line.endswith("prefactor_reactant_to_product"):
            forward = float(line.split()[0])
        elif line.endswith("prefactor_product_to_reactant"):
            reverse = float(line.split()[0])
    if forward is None or reverse is None:
        raise RuntimeError("results.dat has no prefactor")
    return forward, reverse


def atom_line(log: str) -> str:
    for line in log.splitlines():
        if "including" in line and "atoms in the hessian" in line:
            return line.strip()
    return ""


def run_one(client: Path, src: Path, dest: Path, all_free: bool) -> dict:
    copy_geometry(src, dest)
    write_config(dest, all_free)
    env = os.environ.copy()
    proc = subprocess.run(
        [str(client)],
        cwd=dest,
        env=env,
        capture_output=True,
        text=True,
    )
    log_path = dest / "client.log"
    log = log_path.read_text(errors="replace") if log_path.is_file() else proc.stdout + proc.stderr
    results = dest / "results.dat"
    if proc.returncode != 0 or not results.is_file():
        sys.stderr.write(log)
        raise RuntimeError(f"eonclient exited {proc.returncode}")
    forward, reverse = prefactor(results.read_text())
    return {
        "forward": forward,
        "reverse": reverse,
        "atoms": atom_line(log),
        "exit": proc.returncode,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eonclient", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    print(f"start_time {clock()}")
    print("command " + " ".join(Path(part).name if part.startswith("/") else part for part in sys.argv))
    client = args.eonclient.resolve()
    version = subprocess.check_output([str(client), "--version"], text=True, stderr=subprocess.STDOUT).strip()
    print("eonclient_version_begin")
    print(version)
    print("eonclient_version_end")
    filtered = run_one(client, args.geometry, args.work / "filtered", False)
    full = run_one(client, args.geometry, args.work / "all_free", True)
    print(f"filtered_forward {filtered['forward']:.12e}")
    print(f"filtered_reverse {filtered['reverse']:.12e}")
    print(f"filtered_atoms {filtered['atoms']}")
    print(f"all_free_forward {full['forward']:.12e}")
    print(f"all_free_reverse {full['reverse']:.12e}")
    print(f"all_free_atoms {full['atoms']}")
    same = filtered["forward"] == full["forward"] and filtered["reverse"] == full["reverse"]
    print(f"prefactors_match {str(same).lower()}")
    print(f"end_time {clock()}")
    return 1 if same else 0


if __name__ == "__main__":
    raise SystemExit(main())
