#!/usr/bin/env python3
"""Compare metatomic force engines on one Matter geometry.

Backends (key → factory via pyeonclient.backends)::

  metatomic / metatomic_fat   — fat-linked PotType.METATOMIC
  rgpot_metatomic / metatomic_dlopen — RGPOT dlopen libmetatomic_engine.so
  ase                         — metatomic_ase.MetatomicCalculator via ASE wrap

Env::

  EON_PET_MAD_MODEL   path to .pt (required)
  RGPOT_METATOMIC_ENGINE  path to libmetatomic_engine.so (for dlopen)
  EON_PET_MAD_POS     optional pos.con (else cookbook default / d016)

Usage::

  python scripts/compare_metatomic_backends.py
  python scripts/compare_metatomic_backends.py --json /tmp/mta-compare.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


def _load_geometry(pos_con: Path):
    """Return (positions, Z int64, cell, masses) from a .con via readcon or ASE."""
    try:
        from ase.io import read

        atoms = read(str(pos_con))
        pos = np.ascontiguousarray(atoms.get_positions(), dtype=np.float64)
        z = np.ascontiguousarray(atoms.get_atomic_numbers(), dtype=np.int64)
        cell = np.ascontiguousarray(np.asarray(atoms.get_cell()), dtype=np.float64)
        masses = np.ascontiguousarray(atoms.get_masses(), dtype=np.float64)
        return pos, z, cell, masses, bool(np.any(atoms.pbc))
    except Exception as e:
        raise SystemExit(f"failed to read {pos_con}: {e}") from e


def _eval(name: str, pot, pos, z, cell, n_repeat: int = 3) -> dict:
    out: dict = {"backend": name, "ok": False}
    try:
        # warm-up
        E0, F0 = pot.get_ef(pos, z, cell)
        t0 = time.perf_counter()
        E, F = E0, F0
        for _ in range(n_repeat):
            E, F = pot.get_ef(pos, z, cell)
        elapsed = time.perf_counter() - t0
        F = np.asarray(F, dtype=np.float64)
        out.update(
            {
                "ok": True,
                "energy": float(E),
                "max_force": float(np.max(np.linalg.norm(F, axis=1))),
                "force_rms": float(np.sqrt(np.mean(F * F))),
                "n_atoms": int(pos.shape[0]),
                "seconds_per_call": elapsed / n_repeat,
                "n_repeat": n_repeat,
            }
        )
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--model",
        default=os.environ.get("EON_PET_MAD_MODEL", ""),
        help="metatomic .pt model path",
    )
    ap.add_argument(
        "--engine",
        default=os.environ.get("RGPOT_METATOMIC_ENGINE", ""),
        help="libmetatomic_engine.so for dlopen path",
    )
    ap.add_argument(
        "--pos",
        default=os.environ.get("EON_PET_MAD_POS", ""),
        help="pos.con geometry",
    )
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--json", dest="json_out", default="")
    args = ap.parse_args(argv)

    model = Path(args.model).expanduser() if args.model else None
    if not model or not model.is_file():
        # default terra paths
        for cand in (
            Path(
                "/home/rgoswami/Git/Github/TheochemUI/eOn-pyeon-final/"
                "subprojects/gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt"
            ),
            Path(
                "/home/rgoswami/Git/Github/TheochemUI/gpr_optim/"
                "bench_data/petmad/pet-mad-s-v1.5.0.pt"
            ),
        ):
            if cand.is_file():
                model = cand
                break
    if not model or not model.is_file():
        print("error: set --model / EON_PET_MAD_MODEL", file=sys.stderr)
        return 2

    pos_con = Path(args.pos) if args.pos else None
    if not pos_con or not pos_con.is_file():
        for cand in (
            Path(
                "/home/rgoswami/Git/Github/epfl/pixi_envs/atomistic-cookbook/"
                "atomistic-cookbook/examples/eon-pet-neb/min_reactant/pos.con"
            ),
            model.parent / "d016_pos.con",
            model.parent / "minimal_pt.con",
        ):
            if cand.is_file():
                pos_con = cand
                break
    if not pos_con or not pos_con.is_file():
        print("error: set --pos / EON_PET_MAD_POS", file=sys.stderr)
        return 2

    engine = Path(args.engine) if args.engine else None
    if not engine or not engine.is_file():
        for cand in (
            Path("build-pyeon-mta-rgpot/client/libmetatomic_engine.so"),
            Path("build-pyeon-mta/client/libmetatomic_engine.so"),
            Path(
                "/home/rgoswami/Git/Github/TheochemUI/eOn-pyeon-final/"
                "build-pyeon-mta-rgpot/client/libmetatomic_engine.so"
            ),
            Path(
                "/home/rgoswami/Git/Github/TheochemUI/eOn-pyeon-final/"
                "build-pyeon-mta/client/libmetatomic_engine.so"
            ),
        ):
            if cand.is_file():
                engine = cand
                break

    import pyeonclient as pyec
    from pyeonclient.backends import list_backends, make_backend

    print("pyeonclient", pyec.__version__)
    print("built_with_metatomic", pyec.built_with_metatomic())
    print("built_with_rgpot", pyec.built_with_rgpot())
    print("backends", list_backends())
    print("model", model)
    print("pos", pos_con)
    print("engine", engine)

    pos, z, cell, masses, pbc = _load_geometry(pos_con)
    print(f"geometry n={len(z)} pbc={pbc}")

    jobs: list[tuple[str, dict]] = []
    # fat
    jobs.append(
        (
            "metatomic_fat",
            {"model_path": str(model), "device": args.device},
        )
    )
    # dlopen
    jobs.append(
        (
            "rgpot_metatomic",
            {
                "model_path": str(model),
                "engine_path": str(engine) if engine else "",
                "device": args.device,
            },
        )
    )
    # ase metatomic
    jobs.append(("ase_metatomic", {"model_path": str(model), "device": args.device}))

    results = []
    for key, kw in jobs:
        print(f"\n=== {key} ===")
        try:
            if key == "ase_metatomic":
                from metatomic.torch import load_atomistic_model
                from metatomic_ase import MetatomicCalculator

                calc = MetatomicCalculator(
                    load_atomistic_model(str(model)), device=args.device
                )
                pot = make_backend("ase", calculator=calc)
            else:
                pot = make_backend(key if key != "metatomic_fat" else "metatomic", **kw)
            row = _eval(key, pot, pos, z, cell, n_repeat=args.repeat)
        except Exception as e:
            row = {"backend": key, "ok": False, "error": f"{type(e).__name__}: {e}"}
        results.append(row)
        if row.get("ok"):
            print(
                f"  E={row['energy']:.8f}  "
                f"max|F|={row['max_force']:.6e}  "
                f"rmsF={row['force_rms']:.6e}  "
                f"t={row['seconds_per_call']*1e3:.2f} ms/call"
            )
        else:
            print(f"  FAIL {row.get('error')}")

    # pairwise energy deltas among successes
    ok = [r for r in results if r.get("ok")]
    if len(ok) >= 2:
        print("\n=== energy deltas (row - col) ===")
        names = [r["backend"] for r in ok]
        for i, a in enumerate(ok):
            for j, b in enumerate(ok):
                if j <= i:
                    continue
                dE = a["energy"] - b["energy"]
                dF = abs(a["max_force"] - b["max_force"])
                print(
                    f"  {a['backend']} - {b['backend']}: "
                    f"dE={dE:.6e}  d(max|F|)={dF:.6e}"
                )

    payload = {
        "model": str(model),
        "pos": str(pos_con),
        "engine": str(engine) if engine else None,
        "device": args.device,
        "built_with_metatomic": bool(pyec.built_with_metatomic()),
        "built_with_rgpot": bool(pyec.built_with_rgpot()),
        "results": results,
    }
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(payload, indent=2) + "\n")
        print("wrote", args.json_out)
    return 0 if any(r.get("ok") for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
