#!/usr/bin/env python3
"""Plot fat vs RGPOT-dlopen vs ASE-metatomic single-point cost for eOn docs.

Reads the JSON produced by ``scripts/compare_metatomic_backends.py`` (or the
committed copy under ``docs/source/fig/data/``) and writes an SVG bar chart.

Usage::

  python scripts/plot_metatomic_backend_bench.py
  python scripts/plot_metatomic_backend_bench.py \\
      --json docs/source/fig/data/metatomic_backend_bench.json \\
      --out docs/source/fig/metatomic_backend_bench.svg
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# Display order and labels for the docs figure.
_LABELS = {
    "metatomic": "fat",
    "metatomic_fat": "fat",
    "rgpot_metatomic": "RGPOT dlopen",
    "metatomic_dlopen": "RGPOT dlopen",
    "ase_metatomic": "ASE wrap",
}
_ORDER = ("metatomic", "rgpot_metatomic", "ase_metatomic")
# Muted teal / slate palette (readable on light Sphinx theme).
_COLORS = {
    "metatomic": "#0d7377",
    "rgpot_metatomic": "#14919b",
    "ase_metatomic": "#c44536",
}


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--json",
        type=Path,
        default=root / "docs/source/fig/data/metatomic_backend_bench.json",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=root / "docs/source/fig/metatomic_backend_bench.svg",
    )
    args = ap.parse_args(argv)

    payload = json.loads(args.json.read_text())
    by_key = {r["backend"]: r for r in payload["results"] if r.get("ok")}
    # normalize aliases
    if "metatomic_fat" in by_key and "metatomic" not in by_key:
        by_key["metatomic"] = by_key["metatomic_fat"]

    keys = [k for k in _ORDER if k in by_key]
    if len(keys) < 2:
        raise SystemExit(f"need ≥2 ok backends in {args.json}; got {list(by_key)}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    times_ms = [by_key[k]["seconds_per_call"] * 1e3 for k in keys]
    energies = [by_key[k]["energy"] for k in keys]
    labels = [_LABELS.get(k, k) for k in keys]
    colors = [_COLORS.get(k, "#555555") for k in keys]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.2, 3.8),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.28},
    )

    # --- timing ---
    ax = axes[0]
    x = range(len(keys))
    bars = ax.bar(x, times_ms, color=colors, width=0.62, edgecolor="white", linewidth=0.6)
    ax.set_xticks(list(x), labels, fontsize=9)
    ax.set_ylabel("wall time / single-point (ms)", fontsize=10)
    ax.set_title("Force-call cost (CPU, 14 atoms)", fontsize=11, pad=8)
    ax.set_ylim(0, max(times_ms) * 1.18)
    for b, t in zip(bars, times_ms):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height(),
            f"{t:.1f}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.set_axisbelow(True)

    # --- energy agreement ---
    ax = axes[1]
    e0 = energies[0]
    dE = [(e - e0) * 1e6 for e in energies]  # µeV vs fat
    bars = ax.bar(x, dE, color=colors, width=0.62, edgecolor="white", linewidth=0.6)
    ax.axhline(0.0, color="#333333", linewidth=0.8)
    ax.set_xticks(list(x), labels, fontsize=9)
    ax.set_ylabel(r"$\Delta E$ vs fat (µeV)", fontsize=10)
    ax.set_title("Energy agreement (PET-MAD s v1.5)", fontsize=11, pad=8)
    ymax = max(abs(v) for v in dE) if any(dE) else 1.0
    ax.set_ylim(-ymax * 1.4 - 1e-3, ymax * 1.4 + 1e-3)
    for b, v in zip(bars, dE):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height(),
            f"{v:+.1f}" if abs(v) >= 0.05 else "0",
            ha="center",
            va="bottom" if v >= 0 else "top",
            fontsize=9,
            fontweight="bold",
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.45)
    ax.set_axisbelow(True)

    n_atoms = by_key[keys[0]].get("n_atoms", "?")
    model = Path(payload.get("model", "model.pt")).name
    fig.suptitle(
        f"pyeonclient metatomic backends · {model} · n={n_atoms} · CPU",
        fontsize=11,
    )
    # short legend under figure
    fig.text(
        0.5,
        -0.02,
        "fat = PotType.METATOMIC  ·  RGPOT = dlopen libmetatomic_engine  ·  ASE = MetatomicCalculator wrap",
        ha="center",
        va="top",
        fontsize=8,
        color="#444444",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    # also PNG for browsers that prefer raster
    png = args.out.with_suffix(".png")
    fig.savefig(png, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", args.out)
    print("wrote", png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
