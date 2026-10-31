"""Run examples/akmc-pt to its first KMC step and draw the tutorial figures.

From the repository root: ``python docs/source/fig/akmc/akmc_pt.py run DIR``
(needs eonclient and eon), then ``... plot DIR`` (needs matplotlib and ASE).
"""

import configparser
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parents[3] / "examples" / "akmc-pt"


def kmc_steps(run):
    dyn = run / "dynamics.txt"
    return dyn.read_text().splitlines()[2:] if dyn.is_file() else []


def run(rundir, max_calls=400):
    rundir.mkdir(parents=True)
    for name in ("config.ini", "pos.con"):
        shutil.copy(EXAMPLE / name, rundir / name)
    with open(rundir / "server.log", "w") as log:
        for call in range(1, max_calls + 1):
            log.write(f"### server call {call}\n")
            log.flush()
            subprocess.run([sys.executable, "-m", "eon.server"], cwd=rundir,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
            if kmc_steps(rundir):
                return
    raise SystemExit(f"no KMC step after {max_calls} server calls")


def read_calls(run):
    """(results registered so far, state, confidence) after each call."""
    calls, results = [], 0
    for block in (run / "server.log").read_text().split("### server call ")[1:]:
        results += sum(map(int, re.findall(r"Processed (\d+) results", block)))
        state, conf = re.findall(r"in state (\d+) with confidence ([0-9.]+)", block)[-1]
        calls.append((results, int(state), float(conf)))
    return calls


def plot(rundir):
    import ase.io
    import matplotlib
    import numpy as np
    from ase.visualize.plot import plot_atoms

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg = configparser.ConfigParser()
    cfg.read([rundir / "config.ini", rundir / "states" / "0" / "info"])
    target = cfg.getfloat("AKMC", "confidence")
    calls = read_calls(rundir)
    in0 = [c for c in calls if c[1] == 0]
    step_at = next(c[0] for c in calls if c[1] != 0)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.step([c[0] for c in in0], [c[2] for c in in0], where="post")
    ax.axhline(target, ls="--", color="grey")
    ax.text(0, target + 0.02, f"target {target}", color="grey")
    ax.axvline(step_at, ls=":", color="C3")
    ax.text(step_at, 0.3, "KMC step ", color="C3", ha="right")
    ax.set(ylim=(0, 1.05), xlabel="search results registered", ylabel="confidence in state 0")
    fig.tight_layout()
    fig.savefig(HERE / "akmc_pt_confidence.png", dpi=150)

    e0 = cfg.getfloat("MetaData", "reactant energy")
    table = (rundir / "states" / "0" / "processtable").read_text().splitlines()[1:]
    rows = [line.split() for line in table]
    chosen = kmc_steps(rundir)[0].split()[2]
    fig, ax = plt.subplots(figsize=(6, 4))
    for r in sorted(rows, key=lambda r: r[0] == chosen):
        y = [0.0, float(r[1]) - e0, float(r[4]) - e0]
        if r[0] == chosen:
            ax.plot(y, "o-", color="C3", lw=2.5, label=f"process taken, barrier {r[6]} eV")
        else:
            ax.plot(y, "o-", color="grey", lw=1, alpha=0.6)
    ax.set_xticks([0, 1, 2], ["reactant", "saddle", "product"])
    ax.set(ylabel="energy relative to state 0 (eV)",
           title=f"{len(rows)} processes in state 0, equal energies overlap")
    ax.legend(loc="lower center")
    fig.tight_layout()
    fig.savefig(HERE / "akmc_pt_processes.png", dpi=150)

    procdata = rundir / "states" / "0" / "procdata"
    kinds = ("reactant", "saddle", "product")
    frames = [ase.io.read(procdata / f"{k}_{chosen}.con", format="eon") for k in kinds]
    z = frames[0].positions[:, 2]
    shown = np.flatnonzero(z > z.max() - 3.5)
    colors = ["C1" if z[i] > z.max() - 1.0 else "lightgrey" for i in shown]
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.4))
    for ax, kind, atoms in zip(axes, kinds, frames):
        plot_atoms(atoms[shown], ax, radii=1.2, colors=colors)
        ax.set_title(kind)
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(HERE / "akmc_pt_process_taken.png", dpi=150)


if __name__ == "__main__":
    {"run": run, "plot": plot}[sys.argv[1]](Path(sys.argv[2]))
