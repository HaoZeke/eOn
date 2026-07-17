"""NEB and minimization figures (library surface).

Cookbook / notebook code should call :func:`plot_neb` / :func:`plot_min` with
keyword arguments — not argv lists or ``python -m rgpycrumbs``.

Requires optional deps::

    pip install 'pyeonclient[plot]'

Implementation uses the chemparseplot / rgpycrumbs figure stack; the CLI is
not part of the public contract.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np

PlotTypeNeb = Literal["profile", "landscape"]
PlotTypeMin = Literal["profile", "landscape", "convergence"]
StructuresNeb = Literal["none", "all", "crit_points"]
StructuresMin = Literal["none", "endpoints", "all"]


def _require_plot_stack() -> None:
    try:
        import chemparseplot  # noqa: F401
        import rgpycrumbs  # noqa: F401
    except ImportError as exc:  # pragma: no cover - optional extra
        raise ImportError(
            "plot_neb / plot_min require the plot extra. "
            "Install with: pip install 'pyeonclient[plot]'"
        ) from exc


def _flag(name: str, value: bool) -> list[str]:
    # click is_flag pairs: --foo / --no-foo
    bare = name.lstrip("-").replace("_", "-")
    return [f"--{bare}" if value else f"--no-{bare}"]


def plot_neb(
    plot_type: PlotTypeNeb = "profile",
    *,
    con_file: str | Path = "neb.con",
    output: str | Path = "neb.png",
    title: str = "NEB Path",
    figsize: tuple[float, float] = (12.0, 8.0),
    zoom_ratio: float = 0.35,
    plot_structures: StructuresNeb = "all",
    strip_renderer: str = "xyzrender",
    xyzrender_config: str = "paton",
    rotation: str = "90x,0y,0z",
    facecolor: str = "white",
    fontsize_base: int = 14,
    show_pts: bool = True,
    highlight_last: bool = True,
    show_legend: bool = True,
    surface_type: str = "grad_imq",
    landscape_path: Literal["last", "all"] = "all",
    project_path: bool = True,
    rc_mode: Literal["path", "rmsd", "index"] = "path",
    landscape_mode: Literal["path", "surface"] = "surface",
    strip_spacing: float = 2.0,
    strip_dividers: bool = True,
    **extra: Any,
) -> Path:
    """Write a 1D energy profile or 2D RMSD landscape for an eOn NEB run.

    Reads ``neb_*.dat`` / path movies from the current directory (or patterns
    passed via *extra*). Returns the output path.
    """
    _require_plot_stack()
    from rgpycrumbs.eon.plt_neb import main as plt_neb_main

    out = Path(output)
    args: list[str] = [
        "--plot-type",
        plot_type,
        "--con-file",
        str(con_file),
        "--output-file",
        str(out),
        "--title",
        title,
        "--figsize",
        str(figsize[0]),
        str(figsize[1]),
        "--zoom-ratio",
        str(zoom_ratio),
        "--plot-structures",
        plot_structures,
        "--strip-renderer",
        strip_renderer,
        "--xyzrender-config",
        xyzrender_config,
        "--rotation",
        rotation,
        "--facecolor",
        facecolor,
        "--fontsize-base",
        str(fontsize_base),
        "--strip-spacing",
        str(strip_spacing),
        *_flag("show-pts", show_pts),
        *_flag("highlight-last", highlight_last),
        *_flag("show-legend", show_legend),
    ]
    if strip_dividers:
        args.append("--strip-dividers")
    if plot_type == "landscape":
        args.extend(
            [
                "--rc-mode",
                rc_mode,
                "--landscape-mode",
                landscape_mode,
                "--landscape-path",
                landscape_path,
                "--surface-type",
                surface_type,
                *_flag("project-path", project_path),
            ]
        )
    for key, val in extra.items():
        flag = f"--{key.replace('_', '-')}"
        if isinstance(val, bool):
            args.extend(_flag(key.replace("_", "-"), val))
        elif isinstance(val, (list, tuple)):
            args.append(flag)
            args.extend(str(v) for v in val)
        elif val is not None:
            args.extend([flag, str(val)])

    plt_neb_main.main(args, standalone_mode=False)
    return out


def plot_min(
    job_dirs: Sequence[str | Path],
    *,
    plot_type: PlotTypeMin = "profile",
    labels: Sequence[str] | None = None,
    output: str | Path = "min.png",
    surface_type: str = "grad_imq",
    project_path: bool = True,
    plot_structures: StructuresMin = "endpoints",
    strip_renderer: str = "xyzrender",
    xyzrender_config: str = "paton",
    rotation: str = "90x,0y,0z",
    strip_dividers: bool = True,
    **extra: Any,
) -> Path:
    """Write a minimization profile, landscape, or convergence figure."""
    _require_plot_stack()
    from rgpycrumbs.eon.plt_min import main as plt_min_main

    out = Path(output)
    dirs = [Path(d) for d in job_dirs]
    if labels is None:
        labels = [d.name for d in dirs]
    if len(labels) != len(dirs):
        raise ValueError("labels must match job_dirs length")

    args: list[str] = [
        "--plot-type",
        plot_type,
        "-o",
        str(out),
        "--surface-type",
        surface_type,
        "--plot-structures",
        plot_structures,
        "--strip-renderer",
        strip_renderer,
        "--xyzrender-config",
        xyzrender_config,
        "--rotation",
        rotation,
        *_flag("project-path", project_path),
    ]
    if strip_dividers:
        args.append("--strip-dividers")
    for d, lab in zip(dirs, labels, strict=True):
        args.extend(["--job-dir", str(d), "--label", lab])
    for key, val in extra.items():
        flag = f"--{key.replace('_', '-')}"
        if isinstance(val, bool):
            args.extend(_flag(key.replace("_", "-"), val))
        elif isinstance(val, (list, tuple)):
            args.append(flag)
            args.extend(str(v) for v in val)
        elif val is not None:
            args.extend([flag, str(val)])

    plt_min_main.main(args, standalone_mode=False)
    return out


def thin_min_movie(
    job_dir: str | Path,
    *,
    max_frames: int = 64,
    prefix: str = "minimization",
) -> int:
    """Thin a dense eOn minimization movie before landscape surface fits.

    ``write_movies`` records every force evaluation; long LBFGS paths can
    exceed ~150 frames and destabilize gradient surface fits. Keeps first and
    last frames plus evenly spaced intermediates.

    Returns the frame count after thinning (or the original count if no
    thinning was needed).
    """
    import readcon

    job_dir = Path(job_dir)
    movie = None
    for candidate in (job_dir / prefix, job_dir / f"{prefix}.con"):
        if candidate.exists():
            movie = candidate
            break
    if movie is None:
        return 0

    frames = list(readcon.read_con(str(movie)))
    n = len(frames)
    if n <= max_frames:
        return n

    idx = np.unique(np.linspace(0, n - 1, num=max_frames, dtype=int))
    if idx[-1] != n - 1:
        idx = np.unique(np.append(idx, n - 1))
    thinned = [frames[i] for i in idx]
    readcon.write_con(str(movie), thinned)

    dat_path = job_dir / f"{prefix}.dat"
    if dat_path.exists():
        lines = dat_path.read_text().splitlines()
        if lines:
            header, rows = lines[0], lines[1:]
            if len(rows) == n:
                kept = [rows[i] for i in idx]
                dat_path.write_text(header + "\n" + "\n".join(kept) + "\n")
    return len(thinned)


__all__ = [
    "plot_neb",
    "plot_min",
    "thin_min_movie",
]
