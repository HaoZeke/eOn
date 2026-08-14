"""Periodic boundary conditions and distance helpers."""

from __future__ import annotations

from typing import Optional

import numpy as np


def pbc(r, box, ibox: Optional[np.ndarray] = None) -> np.ndarray:
    """Minimum-image convention for displacement(s).

    Parameters
    ----------
    r : (3,) or (N, 3)
        Displacement vector(s).
    box : (3, 3)
        Cell matrix (rows = lattice vectors).
    ibox : (3, 3), optional
        Inverse of *box*; computed if omitted.
    """
    if ibox is None:
        ibox = np.linalg.inv(box)
    vdir = np.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return np.dot(vdir, box)


def per_atom_norm(v, box, ibox: Optional[np.ndarray] = None) -> np.ndarray:
    """Per-row Euclidean norm after PBC (shape ``(N,)``)."""
    diff = pbc(v, box, ibox)
    return np.sqrt(np.sum(diff**2.0, axis=1))


def per_atom_norm_rows(a, b, box, ibox: Optional[np.ndarray] = None) -> np.ndarray:
    """MIC distance from every row of *b* to every row of *a*.

    Shape ``(len(b), len(a))``; entry ``[i, j]`` is the norm of
    ``pbc(a[j] - b[i])``, so one call replaces a loop of
    :func:`per_atom_norm` on ``a - b[i]``.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    # Flattened to (rows, 3): numpy.dot takes a much slower path on a stacked
    # 3-D operand than on a plain matrix.
    diff = pbc((a[None, :, :] - b[:, None, :]).reshape(-1, 3), box, ibox)
    return np.sqrt(np.sum(diff**2.0, axis=1)).reshape(b.shape[0], a.shape[0])


def per_atom_norm_gen(v, box, ibox: Optional[np.ndarray] = None):
    """Yield per-row norms after PBC (legacy generator API)."""
    diff = pbc(v, box, ibox)
    for d in diff:
        yield np.linalg.norm(d)
