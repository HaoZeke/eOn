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


def per_atom_norm_gen(v, box, ibox: Optional[np.ndarray] = None):
    """Yield per-row norms after PBC (legacy generator API)."""
    diff = pbc(v, box, ibox)
    for d in diff:
        yield np.linalg.norm(d)
