"""Cell matrix helpers (row-vector lattice, ASE/vesin convention)."""

from __future__ import annotations

import numpy as np


def length_angle_to_box(boxlengths, angles) -> np.ndarray:
    """Convert cell lengths + angles (degrees) to a 3×3 row-vector box matrix."""
    box = np.zeros((3, 3), dtype=float)
    ang = np.asarray(angles, dtype=float) * (np.pi / 180.0)
    lengths = np.asarray(boxlengths, dtype=float)
    box[0][0] = 1.0
    box[1][0] = np.cos(ang[0])
    box[1][1] = np.sin(ang[0])
    box[2][0] = np.cos(ang[1])
    box[2][1] = (np.cos(ang[2]) - box[1][0] * box[2][0]) / box[1][1]
    box[2][2] = np.sqrt(1.0 - box[2][0] ** 2 - box[2][1] ** 2)
    box[0, :] *= lengths[0]
    box[1, :] *= lengths[1]
    box[2, :] *= lengths[2]
    return box


def box_to_length_angle(box) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`length_angle_to_box`."""
    box = np.asarray(box, dtype=float)
    lengths = np.zeros(3, dtype=float)
    lengths[0] = np.linalg.norm(box[0, :])
    lengths[1] = np.linalg.norm(box[1, :])
    lengths[2] = np.linalg.norm(box[2, :])
    angles = np.zeros(3, dtype=float)
    angles[0] = np.arccos(
        np.dot(box[0, :] / lengths[0], box[1, :] / lengths[1])
    )
    angles[1] = np.arccos(
        np.dot(box[0, :] / lengths[0], box[2, :] / lengths[2])
    )
    angles[2] = np.arccos(
        np.dot(box[1, :] / lengths[1], box[2, :] / lengths[2])
    )
    angles *= 180.0 / np.pi
    return lengths, angles
