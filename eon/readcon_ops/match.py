"""Bijective structure match and rigid alignment.

These functions take any object with ``r``, ``box``, ``names``, ``copy``,
and ``__len__``. A ``readcon.ConFrame`` adapter can satisfy that later.
The minimum-image norm is the eOn kernel so a match here agrees with
:func:`eon.geometry.per_atom_norm`.
"""

from __future__ import annotations

import logging

import numpy

from eon.geometry.pbc import per_atom_norm

logger = logging.getLogger("readcon_ops")


def identical(atoms1, atoms2, epsilon_r: float) -> bool:
    """True when same-element atoms match within ``epsilon_r`` angstroms.

    An atom already matched by index stays taken. A second atom cannot
    claim that site.
    """
    if len(atoms1) != len(atoms2):
        return False

    for i in range(3):
        for j in range(3):
            if abs(atoms1.box[i][j] - atoms2.box[i][j]) > 0.0001:
                logger.warning(
                    "Identical returned false because boxes were not the same"
                )
                return False
    box = atoms1.box
    ibox = numpy.linalg.inv(box)

    mismatch = []
    pan = per_atom_norm(atoms1.r - atoms2.r, box, ibox)
    for i in range(len(pan)):
        if pan[i] > epsilon_r:
            mismatch.append(i)
        elif atoms1.names[i] != atoms2.names[i]:
            return False

    used = {i for i in range(len(atoms1)) if i not in mismatch}
    for i in mismatch:
        pan = per_atom_norm(atoms1.r - atoms2.r[i], box, ibox)
        best = None
        best_d = 1e300
        for j in range(len(pan)):
            if j in used:
                continue
            if (
                pan[j] < epsilon_r
                and pan[j] < best_d
                and atoms1.names[j] == atoms2.names[i]
            ):
                best = j
                best_d = pan[j]
        if best is None:
            return False
        used.add(best)
    return True


def get_rotation_matrix(axis, theta):
    axis = axis / numpy.linalg.norm(axis)
    t = theta
    ct = numpy.cos(t)
    st = numpy.sin(t)
    one_minus = 1.0 - ct
    rx, ry, rz = axis
    rotmat = numpy.zeros((3, 3))
    rotmat[0][0] = one_minus * rx * rx + ct
    rotmat[0][1] = one_minus * ry * rx + rz * st
    rotmat[0][2] = one_minus * rz * rx - ry * st
    rotmat[1][0] = one_minus * rx * ry - rz * st
    rotmat[1][1] = one_minus * ry * ry + ct
    rotmat[1][2] = one_minus * rz * ry + rx * st
    rotmat[2][0] = one_minus * rx * rz + ry * st
    rotmat[2][1] = one_minus * ry * rz - rx * st
    rotmat[2][2] = one_minus * rz * rz + ct
    return rotmat


def rotate(r, axis, center, angle):
    new_r = r.copy()
    if abs(angle) == 0.0:
        return new_r
    rotmat = get_rotation_matrix(axis, angle)
    center = center.copy()
    new_r -= center
    new_r = numpy.dot(new_r, rotmat)
    new_r += center
    return new_r


def internal_motion(a, b):
    """Return ``b`` with translation and rotation removed relative to ``a``."""
    b = b.copy()
    b.r += a.r[0] - b.r[0]
    a0a1 = (a.r[1] - a.r[0]) / numpy.linalg.norm(a.r[1] - a.r[0])
    b0b1 = (b.r[1] - b.r[0]) / numpy.linalg.norm(b.r[1] - b.r[0])
    cross1 = numpy.cross(b0b1, a0a1)
    norm1 = numpy.linalg.norm(cross1)
    if norm1 > 1e-12:
        axis1 = cross1 / norm1
        theta1 = numpy.arccos(numpy.clip((a0a1 * b0b1).sum(), -1.0, 1.0))
        b.r = rotate(b.r, axis1, a.r[0], theta1)
    axis2 = (a.r[2] - a.r[0]) / numpy.linalg.norm(a.r[2] - a.r[0])
    va = a.r[2] - ((a.r[2] - a.r[0]) * axis2).sum() * axis2
    vb = b.r[2] - ((b.r[2] - a.r[0]) * axis2).sum() * axis2
    nva = numpy.linalg.norm(va)
    nvb = numpy.linalg.norm(vb)
    if nva > 1e-12 and nvb > 1e-12:
        va = va / nva
        vb = vb / nvb
        cross2 = numpy.cross(vb, va)
        if numpy.linalg.norm(cross2) > 1e-12:
            theta2 = numpy.arccos(numpy.clip((va * vb).sum(), -1.0, 1.0))
            b.r = rotate(b.r, axis2, a.r[0], theta2)
    return b
