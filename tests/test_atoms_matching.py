"""Structure matching against the serial versions it replaced.

identical, rot_match and get_mappings decide whether two configurations are
the same state, so a faster version that changes a verdict is a defect. Each
case is run through both implementations and the verdicts compared.
"""

from __future__ import annotations

import numpy as np
import pytest

from eon import atoms as at
from tests.reference.atoms_pre_vectorization import (
    get_mappings_scalar,
    identical_scalar,
    rot_match_scalar,
)


def _rand(rng, n, scale=9.0, elements=("Cu",), skew=0.0):
    box = np.diag(np.full(3, float(scale)))
    if skew:
        box = np.triu(box + rng.normal(scale=skew, size=(3, 3)))
    p = at.Atoms(n)
    p.r = rng.random((n, 3)) @ box
    p.box = box
    p.free = np.ones(n)
    p.names = [elements[int(k)] for k in rng.integers(0, len(elements), n)]
    p.mass = np.full(n, 63.5)
    return p


def _fcc(reps=(2, 2, 2), a=3.615, rng=None, jitter=0.0):
    basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    pts = []
    for i in range(reps[0]):
        for j in range(reps[1]):
            for k in range(reps[2]):
                for b in basis:
                    pts.append((np.array([i, j, k], dtype=float) + b) * a)
    r = np.asarray(pts)
    if jitter:
        r = r + rng.normal(scale=jitter, size=r.shape)
    p = at.Atoms(len(r))
    p.r = r
    p.box = np.diag(np.asarray(reps, dtype=float) * a)
    p.free = np.ones(len(r))
    p.names = ["Cu"] * len(r)
    p.mass = np.full(len(r), 63.5)
    return p


def _pairs():
    """(label, a, b) covering same, moved, swapped, permuted and translated."""
    rng = np.random.default_rng(5)
    for n, elems, skew in ((2, ("Cu",), 0.0), (9, ("Cu",), 0.0),
                           (12, ("Cu", "Ag", "Au"), 0.0), (10, ("Cu", "Ag"), 1.0)):
        p = _rand(rng, n, elements=elems, skew=skew)
        yield "n%d same" % n, p, p.copy()

        q = p.copy()
        q.r = q.r + rng.normal(scale=0.01, size=q.r.shape)
        yield "n%d jitter" % n, p, q

        q = p.copy()
        q.r = q.r.copy()
        q.r[0] = q.r[0] + np.array([0.7, 0.3, 0.0])
        yield "n%d one moved" % n, p, q

        q = p.copy()
        q.r = q.r.copy()
        q.r[[0, 1]] = q.r[[1, 0]]
        q.names = list(q.names)
        q.names[0], q.names[1] = q.names[1], q.names[0]
        yield "n%d swapped" % n, p, q

        idx = rng.permutation(n)
        q = p.copy()
        q.r = p.r[idx].copy()
        q.names = [p.names[i] for i in idx]
        yield "n%d permuted" % n, p, q

        q = p.copy()
        q.r = (q.r + np.array([0.12, 0.05, 0.03])) % np.diag(p.box)
        yield "n%d translated" % n, p, q

        q = p.copy()
        q.box = q.box * 1.01
        yield "n%d box changed" % n, p, q

    a = _fcc(rng=rng, jitter=0.01)
    idx = rng.permutation(len(a))
    b = a.copy()
    b.r = a.r[idx].copy()
    yield "fcc permuted", a, b
    yield "fcc vs itself", a, a.copy()


PAIRS = list(_pairs())
IDS = [c[0] for c in PAIRS]
EPSILONS = (0.05, 0.2, 0.6)


@pytest.mark.parametrize("label,a,b", PAIRS, ids=IDS)
@pytest.mark.parametrize("eps", EPSILONS)
def test_identical_matches_scalar(label, a, b, eps):
    assert at.identical(a, b, eps) == identical_scalar(a, b, eps)


@pytest.mark.parametrize("label,a,b", PAIRS, ids=IDS)
@pytest.mark.parametrize("eps", EPSILONS)
def test_rot_match_matches_scalar(label, a, b, eps):
    assert at.rot_match(a, b, eps) == rot_match_scalar(a, b, eps)


@pytest.mark.parametrize("label,a,b", PAIRS, ids=IDS)
@pytest.mark.parametrize("eps", EPSILONS)
def test_get_mappings_matches_scalar(label, a, b, eps):
    for cutoff in (2.0, 4.0):
        assert at.get_mappings(a, b, eps, cutoff) == get_mappings_scalar(
            a, b, eps, cutoff
        )


def test_get_mappings_returns_plain_ints():
    rng = np.random.default_rng(2)
    a = _rand(rng, 6)
    got = at.get_mappings(a, a.copy(), 0.2, 3.0)
    assert got is not None
    assert all(type(k) is int and type(v) is int for k, v in got.items())


def test_get_mappings_grows_a_supplied_mapping():
    rng = np.random.default_rng(4)
    a = _rand(rng, 7)
    b = a.copy()
    for seed in ({0: 0}, {2: 2}):
        assert at.get_mappings(a, b, 0.2, 3.0, dict(seed)) == get_mappings_scalar(
            a, b, 0.2, 3.0, dict(seed)
        )


def test_get_mappings_one_atom_reports_no_mapping():
    """The serial version walked off the end of a.r here instead of answering."""
    rng = np.random.default_rng(6)
    a = _rand(rng, 1)
    assert at.get_mappings(a, a.copy(), 0.2, 3.0) is None
    with pytest.raises(IndexError):
        get_mappings_scalar(a, a.copy(), 0.2, 3.0)


def test_identical_needs_matching_elements_on_a_swap():
    """Two atoms swapped is identical only when they are the same element."""
    rng = np.random.default_rng(8)
    a = _rand(rng, 6, elements=("Cu",))
    b = a.copy()
    b.r = b.r.copy()
    b.r[[0, 1]] = b.r[[1, 0]]
    assert at.identical(a, b, 0.2)
    b.names = list(b.names)
    b.names[0] = "Ag"
    assert not at.identical(a, b, 0.2)
