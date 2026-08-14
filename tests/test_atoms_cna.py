"""CNA kernels against the serial loops they replaced.

The CNA verdict decides which atoms a displacement may pick, so the vectorized
sweep has to return the same array as the pre-vectorization loops for the same
structure, down to the last atom. Cases cover perfect fcc/hcp/bcc, vacancies
and thermal noise, cutoffs past half the cell (vesin reports the same neighbor
under several periodic shifts), empty and one-atom cells, and a zero cutoff.
"""

from __future__ import annotations

import numpy as np
import pytest

from eon import atoms as at
from tests.reference.atoms_pre_vectorization import (
    cna_scalar,
    cnar_scalar,
    cnat_scalar,
    pair_lists_scalar,
)

FCC_BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
BCC_BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]


def _lattice(basis, cell, reps, rng=None, jitter=0.0, vacancies=0):
    cell = np.asarray(cell, dtype=float)
    pts = []
    for i in range(reps[0]):
        for j in range(reps[1]):
            for k in range(reps[2]):
                for b in basis:
                    pts.append((np.array([i, j, k], dtype=float) + b) * cell)
    r = np.asarray(pts)
    if vacancies:
        keep = np.sort(rng.permutation(len(r))[: len(r) - vacancies])
        r = r[keep]
    if jitter:
        r = r + rng.normal(scale=jitter, size=r.shape)
    p = at.Atoms(len(r))
    p.r = r
    p.box = np.diag(np.asarray(reps, dtype=float) * cell)
    p.free = np.ones(len(r))
    p.names = ["Cu"] * len(r)
    p.mass = np.full(len(r), 63.5)
    return p


def _fcc(reps=(3, 3, 3), a=3.615, **kw):
    return _lattice(FCC_BASIS, [a, a, a], reps, **kw)


def _bcc(reps=(3, 3, 3), a=2.87, **kw):
    return _lattice(BCC_BASIS, [a, a, a], reps, **kw)


def _hcp(reps=(3, 3, 2), a=2.55):
    cell = [a, a * np.sqrt(3.0), a * np.sqrt(8.0 / 3.0)]
    basis = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 1.0 / 6.0, 0.5),
        (0.0, 2.0 / 3.0, 0.5),
    ]
    return _lattice(basis, cell, reps)


def _random_cell(rng, n, scale, skew=0.0):
    box = np.diag(np.full(3, float(scale)))
    if skew:
        box = np.triu(box + rng.normal(scale=skew, size=(3, 3)))
    p = at.Atoms(n)
    p.r = rng.random((n, 3)) @ box
    p.box = box
    p.free = np.ones(n)
    p.names = ["Cu"] * n
    p.mass = np.full(n, 63.5)
    return p


def _cases():
    rng = np.random.default_rng(11)
    yield "fcc perfect", _fcc(), 3.0
    yield "fcc wide cutoff", _fcc(), 4.2
    yield "fcc vacancies", _fcc(rng=rng, vacancies=4), 3.0
    yield "fcc thermal", _fcc(rng=rng, jitter=0.08), 3.0
    yield "fcc melted", _fcc(rng=rng, jitter=0.4), 3.2
    yield "bcc", _bcc(), 3.0
    yield "bcc thermal", _bcc(rng=rng, jitter=0.05), 3.0
    yield "hcp", _hcp(), 3.0
    # Cutoff past half the cell: the same neighbor arrives under several shifts
    yield "fcc multi-image", _fcc((2, 2, 2)), 7.0
    yield "gas", _random_cell(rng, 40, 20.0), 2.0
    yield "triclinic", _random_cell(rng, 60, 12.0, skew=1.5), 4.0
    yield "one atom", _random_cell(rng, 1, 8.0), 2.5
    yield "empty", _random_cell(rng, 0, 8.0), 2.5
    yield "zero cutoff", _fcc((2, 2, 2)), 0.0


CASES = list(_cases())
IDS = [c[0] for c in CASES]


@pytest.mark.parametrize("label,p,cutoff", CASES, ids=IDS)
def test_cna_matches_scalar(label, p, cutoff):
    nl = at.neighbor_list(p, cutoff)
    assert np.array_equal(at.cna(p, cutoff), cna_scalar(len(p), nl))


@pytest.mark.parametrize("label,p,cutoff", CASES, ids=IDS)
def test_cnat_matches_scalar(label, p, cutoff):
    nl = at.neighbor_list(p, cutoff)
    assert np.array_equal(at.cnat(p, cutoff), cnat_scalar(len(p), nl))


@pytest.mark.parametrize("label,p,cutoff", CASES, ids=IDS)
def test_cnar_matches_scalar(label, p, cutoff):
    nl = at.neighbor_list(p, cutoff)
    assert at.cnar(p, cutoff) == cnar_scalar(len(p), nl)


@pytest.mark.parametrize("label,p,cutoff", CASES, ids=IDS)
def test_neighbor_list_matches_scalar_build(label, p, cutoff):
    """The CSR pack must hold what the per-atom set build held."""
    n = len(p)
    if n == 0 or cutoff <= 0:
        assert at.neighbor_list(p, cutoff) == [[] for _ in range(n)]
        return
    from vesin import NeighborList

    i, j = NeighborList(cutoff=float(cutoff), full_list=True).compute(
        np.asarray(p.r, dtype=float),
        np.asarray(p.box, dtype=float),
        periodic=True,
        quantities="ij",
    )
    want = pair_lists_scalar(n, np.asarray(i), np.asarray(j))
    assert at.neighbor_list(p, cutoff) == want
    assert at.coordination_numbers(p, cutoff) == [len(x) for x in want]


def test_cna_classifies_the_three_lattices():
    """2 is fcc (421), 1 is hcp (422), cnat 5 is bcc."""
    assert np.all(at.cna(_fcc((4, 4, 4)), 3.0) == 2)
    assert np.all(at.cna(_hcp((4, 4, 3)), 3.0) == 1)
    assert np.all(at.cnat(_bcc((4, 4, 4)), 3.0) == 5)


def test_not_hcp_or_fcc_returns_plain_ints():
    p = _fcc(rng=np.random.default_rng(3), vacancies=4)
    idx = at.not_HCP_or_FCC(p, 3.0)
    assert idx == [i for i, v in enumerate(at.cna(p, 3.0)) if v == 0]
    assert all(type(i) is int for i in idx)


def test_not_tcp_variants_return_plain_ints():
    p = _bcc(rng=np.random.default_rng(5), jitter=0.15)
    values = at.cnat(p, 3.0)
    assert at.not_TCP(p, 3.0) == [i for i, v in enumerate(values) if v in (0, 5)]
    assert at.not_TCP_or_BCC(p, 3.0) == [i for i, v in enumerate(values) if v == 0]
    assert all(type(i) is int for i in at.not_TCP(p, 3.0))


def test_neighbor_csr_and_list_agree():
    p = _fcc(rng=np.random.default_rng(9), jitter=0.05)
    indptr, indices = at.neighbor_csr(p, 3.0)
    nl = at.neighbor_list(p, 3.0)
    assert indptr[-1] == sum(len(x) for x in nl)
    for a, neighbors in enumerate(nl):
        assert indices[indptr[a]:indptr[a + 1]].tolist() == neighbors
