"""Neighbor lists via vesin (PBC-aware).

Historically eOn used a Python sweep-and-prune or an O(N²) brute loop.
Both paths now go through :class:`vesin.NeighborList`, which is the supported
geometry kernel for pair finding. The *brute* flag is retained for API
compatibility but no longer selects a different algorithm.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple, Union

import numpy as np
from vesin import NeighborList as VesinNeighborList

from eon.geometry.pbc import pbc

# Structure-like: needs .r, .box, __len__, and optionally .free
StructureLike = object


def _positions_box(p: StructureLike):
    r = np.asarray(p.r, dtype=float)
    box = np.asarray(p.box, dtype=float)
    if r.ndim != 2 or r.shape[1] != 3:
        raise ValueError(f"positions must be (N, 3), got {r.shape}")
    if box.shape != (3, 3):
        raise ValueError(f"box must be (3, 3), got {box.shape}")
    return r, box


def unique_sorted(keys: np.ndarray) -> np.ndarray:
    """Ascending distinct values of an integer array.

    Sort plus an adjacent-difference mask. :func:`numpy.unique` computes the
    same thing an order of magnitude slower on integer keys of this size
    (0.22 s vs 0.018 s for 162k int64 under numpy 2.5).
    """
    keys = np.sort(keys)
    if keys.size == 0:
        return keys
    first = np.empty(keys.shape, dtype=bool)
    first[0] = True
    np.not_equal(keys[1:], keys[:-1], out=first[1:])
    return keys[first]


def _csr_from_pairs(
    n: int, i: np.ndarray, j: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Pack vesin pair indices into CSR adjacency ``(indptr, indices)``.

    Vesin may return the same neighbor index multiple times under different
    periodic shifts when the cutoff reaches multiple images. eOn stores unique
    atom indices only (no multi-image multiplicity), and self pairs are
    dropped.

    ``indices[indptr[a]:indptr[a + 1]]`` holds the neighbors of atom *a* in
    ascending order, so the flat key ``a * n + neighbor`` is ascending over the
    whole array and a membership test is one :func:`numpy.searchsorted`.
    """
    i = np.asarray(i, dtype=np.int64).reshape(-1)
    j = np.asarray(j, dtype=np.int64).reshape(-1)
    keep = i != j
    i = i[keep]
    j = j[keep]
    # One key per directed pair, sorted, which orders rows then columns.
    keys = unique_sorted(i * np.int64(n) + j)
    rows = keys // np.int64(n)
    cols = keys - rows * np.int64(n)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(rows, minlength=n), out=indptr[1:])
    return indptr, cols


def _pair_lists(n: int, i: np.ndarray, j: np.ndarray) -> List[List[int]]:
    """Convert vesin pair indices to eOn adjacency lists (plain Python ints)."""
    indptr, indices = _csr_from_pairs(n, i, j)
    flat = indices.tolist()
    ptr = indptr.tolist()
    return [flat[ptr[a]:ptr[a + 1]] for a in range(n)]


def neighbor_list(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,  # noqa: ARG001 — API compat; vesin always used
) -> List[List[int]]:
    """Return neighbors within *cutoff* for each atom (PBC, full undirected list).

    Parameters
    ----------
    p
        Structure with ``.r`` (N,3) and ``.box`` (3,3).
    cutoff
        Pair cutoff distance.
    brute
        Ignored; kept so callers using ``config.comp_brute_neighbors`` need no change.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0:
        return [[] for _ in range(n)]
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j = calc.compute(r, box, periodic=True, quantities="ij")
    return _pair_lists(n, np.asarray(i), np.asarray(j))


def neighbor_csr(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,  # noqa: ARG001 — API compat; vesin always used
) -> Tuple[np.ndarray, np.ndarray]:
    """Neighbors within *cutoff* as CSR arrays, without boxing into lists.

    Same content as :func:`neighbor_list` (PBC, full undirected, unique atom
    indices, ascending within a row) held as ``(indptr, indices)`` int64
    arrays. Kernels that walk the adjacency with numpy take this; callers that
    index per atom in Python take :func:`neighbor_list`.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    empty = np.zeros(0, dtype=np.int64)
    if n == 0:
        return np.zeros(1, dtype=np.int64), empty
    if cutoff <= 0:
        return np.zeros(n + 1, dtype=np.int64), empty
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j = calc.compute(r, box, periodic=True, quantities="ij")
    return _csr_from_pairs(n, np.asarray(i), np.asarray(j))


def brute_neighbor_list(p: StructureLike, cutoff: float) -> List[List[int]]:
    """Alias of :func:`neighbor_list` (historical name)."""
    return neighbor_list(p, cutoff, brute=True)


def neighbor_list_vectors(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,
) -> List[List[np.ndarray]]:
    """Neighbor list with minimum-image vectors from center → neighbor."""
    nl = neighbor_list(p, cutoff, brute=brute)
    r, box = _positions_box(p)
    ibox = np.linalg.inv(box)
    out: List[List[np.ndarray]] = []
    for center, neighs in enumerate(nl):
        vecs = []
        for j in neighs:
            vecs.append(pbc(r[j] - r[center], box, ibox))
        out.append(vecs)
    return out


def coordination_numbers(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    indptr, _ = neighbor_csr(p, cutoff, brute)
    return np.diff(indptr).tolist()


def least_coordinated(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    """Indices of free atoms with the lowest coordination number."""
    cn = coordination_numbers(p, cutoff, brute)
    if not cn:
        return []
    maxcoord = max(cn)
    mincoord = min(cn)
    free = np.asarray(getattr(p, "free", np.ones(len(cn))), dtype=bool)
    while mincoord <= maxcoord:
        least = [i for i in range(len(cn)) if cn[i] <= mincoord and free[i]]
        if least:
            return least
        mincoord += 1
    return []
