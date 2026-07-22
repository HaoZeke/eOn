"""Process-atom selection for recycling / active regions."""

from __future__ import annotations

from typing import List

import numpy as np

from eon.geometry.neighbors import neighbor_list
from eon.geometry.pbc import per_atom_norm


def get_process_atoms(r, p, epsilon_r: float = 0.2, nshells: int = 1) -> List[int]:
    """Indices of atoms that move significantly plus neighbor shells.

    Mobile atoms: displacement > *epsilon_r* (PBC), else the single largest
    mover. Neighbor shells use covalent radii from :data:`eon.atoms.elements`
    with cutoff ``(r1+r2)*(1.2)*nshells``, evaluated via a single vesin
    neighbor list at the max such cutoff, then filtered pair-wise.

    Returns plain Python ``int`` indices (safe for recycling metadata
    ``repr``/``eval``).
    """
    # Lazy import: elements table lives with the periodic table data
    from eon.atoms import elements

    r2p = per_atom_norm(p.r - r.r, r.box)
    mobile = [int(i) for i in np.flatnonzero(r2p > epsilon_r)]
    if not mobile:
        mobile = [int(np.argmax(r2p))]

    n = len(r)
    if nshells <= 0 or n == 0:
        return mobile

    # Max pairwise shell cutoff for a single vesin query
    radii = []
    for name in r.names:
        radii.append(float(elements[name]["radius"]))
    radii_a = np.asarray(radii, dtype=float)
    rmax = float(radii_a.max()) if len(radii_a) else 1.0
    shell_cutoff = (2.0 * rmax) * (1.0 + 0.2) * float(nshells)

    # Neighbors of all atoms at shell_cutoff on product geometry (historical: p.r)
    full_nl = neighbor_list(p, shell_cutoff)

    mobile_set = set(mobile)
    neighbor_set = set()
    neighbors: List[int] = []
    for atom in mobile:
        r1 = radii_a[atom]
        for j in full_nl[atom]:
            if j in mobile_set or j in neighbor_set:
                continue
            r2 = radii_a[j]
            max_dist = (r1 + r2) * (1.0 + 0.2) * float(nshells)
            # Distance already constrained by shell_cutoff >= max_dist for all pairs
            # with r1,r2 <= rmax; still enforce pair threshold for smaller atoms.
            # Use stored NL: re-check distance with pbc for correctness when
            # shell_cutoff is loose.
            # Cheap: if max_dist >= shell_cutoff - eps, accept; else measure.
            if max_dist >= shell_cutoff - 1e-12:
                neighbor_set.add(j)
                neighbors.append(int(j))
            else:
                # precise check
                from eon.geometry.pbc import pbc

                d = np.linalg.norm(pbc(p.r[j] - p.r[atom], p.box))
                if d < max_dist:
                    neighbor_set.add(j)
                    neighbors.append(int(j))
    return mobile + neighbors
