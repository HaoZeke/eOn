"""Pre-#368 scalar atoms helpers, extracted from commit ee3012498 (parent of #368).

These are the production implementations that shipped before vectorization.
Golden tests compare live eon.atoms outputs to these functions on the same inputs.
Do not edit algorithmically without re-extracting from git history.
"""
from __future__ import annotations

import numpy

from eon.atoms import elements, pbc, per_atom_norm


def free_r_scalar(atoms) -> numpy.ndarray:
    """Historical Atoms.free_r (ee3012498).

    ``int(...)`` on the free count is required under current NumPy: ``sum`` on a
    float free array yields ``numpy.float64``, which is not a valid shape dim.
    Semantics match the pre-#368 loop (truthy free flags).
    """
    nfree = int(sum(atoms.free))
    temp = numpy.zeros((nfree, 3))
    index = 0
    for i in range(len(atoms.r)):
        if atoms.free[i]:
            temp[index] = atoms.r[i]
            index += 1
    return temp


def get_process_atoms_scalar(r, p, epsilon_r=0.2, nshells=1):
    """Historical get_process_atoms (ee3012498)."""
    mobileAtoms = []
    r2p = per_atom_norm(p.r - r.r, r.box)
    for i in range(len(r)):
        if r2p[i] > epsilon_r:
            mobileAtoms.append(i)
    if len(mobileAtoms) == 0:
        mobileAtoms.append(list(r2p).index(max(r2p)))
    neighborAtoms = []
    for atom in mobileAtoms:
        r1 = elements[r.names[atom]]["radius"]
        for i in range(len(r)):
            if i in mobileAtoms or i in neighborAtoms:
                continue
            r2 = elements[r.names[i]]["radius"]
            maxDist = (r1 + r2) * (1.0 + 0.2) * nshells
            if numpy.linalg.norm(pbc(p.r[atom] - p.r[i], p.box)) < maxDist:
                neighborAtoms.append(i)
    return mobileAtoms + neighborAtoms


def brute_neighbor_list_scalar(p, cutoff):
    """Historical brute_neighbor_list (ee3012498)."""
    nl = []
    ibox = numpy.linalg.inv(p.box)
    for a in range(len(p)):
        nl.append([])
        for b in range(len(p)):
            if b != a:
                dist = numpy.linalg.norm(pbc(p.r[a] - p.r[b], p.box, ibox))
                if dist < cutoff:
                    nl[a].append(b)
    return nl


def cna_scalar(n, nl):
    """Historical cna body (8de45a36b), taking a prebuilt neighbor list."""
    can_values = numpy.zeros(n)
    nr_FCC = numpy.zeros(n)
    nr_HCP = numpy.zeros(n)

    for a2 in range(n):
        nl_a2 = nl[a2]
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2]
            if a1 < a2:
                common = []
                nl_a1 = nl[a1]
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                if len(common) == 4:
                    bonds_nr = 0
                    bonds_sum = 0
                    for j2 in range(1, 4):
                        nl_j2 = nl[common[j2]]
                        for j1 in range(j2):
                            for k in range(len(nl_j2)):
                                if common[j1] == nl_j2[k]:
                                    bonds_nr += 1
                                    bonds_sum += j1 + j2

                    if bonds_nr == 2:
                        if bonds_sum == 6:
                            nr_FCC[a1] += 1
                            nr_FCC[a2] += 1
                        else:
                            nr_HCP[a1] += 1
                            nr_HCP[a2] += 1

    for i in range(n):
        if len(nl[i]) == 12:
            if nr_FCC[i] == 12:
                can_values[i] = 2
            elif (nr_FCC[i] == 6) and (nr_HCP[i] == 6):
                can_values[i] = 1
    return can_values


def cnat_scalar(n, nl):
    """Historical cnat body (8de45a36b), taking a prebuilt neighbor list."""
    can_values = numpy.zeros(n)
    nr_5 = numpy.zeros(n)
    nr_6 = numpy.zeros(n)

    for a2 in range(n):
        nl_a2 = nl[a2]
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2]
            if a1 < a2:
                common = []
                nl_a1 = nl[a1]
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                if len(common) in [5, 6]:
                    bonds_nr = 0
                    bonds_sum = 0
                    for j2 in range(1, len(common)):
                        nl_j2 = nl[common[j2]]
                        for j1 in range(j2):
                            for k in range(len(nl_j2)):
                                if common[j1] == nl_j2[k]:
                                    bonds_nr += 1
                                    bonds_sum += j1 + j2

                    if bonds_nr == 5 and len(common) == 5:
                        nr_5[a1] += 1
                        nr_5[a2] += 1
                    elif bonds_nr == 6 and len(common) == 6:
                        nr_6[a1] += 1
                        nr_6[a2] += 1

    for i in range(n):
        if nr_5[i] == 12:
            if len(nl[i]) == 12:
                can_values[i] = 1
            if (len(nl[i]) == 14) and (nr_6[i] == 2):
                can_values[i] = 2
            if (len(nl[i]) == 15) and (nr_6[i] == 3):
                can_values[i] = 3
            if (len(nl[i]) == 16) and (nr_6[i] == 4):
                can_values[i] = 4
        if nr_5[i] == 0:
            if nr_6[i] == 8 and len(nl[i]) == 14:
                can_values[i] = 5
    return can_values


def cnar_scalar(n, nl):
    """Historical cnar body (8de45a36b), taking a prebuilt neighbor list."""
    cna = {}
    for i in range(n):
        cna[i] = {}

    def codeString(j, k, l):
        return "%d,%d,%d" % (j, k, l)

    for a2 in range(n):
        nl_a2 = nl[a2]
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2]
            if a1 < a2:
                common = []
                nl_a1 = nl[a1]
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                bonds_nr = 0
                bonds_sum = 0
                for j2 in range(1, len(common)):
                    nl_j2 = nl[common[j2]]
                    for j1 in range(j2):
                        for k in range(len(nl_j2)):
                            if common[j1] == nl_j2[k]:
                                bonds_nr += 1
                                bonds_sum += j1 + j2
                code = codeString(len(common), bonds_nr, bonds_sum)
                if code not in cna[a1]:
                    cna[a1][code] = 0
                cna[a1][code] += 1
                if code not in cna[a2]:
                    cna[a2][code] = 0
                cna[a2][code] += 1

    return cna


def pair_lists_scalar(n, i, j):
    """Historical _pair_lists (8de45a36b): per-atom sets, then sorted."""
    nl_sets = [set() for _ in range(n)]
    for a, b in zip(i.tolist(), j.tolist()):
        ai, bi = int(a), int(b)
        if ai == bi:
            continue
        nl_sets[ai].add(bi)
    return [sorted(s) for s in nl_sets]
