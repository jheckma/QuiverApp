"""Conjugacy classes and the (numerical) character table of a MatrixGroup.

Irreducible characters are obtained from the class algebra: the class-sum
matrices N_k (with structure constants c_{ijk}) commute and are simultaneously
diagonalisable, and their common eigenvectors are the irreducible characters
(Burnside / Dixon-Schneider, done here numerically over C).  This needs only
the multiplication table, so it works for any finite matrix group.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .groups import MatrixGroup, _key


@dataclass
class CharacterTable:
    group: MatrixGroup
    classes: list[list[int]]          # indices of elements in each class
    class_reps: list[int]             # one representative index per class
    class_sizes: list[int]
    class_orders: list[int]           # order of a representative element
    inverse_class: list[int]          # index of the class of g^{-1}
    chars: np.ndarray                 # chars[i, k] = chi_i(class k); row 0 = trivial
    chi_Q: np.ndarray                 # defining-rep character on each class

    @property
    def num_irreps(self) -> int:
        return self.chars.shape[0]

    @property
    def dims(self) -> list[int]:
        """Dimensions of the irreps = chi_i(identity class)."""
        e = self._identity_class()
        return [int(round(self.chars[i, e].real)) for i in range(self.num_irreps)]

    def _identity_class(self) -> int:
        for k, o in enumerate(self.class_orders):
            if o == 1:
                return k
        return 0


def _multiplication_table(g: MatrixGroup) -> tuple[dict, np.ndarray, list[int]]:
    """Return (key->index map, mult table, inverse list)."""
    idx = {_key(m): i for i, m in enumerate(g.elements)}
    n = len(g.elements)
    table = np.empty((n, n), dtype=np.int64)
    for i, a in enumerate(g.elements):
        for j, b in enumerate(g.elements):
            table[i, j] = idx[_key(a @ b)]
    inverse = [0] * n
    e_idx = idx[_key(np.eye(3, dtype=complex))]
    for i in range(n):
        for j in range(n):
            if table[i, j] == e_idx:
                inverse[i] = j
                break
    return idx, table, inverse


def _element_order(g: MatrixGroup, table: np.ndarray, e_idx: int, i: int) -> int:
    o, cur = 1, i
    while cur != e_idx:
        cur = table[cur, i]
        o += 1
        if o > len(g.elements):
            raise RuntimeError("element order exceeded group order")
    return o


def _conjugacy_classes(table: np.ndarray, inverse: list[int]) -> list[list[int]]:
    n = table.shape[0]
    seen = [False] * n
    classes = []
    for i in range(n):
        if seen[i]:
            continue
        orbit = set()
        for x in range(n):
            # x i x^{-1}
            cj = table[table[x, i], inverse[x]]
            orbit.add(cj)
        for c in orbit:
            seen[c] = True
        classes.append(sorted(orbit))
    return classes


def build_character_table(g: MatrixGroup) -> CharacterTable:
    idx, table, inverse = _multiplication_table(g)
    n = len(g.elements)
    e_idx = idx[_key(np.eye(3, dtype=complex))]

    classes = _conjugacy_classes(table, inverse)
    r = len(classes)
    class_of = [0] * n
    for k, cls in enumerate(classes):
        for x in cls:
            class_of[x] = k
    reps = [cls[0] for cls in classes]
    sizes = [len(cls) for cls in classes]
    orders = [_element_order(g, table, e_idx, rep) for rep in reps]
    inv_class = [class_of[inverse[rep]] for rep in reps]

    # class-algebra structure constants:
    #   N_k[i, j] = #{ (x in C_i) : x * y = z_j ,  for fixed y in C_k }
    # build via:  c_{k i j} = #{(x,y): x in C_k, y in C_i, xy in C_j} / |C_j|
    # We use the standard class matrices M_k with (M_k)_{ij} = number of
    # (a in C_i, b in C_k) with a b = (rep of C_j).
    Ms = []
    for k in range(r):
        M = np.zeros((r, r))
        Ck = classes[k]
        for j in range(r):
            zj = reps[j]
            # count a in C_i, b in C_k with a*b = zj  ->  a = zj * b^{-1}
            for b in Ck:
                a = table[zj, inverse[b]]
                M[class_of[a], j] += 1
        Ms.append(M)

    # simultaneous diagonalisation via a generic real combination
    coeff = np.array([np.cos(1.0 + 0.731 * k) for k in range(r)])
    Mcomb = sum(c * M for c, M in zip(coeff, Ms))
    _, vecs = np.linalg.eig(Mcomb)

    # For each eigenvector v (a common eigenvector of all M_k), the eigenvalue
    # of M_k is  omega_k = |C_k| chi(C_k) / chi(1).  Recover chi up to scale,
    # then normalize with  sum_k |C_k| |chi(C_k)|^2 = |G|.
    chars = np.zeros((r, r), dtype=complex)
    found = []
    for col in range(vecs.shape[1]):
        v = vecs[:, col]
        omegas = np.array([
            (Ms[k] @ v)[np.argmax(np.abs(v))] / v[np.argmax(np.abs(v))]
            for k in range(r)
        ])
        # chi(C_k) proportional to omega_k / |C_k|; fix chi(1) > 0 real
        ratio = omegas / np.array(sizes)
        # chi(1)^2 = |G| / sum_k |C_k| |ratio_k|^2  (ratio already chi/chi(1))
        denom = sum(sizes[k] * abs(ratio[k]) ** 2 for k in range(r))
        if denom < 1e-9:
            continue
        d = np.sqrt(n / denom)
        chi = ratio * d
        # canonicalise phase so chi(identity) is positive real
        e_class = class_of[e_idx]
        if abs(chi[e_class]) > 1e-9:
            chi = chi * (abs(chi[e_class]) / chi[e_class])
        found.append(chi)

    # deduplicate (numerical) and keep exactly r distinct irreducible characters
    uniq = []
    for chi in found:
        if not any(np.allclose(chi, u, atol=1e-4) for u in uniq):
            uniq.append(chi)
    if len(uniq) != r:
        raise RuntimeError(
            f"character table extraction found {len(uniq)} irreps, expected {r}")

    chars = np.array(uniq)
    # sort: trivial rep first, then by dimension
    e_class = class_of[e_idx]
    order_idx = sorted(range(r), key=lambda i: (
        0 if np.allclose(chars[i], 1.0, atol=1e-4) else 1,
        round(chars[i, e_class].real, 3),
    ))
    chars = chars[order_idx]

    chi_Q = np.array([g.chi_Q(g.elements[reps[k]]) for k in range(r)])

    return CharacterTable(
        group=g,
        classes=classes,
        class_reps=reps,
        class_sizes=sizes,
        class_orders=orders,
        inverse_class=inv_class,
        chars=chars,
        chi_Q=chi_Q,
    )
