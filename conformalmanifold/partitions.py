"""Partition combinatorics for nilpotent orbits and magnetic quivers.

Nilpotent orbits of the classical Lie algebras are labeled by partitions:

  * ``su``  (A-type, sl(N))   -- any partition of N.
  * ``so``  (B/D-type)        -- partition of N; even parts have even multiplicity.
  * ``sp``  (C-type, sp(2n))  -- partition of 2n; odd parts have even multiplicity.

The closure order of nilpotent orbits is the **dominance order** on partitions, and
the covering relations are the elementary degenerations that, in the magnetic-quiver
program, become the edges of the Higgs-branch Hasse diagram (Kraft-Procesi transitions).

The self-contained partition primitives (``transpose``, ``partitions_of``,
``valid_partition``, ``orbit_dim``, ``multiplicities``) are adapted from the
free-fermion-heterotic ``nilpotent/classical.py`` reference; the poset machinery
(``dominates``, ``covers``, ``covering_relations``, ``kp_slice_type``) is new here.

Convention checks locked by the tests:
  * dim_C O_lambda(su) = N^2 - sum_i (lambda^T_i)^2.
  * Kraft-Procesi minimal degenerations in type A are either a Kleinian slice
    ``A_k`` (= C^2/Z_{k+1}, quaternionic dim 1) or a minimal-orbit slice ``a_k``
    (= closure of the minimal nilpotent orbit of sl(k+1), quaternionic dim k).
    The transverse slice to the subregular orbit of sl(N) is C^2/Z_N = A_{N-1}
    (Brieskorn-Slodowy), so e.g. sl(3): (3) -> (2,1) is the slice A_2, and
    (2,1) -> (1^3) is the slice a_2.
"""

from __future__ import annotations

from typing import Iterator


Partition = tuple  # a non-increasing tuple of positive ints


# ---------------------------------------------------------------------------
# basic partition primitives
# ---------------------------------------------------------------------------

def normalize(part) -> Partition:
    """Sort descending and drop zero parts."""
    return tuple(sorted((int(p) for p in part if int(p) > 0), reverse=True))


def partitions_of(n: int, cap: int | None = None) -> Iterator[Partition]:
    """Generate every partition of ``n`` as a non-increasing tuple."""
    if n < 0:
        return
    if n == 0:
        yield ()
        return
    cap = n if cap is None else min(cap, n)
    for first in range(cap, 0, -1):
        for rest in partitions_of(n - first, first):
            yield (first,) + rest


def multiplicities(part) -> dict:
    """Map each part value to how many times it occurs."""
    m: dict = {}
    for k in part:
        m[k] = m.get(k, 0) + 1
    return m


def transpose(part) -> Partition:
    """The transpose (conjugate/dual) partition."""
    part = normalize(part)
    if not part:
        return ()
    return tuple(sum(1 for k in part if k > i) for i in range(part[0]))


def matrix_size(kind: str, part) -> int:
    """Size of the defining representation the partition lives in (= sum of parts)."""
    return sum(part)


def valid_partition(kind: str, part) -> bool:
    """Is ``part`` a valid nilpotent-orbit label for the given classical algebra?

    ``kind`` is ``'su'``/``'A'``, ``'so'``/``'B'``/``'D'`` or ``'sp'``/``'C'``.
    """
    part = normalize(part)
    mult = multiplicities(part)
    k = _kind_key(kind)
    if k == "su":
        return True
    if k == "so":  # even parts have even multiplicity
        return all(mult[v] % 2 == 0 for v in mult if v % 2 == 0)
    if k == "sp":  # odd parts have even multiplicity
        return all(mult[v] % 2 == 0 for v in mult if v % 2 == 1)
    raise ValueError(f"unknown algebra kind {kind!r}")


def _kind_key(kind: str) -> str:
    """Fold the family letters A/B/C/D onto the su/so/sp validity classes."""
    k = kind.lower()
    if k in ("su", "a", "gl", "u"):
        return "su"
    if k in ("so", "b", "d", "o"):
        return "so"
    if k in ("sp", "c", "usp"):
        return "sp"
    raise ValueError(f"unknown algebra kind {kind!r}")


def orbits(kind: str, N: int) -> list:
    """Every valid nilpotent-orbit partition of the size-``N`` defining rep."""
    return [p for p in partitions_of(N) if valid_partition(kind, p)]


def orbit_label(part) -> str:
    """Compact exponent notation, e.g. ``[3,2^2,1]``."""
    mult = multiplicities(normalize(part))
    pieces = [f"{k}^{mult[k]}" if mult[k] > 1 else f"{k}" for k in sorted(mult, reverse=True)]
    return "[" + ",".join(pieces) + "]"


# ---------------------------------------------------------------------------
# orbit dimension
# ---------------------------------------------------------------------------

def orbit_dim(kind: str, part) -> int:
    """Complex dimension of the nilpotent orbit O_part (Collingwood-McGovern 6.1.3-4).

    su:  N^2 - sum_i (t_i)^2                         with t = transpose(part).
    so:  (N^2 - N)/2 - (sum t_i^2 - #odd parts)/2.
    sp:  n(2n+1)   - (sum t_i^2 + #odd parts)/2      with N = 2n.
    """
    part = normalize(part)
    N = matrix_size(kind, part)
    t = transpose(part)
    s = sum(x * x for x in t)
    odd = sum(1 for k in part if k % 2 == 1)
    k = _kind_key(kind)
    if k == "su":
        return N * N - s
    if k == "so":
        return (N * N - N) // 2 - (s - odd) // 2
    if k == "sp":
        n = N // 2
        return n * (2 * n + 1) - (s + odd) // 2
    raise ValueError(f"unknown algebra kind {kind!r}")


def orbit_dim_quaternionic(kind: str, part) -> int:
    """Quaternionic (half-)dimension of the orbit closure = magnetic-quiver Coulomb dim."""
    d = orbit_dim(kind, part)
    if d % 2:
        raise ValueError(f"orbit dimension {d} is odd for {kind} {part!r}")
    return d // 2


# ---------------------------------------------------------------------------
# dominance order and covering relations (the orbit-closure poset)
# ---------------------------------------------------------------------------

def _partial_sums(part, length: int) -> list:
    out, acc = [], 0
    for i in range(length):
        acc += part[i] if i < len(part) else 0
        out.append(acc)
    return out


def dominates(a, b) -> bool:
    """``a >= b`` in dominance order: sum(a[:k]) >= sum(b[:k]) for all k (same total)."""
    a, b = normalize(a), normalize(b)
    if sum(a) != sum(b):
        raise ValueError("dominance order compares partitions of the same integer")
    length = max(len(a), len(b))
    pa, pb = _partial_sums(a, length), _partial_sums(b, length)
    return all(x >= y for x, y in zip(pa, pb))


def strictly_dominates(a, b) -> bool:
    a, b = normalize(a), normalize(b)
    return a != b and dominates(a, b)


def covers(a, b, universe=None) -> bool:
    """Does ``a`` cover ``b`` (a > b with nothing strictly between) in the poset?

    ``universe`` restricts the poset to a set of valid partitions (e.g. all so-orbits);
    if ``None`` the full partition dominance lattice is used.
    """
    a, b = normalize(a), normalize(b)
    if not strictly_dominates(a, b):
        return False
    if universe is None:
        universe = list(partitions_of(sum(a)))
    for c in universe:
        c = normalize(c)
        if c == a or c == b:
            continue
        if strictly_dominates(a, c) and strictly_dominates(c, b):
            return False
    return True


def covering_relations(kind: str, N: int) -> list:
    """All covering pairs (a, b) with a ⋗ b among the valid orbits of the type."""
    universe = orbits(kind, N)
    rels = []
    for a in universe:
        for b in universe:
            if a != b and covers(a, b, universe):
                rels.append((a, b))
    return rels


# ---------------------------------------------------------------------------
# Kraft-Procesi transition type of a covering (type A / su)
# ---------------------------------------------------------------------------

def _changed_rows(lam, mu) -> tuple:
    """Indices where two partitions differ; for a single-box move this is (i, j)."""
    length = max(len(lam), len(mu))
    lam = _pad(lam, length)
    mu = _pad(mu, length)
    diff = [r for r in range(length) if lam[r] != mu[r]]
    return diff


def _pad(part, length: int) -> tuple:
    part = normalize(part)
    return tuple(part[i] if i < len(part) else 0 for i in range(length))


def kp_slice_type(lam, mu) -> dict:
    """Kraft-Procesi transverse-slice type of a type-A covering ``lam ⋗ mu``.

    Returns ``{'letter': 'a'|'A', 'k': int, 'dim_H': int, 'label': str}``.

    * A single box moves from row ``i`` to row ``j`` (i<j).  If the move is between
      **adjacent rows** (j = i+1) the slice is the Kleinian ``A_k = C^2/Z_{k+1}``
      (quaternionic dim 1), with ``k`` read off the transpose column gap.
    * Otherwise the move is adjacent in the transpose (columns), giving the minimal
      slice ``a_k`` (closure of the minimal orbit of sl(k+1), quaternionic dim k).

    (``a_1`` and ``A_1`` coincide as ``C^2/Z_2``.)
    """
    lam, mu = normalize(lam), normalize(mu)
    if not covers(lam, mu):
        raise ValueError(f"{lam!r} does not cover {mu!r} in dominance order")
    rows = _changed_rows(lam, mu)
    trows = _changed_rows(transpose(lam), transpose(mu))
    i, j = rows[0], rows[-1]
    ip, jp = trows[0], trows[-1]
    dim_h = (orbit_dim("su", lam) - orbit_dim("su", mu)) // 2
    if (j - i) == 1:  # adjacent rows -> Kleinian slice
        k = jp - ip
        return {"letter": "A", "k": k, "dim_H": 1, "label": f"A{_sub(k)}"}
    # adjacent in the transpose -> minimal (row) slice
    k = j - i
    return {"letter": "a", "k": k, "dim_H": dim_h, "label": f"a{_sub(k)}"}


def _sub(k: int) -> str:
    """Render an integer subscript with Unicode subscript digits."""
    subs = "₀₁₂₃₄₅₆₇₈₉"
    return "".join(subs[int(d)] for d in str(k))
