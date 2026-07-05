"""5d SCFT / AdS6 readings of a toric Calabi-Yau diagram.

The same toric polygon that `toric.py` reads as a 4d N=1 quiver (D3-branes at the
singularity) also engineers a **5d SCFT** via M-theory on the CY3 / the IIB (p,q)
5-brane web.  This module collects the 5d-side data (the AdS6/CFT5 beat), kept
separate from the 4d quiver reading so it can grow into its own home (defect
group + pairing, SymTFT, rank/flavor, S^5 free energy, ...).

So far:

* **rank** (Coulomb-branch dimension) of the 5d SCFT = number of *interior*
  lattice points of the toric diagram (= number of compact divisors in a
  crepant resolution).  Computed in `api` straight from the diagram signature.

* **1-form symmetry** -- the 1-form part of the defect group, and the direct
  input to the 6d SymTFT that the AdS6 dual encodes.  For a toric CY3 it is the
  torsion group

      Z^2 / < edge vectors between consecutive polygon corners >,

  i.e. the cokernel of the 2 x V integer matrix whose columns are the
  corner-to-corner edge vectors of the (CCW) toric polygon.  Crucially the
  *corners* are used, so an edge of lattice length L contributes a length-L
  vector -- that is what distinguishes, e.g., C^3/(Z_n x Z_n) (the size-n
  triangle -> Z_n x Z_n) from C^3 (the unit triangle -> trivial).

  Validated against standard results: C^3 and the conifold -> trivial;
  dP0 = local P^2 (the E_0 SCFT) -> Z_3; C^3/(Z_n x Z_m) -> Z_n x Z_m;
  C^2/Z_n x C -> Z_n; F0 = local P^1 x P^1 -> Z_2; dP1, dP3 -> trivial.

Scope: this is the (electric) **1-form** symmetry.  The full defect group also
carries a magnetic part and a linking pairing -- the actual seed of the SymTFT
action -- which is the natural next addition on top of this.
"""

from __future__ import annotations

from functools import reduce
from math import gcd


def one_form_symmetry(hull) -> list[int]:
    """Invariant factors (each > 1) of the 5d SCFT's 1-form symmetry, read from
    the toric diagram's corner polygon `hull` (CCW lattice vertices, as returned
    by `toric.convex_hull`).

    Returns the list of invariant factors of  Z^2 / <edge vectors>:
        []      -> trivial
        [3]     -> Z_3
        [2, 2]  -> Z_2 x Z_2
    """
    pts = [(int(x), int(y)) for (x, y) in hull]
    V = len(pts)
    if V < 3:
        return []

    edges = [(pts[(i + 1) % V][0] - pts[i][0],
              pts[(i + 1) % V][1] - pts[i][1]) for i in range(V)]

    # Invariant factors d1 | d2 of the 2 x V integer matrix of edge vectors:
    #   d1       = gcd of all entries,
    #   d1 * d2  = gcd of all 2x2 minors,
    # and  Z^2 / <columns>  =  Z_d1 (+) Z_d2.
    d1 = reduce(gcd, (abs(c) for e in edges for c in e), 0)
    if d1 == 0:                       # degenerate (all corners coincide)
        return []
    minors = (abs(edges[i][0] * edges[j][1] - edges[i][1] * edges[j][0])
              for i in range(V) for j in range(i + 1, V))
    big = reduce(gcd, minors, 0)
    d2 = big // d1
    return [d for d in (d1, d2) if d > 1]


def non_isolated_note(edge_lengths) -> str:
    """Caveat for NON-isolated toric singularities, '' when isolated.

    A polygon edge of lattice length > 1 (extra lattice points on the edge)
    means a line of A-type singularities running off to infinity, carrying a 7d
    "flavor" gauge sector.  The quantities reported here are still the standard
    link/toric computations, but the theory is then not an *isolated*
    interacting 5d SCFT -- part of the data (e.g. some of the flavor rank, and
    center charges under the 1-form symmetry) belongs to that non-compact
    flavor sector.  Mirrors the non-isolated caveat on the 4d side (toric.py).
    """
    if any(int(l) > 1 for l in edge_lengths):
        return ("non-isolated singularity (a boundary edge has interior lattice "
                "points): a non-compact A-type singular line carries a 7d flavor "
                "sector, so these are link readings of a non-isolated geometry, "
                "not of an isolated interacting 5d SCFT")
    return ""


def flavor_rank(boundary_points: int) -> int:
    """Rank of the 5d SCFT's flavor symmetry (= number of mass parameters):

        flavor rank = (# boundary lattice points = # external (p,q) legs) - 3.

    Examples: C^3 and dP0 = local P^2 (B=3) -> 0 (E_0 / free: no flavor);
    local dP_n (B = n+3) -> n (the E_n flavor symmetry has rank n);
    SU(2) with N_f flavors -> N_f + 1 (the E_{N_f+1} symmetry).  This is a *rank*
    count -- robust even where the precise flavor group is subtle/enhanced.
    """
    return max(boundary_points - 3, 0)


def abelian_label(factors: list[int]) -> str:
    """Plain-text label for a finite abelian group given its invariant factors:
    [] -> 'trivial', [3] -> 'Z_3', [2, 2] -> 'Z_2 x Z_2'."""
    if not factors:
        return "trivial"
    return " x ".join(f"Z_{d}" for d in factors)
