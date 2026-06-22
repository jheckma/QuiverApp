"""Superconformal (SCFT) data for the quiver gauge theories.

For N D3-branes probing a Calabi-Yau three-fold singularity the low-energy
theory is a 4D N=1 superconformal field theory.  Where `conformal.py` reports
the *dimension of the conformal manifold*, this module reports the rest of the
protected superconformal data:

  * the gauge group  prod_i U(N_i),  N_i = N * dim R_i,  and the bifundamental
    matter content;
  * the exact superconformal R-charges of the fields (and of the superpotential,
    which is marginal: R_W = 2);
  * the 't Hooft anomalies  Tr R,  Tr R^3,  hence the central charges
    (Anselmi-Freedman-Grisaru-Johansen)

        a = (3/32)(3 Tr R^3 - Tr R),     c = (1/32)(9 Tr R^3 - 5 Tr R),

    and the holographic statement  a = c  at leading order in N.

Two routes that agree on their overlap (the C^3/Gamma orbifold points):

  ORBIFOLD  C^3/Gamma  -- exact and closed form.  The theory is the Gamma-
    projection of N=4  U(N|Gamma|...)  SYM, so every bifundamental inherits the
    N=4 value  R = 2/3  and the cubic superpotential  W = Tr Phi^1[Phi^2,Phi^3]
    has  R = 2.  The McKay relation  sum_j a_ij d_j = 3 d_i  (because the
    defining rep Q is 3-dimensional) forces  Tr R = 0, so  a = c  identically,
    and

        a = c = (1/4)(sum_i d_i^2) N^2 = (|Gamma|/4) N^2,

    i.e. |Gamma| times the N=4 value -- the field-theory shadow of
    Vol(S^5/Gamma) = Vol(S^5)/|Gamma| in the dual AdS_5 x S^5/Gamma.

  TORIC CY3  -- Martelli-Sparks-Yau volume (Z-) minimisation over the Reeb
    vector (hep-th/0503183, hep-th/0503184).  At leading order in N

        a = c = pi^3 N^2 / (4 Vol(X_5)),

    with  Vol(X_5)  the minimised Sasaki-Einstein volume read off the toric
    diagram, and the per-corner R-charges  R_a = 2 t_a / sum_b t_b  (so
    sum_a R_a = 2).  Reproduces  a = N^2/4  for S^5 (C^3),  27 N^2/64  for the
    conifold, and  |Gamma| N^2/4  at the C^3/Gamma orbifold points -- matching
    the orbifold route above.

Only dependency is the standard library plus numpy (already required).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from .groups import MatrixGroup
from .quiver import McKayQuiver


# ===========================================================================
# Central charges from 't Hooft anomalies (shared helper)
# ===========================================================================
def central_charges(tr_R, tr_R3):
    """a, c from the R-symmetry 't Hooft anomalies (AFGJ).

        a = (3/32)(3 Tr R^3 - Tr R),   c = (1/32)(9 Tr R^3 - 5 Tr R).

    `tr_R`, `tr_R3` may be Fractions/ints (exact) or floats; the return type
    follows the input.
    """
    if isinstance(tr_R, Fraction) or isinstance(tr_R3, Fraction):
        a = Fraction(3, 32) * (3 * tr_R3 - tr_R)
        c = Fraction(1, 32) * (9 * tr_R3 - 5 * tr_R)
    else:
        a = (3.0 / 32.0) * (3 * tr_R3 - tr_R)
        c = (1.0 / 32.0) * (9 * tr_R3 - 5 * tr_R)
    return a, c


def _frac_json(x: Fraction) -> dict:
    """JSON-friendly view of an exact rational coefficient of N^2."""
    return {"frac": (str(x) if x.denominator != 1 else str(x.numerator)),
            "val": float(x)}


# ===========================================================================
# Orbifold  C^3/Gamma  -- exact superconformal data on the McKay quiver
# ===========================================================================
@dataclass
class OrbifoldSCFT:
    gauge_group: str          # "U(N) x U(2N) x ..."
    ranks: list[int]          # the d_i (N_i = d_i * N)
    num_bifund: int           # total bifundamental chiral multiplets
    R_chiral: Fraction        # superconformal R-charge of every bifundamental
    R_superpotential: int     # = 2 (marginal cubic W)
    tr_R: Fraction            # coefficient of N^2 in Tr R
    tr_R3: Fraction           # coefficient of N^2 in Tr R^3
    a_coeff: Fraction         # a = a_coeff * N^2
    c_coeff: Fraction         # c = c_coeff * N^2
    a_eq_c: bool
    a_over_aN4: Fraction      # a / a_{N=4} = |Gamma|   (faithful)
    note: str = ""


def orbifold_scft(group: MatrixGroup, quiver: McKayQuiver) -> OrbifoldSCFT:
    """Exact superconformal data of the C^3/Gamma orbifold SCFT.

    Anomalies are assembled from the McKay quiver directly (gauginos R = +1,
    bifundamental fermions R = R_chiral - 1 = -1/3) and returned as exact
    rational coefficients of N^2, with N_i = N * d_i.
    """
    dims = [int(d) for d in quiver.dims]
    A = quiver.adjacency
    n = len(dims)
    adj = [[int(round(float(A[i][j].real if np.iscomplexobj(A) else A[i][j])))
            for j in range(n)] for i in range(n)]

    R_chi = Fraction(2, 3)              # inherited N=4 value
    R_ferm = R_chi - 1                  # = -1/3, the bifundamental fermion

    # P = sum_i d_i^2 = |Gamma|;   Qd = sum_ij a_ij d_i d_j  (= 3 P, faithful Q)
    P = sum(d * d for d in dims)
    Qd = sum(adj[i][j] * dims[i] * dims[j] for i in range(n) for j in range(n))

    # Tr over Weyl fermions, as a coefficient of N^2:
    #   gauginos     -> sum_i N_i^2 (R=+1)            ->  P
    #   bifund ferms -> sum_ij a_ij N_i N_j (R_ferm)  ->  Qd * R_ferm
    tr_R = Fraction(P) + Qd * R_ferm
    tr_R3 = Fraction(P) + Qd * (R_ferm ** 3)
    a_coeff, c_coeff = central_charges(tr_R, tr_R3)

    labels = [quiver.node_rank_label(i) for i in range(n)]
    gauge = " × ".join(labels)

    note = ("Leading order in N (planar): gauge group taken as prod U(N_i); the "
            "overall and relative U(1)s are IR-free / become global, an O(1) "
            "correction. a = c is exact (Tr R = 0).")
    if not a_eq_c_flag(a_coeff, c_coeff):
        note = ("Tr R != 0 for this action (Q not a faithful 3d rep with "
                "chi_Q(1)=3); a and c differ -- reported as computed.")

    return OrbifoldSCFT(
        gauge_group=gauge,
        ranks=dims,
        num_bifund=quiver.num_arrows,
        R_chiral=R_chi,
        R_superpotential=2,
        tr_R=tr_R,
        tr_R3=tr_R3,
        a_coeff=a_coeff,
        c_coeff=c_coeff,
        a_eq_c=a_eq_c_flag(a_coeff, c_coeff),
        a_over_aN4=a_coeff * 4,
        note=note,
    )


def a_eq_c_flag(a_coeff, c_coeff) -> bool:
    return a_coeff == c_coeff


def orbifold_scft_json(group: MatrixGroup, quiver: McKayQuiver) -> dict:
    s = orbifold_scft(group, quiver)
    return {
        "method": "orbifold (exact): N=4 R-charges + 't Hooft anomalies",
        "gauge_group": s.gauge_group,
        "ranks": s.ranks,
        "num_bifundamentals": s.num_bifund,
        "R_chiral": {"frac": str(s.R_chiral), "val": float(s.R_chiral)},
        "R_superpotential": s.R_superpotential,
        "tr_R": _frac_json(s.tr_R),
        "tr_R3": _frac_json(s.tr_R3),
        "tr_R_zero": s.tr_R == 0,
        "a": _frac_json(s.a_coeff),       # coefficient of N^2
        "c": _frac_json(s.c_coeff),       # coefficient of N^2
        "a_eq_c": s.a_eq_c,
        "a_over_aN4": _frac_json(s.a_over_aN4),
        "note": s.note,
    }


# ===========================================================================
# Toric CY3 -- Martelli-Sparks-Yau Sasaki-Einstein volume minimisation
# ===========================================================================
def _det3(u, v, w) -> float:
    return float(np.linalg.det(np.array([u, v, w], dtype=float)))


def _cross(u, v):
    return np.cross(np.asarray(u, float), np.asarray(v, float))


def minimize_volume(corners):
    """Z-minimise the toric Sasaki-Einstein volume over the Reeb vector.

    `corners` : CCW lattice-polygon vertices (the corners of the toric diagram;
                collinear boundary points dropped) -- the extremal rays of the
                Calabi-Yau cone, lifted to V_a = (1, x_a, y_a).

    Returns (g_min, reeb, t, converged) where, with b = (3, b2, b3),

        g(b) = sum_a  n_a / (D1_a(b) D2_a(b)),
        n_a  = det(V_{a-1}, V_a, V_{a+1}),
        D1_a = (V_{a-1} x V_a) . b,   D2_a = (V_a x V_{a+1}) . b,

    minimised over (b2, b3) (b1 = 3 fixed by the CY condition).  The minimised
    Sasaki-Einstein volume is Vol(X_5) = (pi^3 / 3) g_min, and `t` are the
    per-corner contributions t_a = n_a/(D1_a D2_a) at the minimum (R_a = 2 t_a /
    g_min, sum_a R_a = 2).  The function is strictly convex on the Reeb cone
    (MSY), so damped Newton from the lattice centroid converges.
    """
    pts = [(float(x), float(y)) for (x, y) in corners]
    D = len(pts)
    if D < 3:
        raise ValueError("need at least a triangle (3 corners) for a CY3 cone")
    V = [(1.0, x, y) for (x, y) in pts]
    n = [_det3(V[a - 1], V[a], V[(a + 1) % D]) for a in range(D)]
    c1 = [_cross(V[a - 1], V[a]) for a in range(D)]          # for D1_a
    c2 = [_cross(V[a], V[(a + 1) % D]) for a in range(D)]    # for D2_a

    cx = sum(p[0] for p in pts) / D
    cy = sum(p[1] for p in pts) / D
    b = np.array([3.0, 3.0 * cx, 3.0 * cy])                  # symmetric start

    def denoms(bvec):
        p = np.array([c1[a] @ bvec for a in range(D)])
        q = np.array([c2[a] @ bvec for a in range(D)])
        return p, q

    converged = False
    for _ in range(200):
        p, q = denoms(b)
        if np.any(p <= 0) or np.any(q <= 0):
            # outside the Reeb cone (shouldn't happen from the centroid start);
            # nudge back toward the centroid and retry.
            b[1:] = 0.5 * b[1:] + 0.5 * np.array([3.0 * cx, 3.0 * cy])
            continue
        t = np.array([n[a] / (p[a] * q[a]) for a in range(D)])
        grad = np.zeros(2)
        H = np.zeros((2, 2))
        for a in range(D):
            pk = np.array([c1[a][1], c1[a][2]])
            qk = np.array([c2[a][1], c2[a][2]])
            u = pk / p[a] + qk / q[a]
            grad += -t[a] * u
            H += t[a] * (np.outer(u, u)
                         + np.outer(pk, pk) / p[a] ** 2
                         + np.outer(qk, qk) / q[a] ** 2)
        if np.linalg.norm(grad) < 1e-13:
            converged = True
            break
        step = np.linalg.solve(H, grad)
        # damped Newton: backtrack so the step stays inside the Reeb cone
        s = 1.0
        for _ls in range(40):
            nb = b.copy()
            nb[1:] = b[1:] - s * step
            pp, qq = denoms(nb)
            if np.all(pp > 0) and np.all(qq > 0):
                b = nb
                break
            s *= 0.5
        else:
            converged = True            # cannot move; treat as a minimum
            break

    p, q = denoms(b)
    t = np.array([n[a] / (p[a] * q[a]) for a in range(D)])
    g_min = float(t.sum())
    return g_min, b, t, converged


def toric_field_R_charges(corners, fields, superpotential, corner_R) -> dict:
    """Per-field superconformal R-charges for a toric quiver (Butti-Zaffaroni).

    Each field sits at a crossing of two zig-zag legs; its R-charge is the sum of
    the per-corner MSY R-charges over the corners in the arc between those legs.
    The arc (one of two, summing to 2) is fixed by requiring every superpotential
    term to be marginal (sum R = 2); the assignment is then independently checked
    against the gauge-node NSVZ condition (sum_{fields at node}(1 - R) = 2).

    `corners`          : CCW polygon vertices (the toric-diagram corners).
    `fields`           : the inverse-algorithm fields, each with "label",
                         "zigzag" = [leg_k, leg_l], "src", "tgt".
    `superpotential`   : the two-term W, each {"fields": [labels...]}.
    `corner_R`         : per-corner R-charges [{"corner": [x, y], "R": val}, ...].
    """
    from math import gcd

    from .inverse import zigzag_windings

    hull = [(int(x), int(y)) for (x, y) in corners]
    n = len(hull)
    Rof = {tuple(c["corner"]): c["R"] for c in corner_R}
    # which polygon edge each zig-zag leg belongs to (CCW, with multiplicity g)
    leg_edge = []
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        g = gcd(abs(b[0] - a[0]), abs(b[1] - a[1])) or 1
        leg_edge += [i] * g
    B = len(leg_edge)
    # gap[j] sits between leg j and leg j+1: the shared corner's R if the edge
    # index changes there, else 0.  (sum of gaps = sum of corner R = 2)
    gap = [Rof[hull[(leg_edge[(j + 1) % B])]] if leg_edge[j] != leg_edge[(j + 1) % B]
           else 0.0 for j in range(B)]

    def arc(k, l):
        s, j = 0.0, k
        while j != l:
            s += gap[j]
            j = (j + 1) % B
        return s

    # candidate arc-sums; pick the minimal arc (verified below to be marginal)
    R = {}
    for f in fields:
        r = arc(f["zigzag"][0], f["zigzag"][1])
        R[f["label"]] = min(r, 2.0 - r)

    # certificate 1: every superpotential term is marginal (sum R = 2)
    w_marginal = all(abs(sum(R[x] for x in term["fields"]) - 2.0) < 1e-6
                     for term in superpotential)
    # certificate 2: gauge-node NSVZ condition  sum_{fields at node}(1 - R) = 2
    node = {}
    for f in fields:
        for gnode in (f["src"], f["tgt"]):
            node[gnode] = node.get(gnode, 0.0) + (1.0 - R[f["label"]])
    nsvz_ok = all(abs(v - 2.0) < 1e-6 for v in node.values())

    field_R = [{"label": f["label"], "src": f["src"], "tgt": f["tgt"],
                "R": round(R[f["label"]], 5)} for f in fields]
    # distinct values, clustering numeric-minimisation noise (tol 1e-3)
    distinct = []
    for v in sorted(R[f["label"]] for f in fields):
        if not distinct or abs(v - distinct[-1]) > 1e-3:
            distinct.append(round(v, 4))
    return {
        "field_R": field_R,
        "distinct_R": distinct,
        "W_marginal": bool(w_marginal),     # every term sum R = 2
        "nsvz_ok": bool(nsvz_ok),           # every gauge node sum(1-R) = 2
    }


def toric_scft_json(corners) -> dict:
    """Superconformal data of the toric CY3 SCFT, from its toric diagram.

    `corners` : CCW lattice-polygon vertices.  Returns central charge a = c (as
    a coefficient of N^2), the minimising Reeb vector and SE volume ratio, and
    the per-corner superconformal R-charges (sum = 2)."""
    g, b, t, ok = minimize_volume(corners)
    a_over_n2 = 3.0 / (4.0 * g)                 # = pi^3 / (4 Vol);  S^5 -> 1/4
    R = [2.0 * ta / g for ta in t]
    pts = [(int(round(x)), int(round(y))) for (x, y) in corners]
    return {
        "method": "Martelli–Sparks–Yau Sasaki–Einstein volume "
                  "minimisation",
        "a": {"val": a_over_n2},                # coefficient of N^2
        "c": {"val": a_over_n2},                # a = c at leading order in N
        "a_eq_c": True,
        "a_over_aN4": 3.0 / g,                  # = Vol(S^5)/Vol(X_5)
        "reeb": [round(float(x), 5) for x in b],
        "vol_over_volS5": g / 3.0,              # Vol(X_5)/Vol(S^5)
        "corner_R": [{"corner": list(pts[a]), "R": round(R[a], 5)}
                     for a in range(len(pts))],
        "R_sum": round(sum(R), 5),
        "converged": bool(ok),
        "note": "Leading order in N: a = c for any toric Sasaki–Einstein "
                "dual. R-charges shown are the extremal (per-corner) values; "
                "field R-charges are non-negative integer sums of these.",
    }
