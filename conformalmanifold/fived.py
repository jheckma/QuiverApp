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

* **defect group + pairing + global forms** (`defect_group`).  Lines (M2 on
  non-compact 2-cycles mod screening) and surface defects (M5 on non-compact
  4-cycles mod screening) organize into  D = Gamma_e (+) Gamma_m  with
  Gamma_m ~= Gamma_e = Gamma  by the perfectness of the torsion linking form on
  the link.  Following Albertini-Del Zotto-Garcia Etxebarria-Hosseini
  (arXiv:2005.12831): Gamma = Tor coker(Q2), with Q2 the intersection matrix of
  curve classes x compact divisors of a crepant resolution (their eq. 3.9c),
  and the linking pairing = "inverse of Q2 mod 1" on torsion generators (their
  eqs. 3.15-3.17).  Both are computed here WITHOUT choosing a triangulation:
  curve classes = integer relations among the rays (v_i, 1) over the lattice
  points (GLSM charges), so flop-invariance is manifest.  On canonically paired
  (Smith-normal-form) generators the pairing collapses to the hyperbolic form
  l = 1/d per factor; published generator conventions differ by a unit/sign
  (ADEH quote -1/3 for dP0, -1/gcd(p,q) for Y^{p,q}), which does not affect
  the global-form (polarization) enumeration.  Absolute 5d theories correspond
  to maximal isotropic ("Lagrangian") subgroups of D; for isolated toric
  singularities Gamma is always cyclic (ADEH eq. 3.24-3.25) and the number of
  global forms is sigma(d) = sum of divisors of d (dP0: 4, F0/SU(2)_0: 3).

  Validated against: dP0 = E0 (D = Z3 (+) Z3, ADEH eq. 6.10), Y^{p,q}
  (Gamma = Z_gcd(p,q), ADEH eq. 4.4), the corner-triangle gcd shortcut (their
  eqs. 3.24-3.25), and agreement with the boundary formula `one_form_symmetry`
  on all isolated cases (named + randomized sweeps).

Scope: for NON-isolated singularities the field-theory defect group can be
*smaller* than the link reading (screening by the 7d flavor sector on the
singular line); ADEH's modified prescription (their sec. 5.2) is conjectural,
so here we keep reporting the link group with a note and defer the screened
group.  The cubic (anomaly) SymTFT term (Apruzzi et al, arXiv:2112.02092,
eq. 5.24) is extracted but not yet implemented.
"""

from __future__ import annotations

from fractions import Fraction
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


def _prim_len(vx, vy):
    """(primitive vector, lattice length) of an integer vector."""
    g = gcd(abs(int(vx)), abs(int(vy)))
    if g == 0:
        return (0, 0), 0
    return (vx // g, vy // g), g


def gravity_dual(hull) -> dict:
    """The Type IIB supergravity dual of the 5d SCFT engineered by this toric
    diagram's (p,q) 5-brane web, in the **D'Hoker-Gutperle-Uhlemann** family:
    a warped

        ds^2 = f_6^2 ds^2_{AdS_6} + f_2^2 ds^2_{S^2} + 4 rho^2 |dz|^2

    over a Riemann surface Sigma (a disk / upper-half-plane, genus 0), with the
    entire solution fixed by two locally holomorphic functions

        A_+-(z) = A_+-^0 + sum_{l=1}^L  Z_+-^l  ln(z - r_l),   r_l in dSigma,

    whose logarithmic branch points r_l on the boundary are the external
    5-brane stacks of the web.  The residues Z_+-^l encode the (p_l, q_l)
    charges of the l-th external leg.

    Dictionary to the toric diagram (this is the exact, computable part):
    the external legs of the (p,q) web are the OUTWARD primitive normals of the
    toric polygon's edges; an edge of lattice length n_l is a stack of n_l
    (p_l, q_l) 5-branes, so the l-th DGU puncture carries total charge
    (p_l, q_l) = n_l * (outward primitive normal).  5-brane charge conservation
    sum_l (p_l, q_l) = 0 is automatic because the polygon closes.

    Refs: D'Hoker-Gutperle-Uhlemann arXiv:1611.09411, 1703.08186, 1706.00433;
    free energy: Fluder-Uhlemann arXiv:1806.08374, Kaidi-Uhlemann.

    Returns a descriptive dict (the ansatz, Sigma data, the external-leg pole
    table with charges, the holomorphic-function form, charge-conservation
    check, and the free-energy reading); NOT a full supergravity solver."""
    from .toric import convex_hull

    hull = convex_hull([(int(x), int(y)) for (x, y) in hull])
    V = len(hull)
    poles = []
    sx = sy = 0
    for i in range(V):
        (ax, ay), (bx, by) = hull[i], hull[(i + 1) % V]
        # CCW hull -> outward normal is the RIGHT normal (e_y, -e_x)
        prim, n = _prim_len(by - ay, -(bx - ax))
        p, q = prim[0] * n, prim[1] * n
        sx += p
        sy += q
        poles.append({
            "index": i + 1,
            "pq": [p, q],                    # total charge of the stack
            "brane_type": list(prim),        # the primitive (p,q) 5-brane
            "num_branes": n,                 # # 5-branes in the stack = edge length
            "edge": [list(hull[i]), list(hull[(i + 1) % V])],
            "position": f"r_{i + 1} on dSigma",
            # residue dictionary (Fluder-Uhlemann arXiv:2011.00006):
            #   Z_+-^l = (3/4) alpha' (+- q_l + i p_l)
            "residue": f"Z_pm^{i+1} = (3/4)a'(±{q} + i·{p})",
        })
    L = V
    charge_conserved = (sx, sy) == (0, 0)
    # holomorphic functions, with residues shown as the (p,q) charges
    a_terms = " + ".join(f"Z_pm^{p['index']} ln(z - r_{p['index']})" for p in poles)
    return {
        "family": "D'Hoker-Gutperle-Uhlemann warped AdS_6 x S^2 (Type IIB)",
        "ansatz": "ds^2 = f_6^2 ds^2_{AdS_6} + f_2^2 ds^2_{S^2} + 4 rho^2 |dz|^2",
        "riemann_surface": {
            "type": "disk (upper half-plane)",
            "genus": 0,
            "num_punctures": L,
            "boundary": "the L external 5-brane stacks sit at r_1..r_L on dSigma",
        },
        "num_external_stacks": L,
        "five_brane_poles": poles,
        "holomorphic_functions": f"A_pm(z) = A_pm^0 + {a_terms}",
        "residue_dictionary": ("Z_pm^l fixes the l-th external leg's (p_l,q_l); "
                               "(p_l,q_l) = (edge lattice length) x (outward "
                               "primitive normal of the polygon edge)"),
        "charge_conservation": {
            "sum_pq": [sx, sy],
            "conserved": charge_conserved,
            "note": "sum_l (p_l,q_l)=0 (the polygon closes) -- required for a web",
        },
        "free_energy": _ads6_free_energy(hull, poles),
        "refs": ["arXiv:1611.09411", "arXiv:1703.08186", "arXiv:1706.00433",
                 "arXiv:1806.08374"],
    }


def _ads6_free_energy(hull, poles) -> dict:
    """The holographic S^5 free energy of the DGU dual, computed from the
    external-leg charges (Fluder-Uhlemann arXiv:1806.08374, 2011.00006).

    The general result is a double sum over the L poles,

        F_{S^5} = 64 pi^2 / (9 (2 pi a')^4) [ zeta(3) sum_{l,k} (Z^[lk])^2
                                              - (tri-log term over 4 distinct poles) ],

    Z^[lk] = Z_+^l Z_-^k - Z_+^k Z_-^l = -(9 i/8) a'^2 (p_l q_k - p_k q_l).
    The charge-only (zeta(3)) piece evaluates in closed form to

        F_{S^5} = -(9/8) zeta(3)/pi^2 * S,   S = sum_{l<k} (p_l q_k - p_k q_l)^2,

    which is EXACT for a 3-leg web (no 4-distinct-pole tri-log term) and the
    leading piece otherwise (the tri-log term depends on the boundary moduli
    r_l, fixed only by the global regularity analysis).  Verified: T_N (legs
    (N,0),(0,N),(-N,-N)) gives S = 3N^4 -> F = -27 zeta(3)/(8 pi^2) N^4, the
    published value.  Quartic in the charges: a web scaled by N scales as N^4
    (NOT N^{5/2}; that scaling is the massive-IIA Brandhuber-Oz dual, a
    different construction)."""
    from math import pi
    try:
        from math import zeta                      # py3.12+: not in math
    except ImportError:
        zeta3 = 1.2020569031595942
    else:                                           # pragma: no cover
        zeta3 = zeta(3)
    zeta3 = 1.2020569031595942
    L = len(poles)
    pq = [p["pq"] for p in poles]
    S = 0
    for a in range(L):
        for b in range(a + 1, L):
            det = pq[a][0] * pq[b][1] - pq[a][1] * pq[b][0]
            S += det * det
    coeff = -9 * zeta3 / (8 * pi ** 2)              # F = coeff * S  (units a'^0)
    F = coeff * S
    exact = (L == 3)
    return {
        "quantity": "F_{S^5}  (= -log Z_{S^5}, in units where 2 pi a' = 1)",
        "formula": "F_{S^5} = -(9/8) zeta(3)/pi^2 * sum_{l<k}(p_l q_k - p_k q_l)^2"
                   + ("" if exact else "  + tri-log(r_l) term (4-pole)"),
        "wedge_sum_S": S,                           # SL(2,Z)-invariant charge sum
        "value": round(F, 6),
        "exact": exact,                             # exact for 3-leg webs
        "scaling": "quartic in the 5-brane charges (a uniformly-N web ~ N^4; "
                   "e.g. T_N: F = -27 zeta(3)/(8 pi^2) N^4)",
        "note": ("" if exact else
                 "L>=4: the shown value is the charge-only zeta(3) piece; the "
                 "full F has an additional single-valued-trilog term depending "
                 "on the boundary pole positions r_l (regularity moduli)"),
        "refs": ["arXiv:1806.08374", "arXiv:2011.00006"],
    }


# ----------------------------------------------------------------------------
# integer linear algebra (pure stdlib): Smith normal form with transforms
# ----------------------------------------------------------------------------
def smith_normal_form(Ain):
    """Smith normal form over Z with transforms: returns (P, D, Q) with
    P @ A @ Q = D, P and Q unimodular, D diagonal with d1 | d2 | ...  Small
    dense matrices only (polygon-sized)."""
    A = [list(map(int, row)) for row in Ain]
    m, n = len(A), (len(A[0]) if A else 0)
    P = [[int(i == j) for j in range(m)] for i in range(m)]
    Q = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]
        P[i], P[j] = P[j], P[i]

    def swap_cols(i, j):
        for r in A: r[i], r[j] = r[j], r[i]
        for r in Q: r[i], r[j] = r[j], r[i]

    def add_row(src, dst, c):          # row_dst += c * row_src
        A[dst] = [a + c * b for a, b in zip(A[dst], A[src])]
        P[dst] = [a + c * b for a, b in zip(P[dst], P[src])]

    def add_col(src, dst, c):
        for r in A: r[dst] += c * r[src]
        for r in Q: r[dst] += c * r[src]

    t = 0
    while t < min(m, n):
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        swap_rows(t, piv[0]); swap_cols(t, piv[1])
        again = True
        while again:
            again = False
            for i in range(t + 1, m):
                if A[i][t]:
                    q = A[i][t] // A[t][t]
                    add_row(t, i, -q)
                    if A[i][t]:
                        swap_rows(t, i); again = True
            for j in range(t + 1, n):
                if A[t][j]:
                    q = A[t][j] // A[t][t]
                    add_col(t, j, -q)
                    if A[t][j]:
                        swap_cols(t, j); again = True
        ok = False
        while not ok:                   # enforce d_t | (everything below-right)
            ok = True
            for i in range(t + 1, m):
                bad = next((j for j in range(t + 1, n) if A[i][j] % A[t][t]), None)
                if bad is not None:
                    add_row(i, t, 1)
                    for jj in range(t + 1, n):
                        if A[t][jj]:
                            q = A[t][jj] // A[t][t]
                            add_col(t, jj, -q)
                            if A[t][jj]:
                                swap_cols(t, jj)
                    ok = False
                    break
        if A[t][t] < 0:
            A[t] = [-x for x in A[t]]
            P[t] = [-x for x in P[t]]
        t += 1
    return P, A, Q


def _integer_kernel(U):
    """Z-basis (as rows) of the saturated lattice {x in Z^n : U x = 0}."""
    P, D, Q = smith_normal_form(U)
    m, n = len(U), len(U[0])
    r = sum(1 for k in range(min(m, n)) if D[k][k])
    cols = list(zip(*Q))
    return [list(cols[j]) for j in range(r, n)]


def sigma_divisors(d: int) -> int:
    """sigma(d) = sum of divisors of d = # of maximal isotropic (Lagrangian)
    subgroups of Z_d (+) Z_d with the standard symplectic pairing = # of
    global forms of a 5d theory with cyclic defect-group factor Z_d."""
    return sum(k for k in range(1, d + 1) if d % k == 0)


# ----------------------------------------------------------------------------
# the full defect group (ADEH, arXiv:2005.12831 secs. 3.2-3.3)
# ----------------------------------------------------------------------------
def defect_group(hull) -> dict:
    """Full defect group D = Gamma_e (+) Gamma_m of the 5d theory, with the
    canonical Dirac/linking pairing and the global-form count.

    Returns a dict:
      isolated         bool -- polygon edges all primitive?
      factors          invariant factors of Gamma (electric = magnetic half)
      label            e.g. 'Z_3 (+) Z_3' for D, or 'trivial'
      pairing          ['1/3', ...] canonical hyperbolic coefficients (isolated
                       only; generator-convention note in the module docstring)
      num_global_forms sigma(d) for isolated (cyclic) Gamma; None otherwise
      note             '' or the honest-scope caveat
    """
    from .resolution import lattice_points, interior_lattice_points
    from .toric import convex_hull, polygon_signature

    hull = convex_hull([(int(x), int(y)) for (x, y) in hull])
    _, _, _, edges = polygon_signature(hull)
    isolated = all(int(e) == 1 for e in edges)
    link_factors = one_form_symmetry(hull)

    if not isolated:
        return {
            "isolated": False,
            "factors": link_factors,
            "label": _defect_label(link_factors),
            "pairing": None,
            "num_global_forms": None,
            "note": ("link reading only: for a non-isolated singularity the "
                     "field-theory defect group can be smaller (screening by "
                     "the 7d flavor sector on the singular line; ADEH sec. 5.2 "
                     "prescription not yet implemented)"),
        }

    # Route A (ADEH eq. 3.9c), triangulation-free:
    #   rays u_i = (v_i, 1) over ALL lattice points; curve classes = integer
    #   relations among the rays; Q2 = relations restricted to interior columns.
    lat = lattice_points(hull)
    interior = set(interior_lattice_points(hull))
    int_idx = [i for i, p in enumerate(lat) if p in interior]
    if int_idx:
        rays = [[p[0] for p in lat], [p[1] for p in lat], [1] * len(lat)]
        K = _integer_kernel(rays)                     # rows = curve classes
        A = [[row[i] for row in K] for i in int_idx]  # I x (n-3) map matrix
        _, D, _ = smith_normal_form(A)
        ds = [D[k][k] for k in range(min(len(A), len(A[0])))]
        factors = [d for d in ds if d > 1]
    else:
        factors = []                                  # I = 0 (isolated): trivial

    note = ""
    if sorted(factors) != sorted(link_factors):       # defensive; never seen
        note = (f"intersection-matrix group {abelian_label(factors)} disagrees "
                f"with the boundary reading {abelian_label(link_factors)} -- "
                "please report this diagram")

    return {
        "isolated": True,
        "factors": factors,
        "label": _defect_label(factors),
        "pairing": [str(Fraction(1, d)) for d in factors],
        "num_global_forms": _num_global_forms(factors),
        "note": note,
    }


NOT_COMPUTED = "not computed — no published closed form"


def cubic_anomaly(hull) -> dict:
    """Cubic 't Hooft anomaly of the 5d 1-form symmetry: the coefficient of
    int c_2^3 in the 6d SymTFT (the anomaly of Gukov-Hsin-Pei,
    arXiv:2010.15890).  It refines the global-form count: polarizations that
    would gauge an anomalous subgroup are obstructed.

    Only PUBLISHED closed forms are reported (Apruzzi-Bonetti-Garcia
    Etxebarria-Hosseini-Schafer-Nameki, arXiv:2112.02092); the general web
    evaluation needs case-by-case counterterm fixing (their sec. 5.3), so
    everything else honestly returns "not computed -- no published closed
    form" (NOT_COMPUTED above):

      * trivial 1-form symmetry            -> 0
      * su(p)_q webs (Y^{p,q} trapezoids)  -> q p (p-1)(p-2) / (6 gcd(p,q)^3)
                                              [eq. 5.42; su(2)_0 -> 0, eq. 5.33]
      * B_N triangles (B_3 = local P^2)    -> (N-1)(N-2) / (6 (N^2-3N+3))
                                              [eqs. 4.34 = 5.48, two routes]

    Conventions: coefficients as published, for the canonical orientation
    (0 <= q < p); an orientation flip (parity) reverses the sign.  Mixed
    anomalies with the instanton U(1) are not included.

    Returns {status: 'computed'|'not_computed', coefficient: str|None,
             family: str|None, note: str}.
    """
    from math import isqrt
    from .toric import convex_hull, polygon_signature, gl2z_equiv

    hull = convex_hull([(int(x), int(y)) for (x, y) in hull])
    area2, B, I, edge_lengths = polygon_signature(hull)
    not_computed = {"status": "not_computed", "coefficient": None,
                    "family": None, "note": NOT_COMPUTED}
    if any(int(e) != 1 for e in edge_lengths):        # non-isolated: out of scope
        return not_computed

    factors = one_form_symmetry(hull)
    if not factors:
        return {"status": "computed", "coefficient": "0",
                "family": "trivial 1-form symmetry", "note": ""}
    d = factors[0]                                     # isolated toric: cyclic

    # su(p)_q family: trapezoid {(-1,0),(0,0),(1,p-q),(0,p)}, 2*area = 2p
    if area2 % 2 == 0 and area2 >= 4:
        p = area2 // 2
        for q in range(p):
            g = gcd(p, q) if q else p                  # gcd(p, 0) = p
            if g != d:
                continue
            ref = convex_hull([(-1, 0), (0, 0), (1, p - q), (0, p)])
            if gl2z_equiv(hull, ref):
                coeff = Fraction(q * p * (p - 1) * (p - 2), 6 * d ** 3)
                return {"status": "computed", "coefficient": str(coeff),
                        "family": f"su({p})_{q} web  (Y^{{{p},{q}}})",
                        "note": ""}

    # B_N family: triangle {(N-1,0),(1,N-1),(0,1)}, 2*area = N^2-3N+3 = |Gamma|
    s = isqrt(4 * area2 - 3)
    if s * s == 4 * area2 - 3 and (3 + s) % 2 == 0:
        N = (3 + s) // 2
        if N >= 3 and d == area2:
            ref = convex_hull([(N - 1, 0), (1, N - 1), (0, 1)])
            if gl2z_equiv(hull, ref):
                coeff = Fraction((N - 1) * (N - 2), 6 * (N * N - 3 * N + 3))
                fam = f"B_{N} web" + ("  (local P^2, E_0)" if N == 3 else "")
                return {"status": "computed", "coefficient": str(coeff),
                        "family": fam, "note": ""}

    return not_computed


def _defect_label(factors: list[int]) -> str:
    """Label for D = Gamma (+) Gamma, e.g. 'Z_3 (+) Z_3'."""
    if not factors:
        return "trivial"
    g = abelian_label(factors)
    return f"{g} (+) {g}"


def _num_global_forms(factors: list[int]) -> int:
    """Number of polarizations of D = Gamma (+) Gamma.  Isolated toric Gamma is
    cyclic (single invariant factor), where the count is sigma(d)."""
    if not factors:
        return 1
    if len(factors) == 1:
        return sigma_divisors(factors[0])
    return None  # non-cyclic Gamma cannot arise for isolated toric; be honest
