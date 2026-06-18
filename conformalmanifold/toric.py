"""Toric Calabi-Yau quivers and their conformal manifolds.

Companion to the C^3/Gamma orbifold pipeline (`pipeline.run`).  Where that module
covers the *orbifold* singularities C^3/Gamma, this one covers the broader class
of **isolated toric** Calabi-Yau three-fold singularities -- the conifold, the
Y^{p,q} and L^{a,b,c} families, the del Pezzo cones, and (as a consistency
overlap) the abelian orbifolds C^3/Z_K.  Any other toric singularity can be
added by its toric diagram via `from_diagram` / `ToricDiagram` (dimension and
gauge-group count from the polygon alone, no hand-built tiling required).

Each geometry is specified by

  * a **toric diagram**: a convex lattice polygon in Z^2 (the Newton polygon of
    the Kasteleyn determinant of the brane tiling), defined up to GL(2,Z) +
    translation.  This is the clean *label* distinguishing the geometries.
  * a **quiver gauge theory**: gauge nodes (one U(N d_i) per node), oriented
    bifundamental arrows, and a *toric* (two-term) superpotential.

The dimension of the conformal manifold is given as **two formulations of the
same toric count** (the library tests assert they agree):

  (1) **Field theory.**  Leigh-Strassler / NSVZ marginal-coupling counting on
      the quiver, using only the *toric* (two-term) superpotential couplings:

          dim_C M_conf = (n_gauge + n_W) - rank(M)

      where M is the incidence matrix of the linearised beta-function conditions
      (one row per gauge node, one per superpotential term; one column per field;
      the entry is 1 iff the field touches that node / appears in that term).
      This is exactly the count `pipeline`/`conformal` perform for the orbifolds,
      written directly on the quiver.

  (2) **Geometry.**  The number of boundary lattice points of the toric diagram,
      minus one:

          dim_C M_conf = B - 1,        B = # lattice points on the polygon boundary

      B is the number of external legs of the (p,q) 5-brane web dual to the
      toric diagram.  For C^3/Z_K(a,b,c) the triangle has
      B = gcd(K,a) + gcd(K,b) + gcd(K,c), so B - 1 reproduces the orbifold
      character formula `sum_g fix_Q(g) - 1` (see `conformal.py`).

(1) and (2) are two views of the *same* combinatorial quantity, not independent
derivations: both depend only on the toric data, and where the true conformal
manifold has extra (non-toric) marginal directions they miss them together.

CONVENTION / SCOPE.  This is the conformal-manifold dimension counted in the
gauge + *toric*-superpotential coupling space -- the "beta-deformation sector",
the same convention the orbifold pipeline uses.  For a *generic* toric Calabi-Yau
this is the full conformal manifold (Benvenuti-Hanany, hep-th/0502043: a general
Y^{p,q} is 3-dimensional, matching `Y(p,q)` below).  At *symmetry-enhanced*
special points the true conformal manifold is larger because extra, non-toric
marginal operators open up; the two known cases here are flagged with a `note`:

    * conifold = Y^{1,0}: this count gives 3, the true dimension is 5
      (Benvenuti-Hanany) -- enhanced SU(2)xSU(2) flavour symmetry;
    * C^3 / N=4 SYM:      this count gives 2, the true dimension is 3
      (Leigh-Strassler tau, beta, h).

Only dependency is the standard library plus numpy (already required).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd


# ===========================================================================
# Toric-diagram geometry (pure integer lattice; GL(2,Z) + translation invariant)
# ===========================================================================
def convex_hull(points):
    """Andrew's monotone chain; CCW hull vertices, collinear points dropped."""
    pts = sorted(set(map(tuple, points)))
    if len(pts) <= 2:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def boundary_lattice_points(hull):
    """Set of all lattice points on the polygon boundary (vertices + edge points).

    The count B = len(boundary_lattice_points(hull)) is the number of external
    legs of the (p,q) web and the geometric input to `dim = B - 1`.
    """
    pts = set()
    n = len(hull)
    for i in range(n):
        a = hull[i]
        b = hull[(i + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        g = gcd(abs(dx), abs(dy)) or 1
        for t in range(g):
            pts.add((a[0] + dx * t // g, a[1] + dy * t // g))
    return pts


def normalized_area(hull):
    """2 * Euclidean area = lattice (normalized) area; a GL(2,Z) invariant.

    Equals the number of gauge nodes of the brane tiling for a minimal toric
    phase (twice the area of the toric diagram)."""
    n = len(hull)
    s = 0
    for i in range(n):
        x1, y1 = hull[i]
        x2, y2 = hull[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s)


def interior_lattice_points_count(hull):
    """# interior lattice points I, via Pick: area2 = 2I + B - 2."""
    area2 = normalized_area(hull)
    B = len(boundary_lattice_points(hull))
    return (area2 - B + 2) // 2


def polygon_signature(hull):
    """GL(2,Z)+translation invariant fingerprint: (2*area, B, I, sorted edge
    primitive lengths).  Two diagrams with the same signature are almost always
    GL(2,Z)-equivalent -- enough to label/distinguish the geometries here."""
    area2 = normalized_area(hull)
    bpts = boundary_lattice_points(hull)
    B = len(bpts)
    I = (area2 - B + 2) // 2
    edges = []
    n = len(hull)
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        edges.append(gcd(abs(b[0] - a[0]), abs(b[1] - a[1])))
    return (area2, B, I, tuple(sorted(edges)))


def dim_from_polygon(hull):
    """Closed-form conformal-manifold dimension dim_C M_conf = B - 1."""
    return len(boundary_lattice_points(hull)) - 1


# ===========================================================================
# Exact integer rank (fraction-free is overkill; rationals are plenty and exact)
# ===========================================================================
def _rank_exact(rows):
    """Rank of an integer matrix (list of int lists) computed exactly over Q."""
    M = [[Fraction(x) for x in r] for r in rows]
    if not M:
        return 0
    nrows, ncols = len(M), len(M[0])
    rank = 0
    pcol = 0
    for pcol in range(ncols):
        piv = None
        for r in range(rank, nrows):
            if M[r][pcol] != 0:
                piv = r
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        pv = M[rank][pcol]
        M[rank] = [x / pv for x in M[rank]]
        for r in range(nrows):
            if r != rank and M[r][pcol] != 0:
                f = M[r][pcol]
                M[r] = [a - f * b for a, b in zip(M[r], M[rank])]
        rank += 1
        if rank == nrows:
            break
    return rank


# ===========================================================================
# Toric quiver + Leigh-Strassler / NSVZ conformal-manifold counting
# ===========================================================================
@dataclass
class ToricQuiver:
    """A toric quiver gauge theory and the toric diagram it lives on.

    arrows : dict   field_label -> (src_node, tgt_node)
    W      : list   of (sign:+-1, cycle:tuple[field_label,...]); a *toric* (two-
             term) superpotential -- every field appears in exactly one + and one
             - monomial.
    nodes  : list   gauge-node labels.
    diagram: list   CCW lattice-polygon vertices of the toric diagram (optional;
             used for the geometric B-1 cross-check and the GL(2,Z) signature).
    """
    label: str
    family: str
    description: str
    nodes: list
    arrows: dict
    W: list
    diagram: list = field(default_factory=list)
    note: str = ""   # e.g. symmetry-enhancement caveat on the toric count

    # ---- size -------------------------------------------------------------
    @property
    def num_nodes(self) -> int:
        return len(self.nodes)

    @property
    def num_fields(self) -> int:
        return len(self.arrows)

    @property
    def num_w_terms(self) -> int:
        return len(self.W)

    # ---- field-theory dimension (authoritative) ---------------------------
    def incidence_matrix(self):
        """Rows = gauge nodes then W terms; cols = fields; entry 1 iff incident."""
        fields = sorted(self.arrows)
        fidx = {f: i for i, f in enumerate(fields)}
        rows = []
        for n in self.nodes:                       # one row per gauge node
            r = [0] * len(fields)
            for f, (s, t) in self.arrows.items():
                if s == n or t == n:
                    r[fidx[f]] += 1
            rows.append(r)
        for _sign, cyc in self.W:                  # one row per W term
            r = [0] * len(fields)
            for f in cyc:
                r[fidx[f]] += 1
            rows.append(r)
        return rows

    def dim_conf_ls(self) -> int:
        """dim_C M_conf via Leigh-Strassler / NSVZ counting (exact)."""
        rows = self.incidence_matrix()
        rank = _rank_exact(rows)
        return (self.num_nodes + self.num_w_terms) - rank

    # ---- geometric dimension (closed form) --------------------------------
    def boundary_points(self) -> int:
        return len(boundary_lattice_points(convex_hull(self.diagram))) if self.diagram else 0

    def dim_conf_geometric(self):
        """dim_C M_conf = B - 1 from the toric diagram (None if no diagram)."""
        if not self.diagram:
            return None
        return dim_from_polygon(convex_hull(self.diagram))

    def dim_conf(self) -> int:
        """The conformal-manifold dimension (field-theory value)."""
        return self.dim_conf_ls()

    def signature(self):
        return polygon_signature(convex_hull(self.diagram)) if self.diagram else None

    # ---- consistency ------------------------------------------------------
    def validate(self) -> list[str]:
        """Return a list of consistency errors ([] if the quiver is a good toric
        quiver): closed superpotential loops, the toric two-term condition, and
        anomaly cancellation (in-degree == out-degree, all ranks 1)."""
        errs = []
        # every superpotential field must be an actual arrow
        for sign, cyc in self.W:
            for f in cyc:
                if f not in self.arrows:
                    errs.append(f"W term references unknown field {f!r}: {cyc}")
        if errs:
            return errs
        # every W term is a closed oriented loop
        for sign, cyc in self.W:
            ok = all(self.arrows[cyc[i]][1] == self.arrows[cyc[(i + 1) % len(cyc)]][0]
                     for i in range(len(cyc)))
            if not ok:
                errs.append(f"W term not a closed loop: {cyc}")
        # toric two-term: each field appears in exactly one + and one - monomial
        plus, minus = {}, {}
        for sign, cyc in self.W:
            for f in cyc:
                (plus if sign > 0 else minus)[f] = (plus if sign > 0 else minus).get(f, 0) + 1
        for f in self.arrows:
            if plus.get(f, 0) != 1 or minus.get(f, 0) != 1:
                errs.append(f"toric 2-term fails for {f}: +{plus.get(f,0)} -{minus.get(f,0)}")
        # anomaly free: in-degree == out-degree at every node (all ranks 1)
        for n in self.nodes:
            indeg = sum(1 for (s, t) in self.arrows.values() if t == n)
            outdeg = sum(1 for (s, t) in self.arrows.values() if s == n)
            if indeg != outdeg:
                errs.append(f"anomaly at node {n}: in {indeg} != out {outdeg}")
        return errs

    def __str__(self) -> str:
        d = self.dim_conf_ls()
        lines = [
            f"Toric CY3  {self.label}   [{self.family}]",
            f"  {self.description}",
            f"  quiver: {self.num_nodes} nodes, {self.num_fields} fields, "
            f"{self.num_w_terms} superpotential terms",
        ]
        if self.diagram:
            hull = convex_hull(self.diagram)
            area2, B, I, edges = polygon_signature(hull)
            lines.append(f"  toric diagram: vertices {hull}")
            lines.append(f"    signature (2*area, B_boundary, I_interior, edges) "
                         f"= ({area2}, {B}, {I}, {list(edges)})")
            lines.append(f"  dim_C M_conf = B - 1 = {self.dim_conf_geometric()}  "
                         f"(geometry)")
        lines.append(f"  dim_C M_conf = (n_G + n_W) - rank(M) = "
                     f"({self.num_nodes} + {self.num_w_terms}) - "
                     f"{_rank_exact(self.incidence_matrix())} = {d}  (field theory)")
        if self.note:
            lines.append(f"  note: {self.note}")
        return "\n".join(lines)


# ===========================================================================
# Library of named toric Calabi-Yau three-folds
# ===========================================================================
# Labeling scheme (canonical string labels, parsed by `make_toric`):
#   C3                 flat C^3                         (N=4 SYM)        dim 2
#   conifold           C(T^{1,1}) Klebanov-Witten                       dim 3
#   Y(p,q)             Y^{p,q}  (0<q<p)                                  dim 3
#   L(1,5,2)           L^{1,5,2}  (explicit non-Y example)              dim 3
#   dP0                del Pezzo 0 = C^3/Z_3(1,1,1) = C(P^2)            dim 2
#   dP1                del Pezzo 1 = Y^{2,1}                            dim 3
#
# The *explicit-quiver* library (ToricQuiver, below) is restricted to ISOLATED
# toric CY3 (plus smooth C^3), where the field-theory count and the geometric
# B-1 closed form agree exactly -- the tests assert LS == B-1 there.  Two
# NON-ISOLATED cases are deliberately kept out of the explicit-quiver library
# because their N=2 symmetry enhancement breaks that identity:
#   * A_{n-1} = C^2/Z_n x C   (line of A-singularities) -- dim = n+1 != incidence
#   * SPP, xy = z w^2         (line of A_1 singularities)
# (The C^2/Gamma x C orbifolds are handled by the main orbifold pipeline via the
# character formula.)  The diagram-only catalog (ToricDiagram, further below) is
# NOT so restricted: it includes non-isolated geometries too (e.g. some
# C^3/(Z_n x Z_m)), where only the geometric dim = B-1 is reported -- no LS
# cross-check is asserted, so no false field-theory identity is claimed.
# ---------------------------------------------------------------------------

def c3() -> ToricQuiver:
    """C^3: N=4 SYM read as N=1.  One node, three adjoints, cubic W."""
    arrows = {"X": (0, 0), "Y": (0, 0), "Z": (0, 0)}
    W = [(+1, ("X", "Y", "Z")), (-1, ("X", "Z", "Y"))]
    diagram = [(0, 0), (1, 0), (0, 1)]
    return ToricQuiver("C3", "base", "Flat C^3 (N=4 SYM as N=1).",
                       [0], arrows, W, diagram,
                       note="toric/beta-sector count = 2; the full N=4 conformal "
                            "manifold is 3 (Leigh-Strassler tau, beta, h).")


def conifold() -> ToricQuiver:
    """The conifold C(T^{1,1}); Klebanov-Witten SU(N)xSU(N) with quartic W."""
    arrows = {"A1": (0, 1), "A2": (0, 1), "B1": (1, 0), "B2": (1, 0)}
    W = [(+1, ("A1", "B1", "A2", "B2")), (-1, ("A1", "B2", "A2", "B1"))]
    diagram = [(0, 0), (1, 0), (1, 1), (0, 1)]
    return ToricQuiver("conifold", "conifold",
                       "Conifold cone over T^{1,1} (Klebanov-Witten); = Y^{1,0}.",
                       [0, 1], arrows, W, diagram,
                       note="toric/beta-sector count = 3; the full conifold "
                            "conformal manifold is 5 (Benvenuti-Hanany, "
                            "hep-th/0502043) due to enhanced SU(2)xSU(2).")


def _ypq_arrangement(p: int, q: int):
    """Even (Bresenham) spread of q 'V'-units among p slots; rest are 'Z'."""
    out, acc = [], 0
    for _ in range(p):
        acc += q
        if acc >= p:
            acc -= p
            out.append("V")
        else:
            out.append("Z")
    assert out.count("V") == q, (p, q, out)
    return out


def ypq(p: int, q: int) -> ToricQuiver:
    """Y^{p,q} quiver (0 < q < p).  2p nodes, 4p+2q fields, toric superpotential.

    Construction follows Benvenuti-Franco-Hanany-Martelli-Sparks (hep-th/0411264):
    p 'units' on a ring of 2p nodes; q cubic (V-type) units and p-q quartic
    (Z-type) units, adjacent units sharing a U doublet.
    """
    if not (1 <= q <= p - 1):
        raise ValueError("need 1 <= q <= p-1")
    arr = _ypq_arrangement(p, q)
    M = 2 * p
    nodes = list(range(M))
    arrows, W = {}, []

    def add(label, s, t):
        arrows[label] = (s % M, t % M)

    for i in range(p):
        a, b, c, d = 2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3
        add(f"U{i}_1", a, b)
        add(f"U{i}_2", a, b)
        if arr[i] == "V":
            add(f"V{i}_1", b, c)
            add(f"V{i}_2", b, c)
            add(f"YA{i}", c, a)
            add(f"YB{i}", d, b)
            W.append((+1, (f"U{i}_1", f"V{i}_2", f"YA{i}")))
            W.append((-1, (f"U{i}_2", f"V{i}_1", f"YA{i}")))
            W.append((+1, (f"V{i}_1", f"U{(i+1)%p}_2", f"YB{i}")))
            W.append((-1, (f"V{i}_2", f"U{(i+1)%p}_1", f"YB{i}")))
        else:
            add(f"Z{i}", b, c)
            add(f"YQ{i}", d, a)
            W.append((+1, (f"U{i}_1", f"Z{i}", f"U{(i+1)%p}_2", f"YQ{i}")))
            W.append((-1, (f"U{i}_2", f"Z{i}", f"U{(i+1)%p}_1", f"YQ{i}")))

    diagram = _ypq_polygon(p, q)
    return ToricQuiver(f"Y({p},{q})", "Y^{p,q}",
                       f"Cone over the Sasaki-Einstein Y^{{{p},{q}}}.",
                       nodes, arrows, W, diagram)


def _ypq_polygon(p: int, q: int):
    """Toric diagram of Y^{p,q} = L^{p-q,p+q,p} (the L reference quadrilateral).

    The smooth-L construction needs gcd(p,q)=1; for non-coprime (p,q) (where
    Y^{p,q} is a Z_{gcd} quotient) return [] -- the quiver / LS dimension are
    still well defined, only the explicit lattice polygon is omitted."""
    try:
        return labc_polygon(p - q, p + q, p)
    except ValueError:
        return []


def labc_polygon(a: int, b: int, c: int):
    """L^{a,b,c} quadrilateral toric diagram (Franco et al.); needs gcd(b,c)=1.

    w1=(1,0), w2=(a*k, b), w3=(-a*l, c), w4=(0,0)  with c*k + b*l = 1.
    """
    def egcd(m, n):
        if n == 0:
            return (m, 1, 0)
        g, x, y = egcd(n, m % n)
        return (g, y, x - (m // n) * y)
    g, k, l = egcd(c, b)
    if g != 1:
        raise ValueError(f"gcd(b,c)={g} != 1 (need coprime for smooth L^{{{a},{b},{c}}})")
    return convex_hull([(1, 0), (a * k, b), (-a * l, c), (0, 0)])


# Explicit non-Y L^{1,5,2} quiver, reconstructed from the Hanany-Vegh dimer
# (hep-th/0511063 eq. 5.22) and validated; nodes 0..5, 16 fields, 10 W terms.
_L152_ARROWS = {
    0: (1, 0), 1: (2, 1), 2: (3, 2), 3: (0, 3), 4: (0, 4), 5: (1, 0),
    6: (4, 1), 7: (0, 2), 8: (2, 5), 9: (5, 0), 10: (4, 3), 11: (3, 5),
    12: (5, 4), 13: (4, 1), 14: (5, 4), 15: (1, 5),
}
_L152_W = [
    (+1, (0, 3, 2, 1)), (+1, (4, 6, 5)), (+1, (7, 8, 9)), (+1, (10, 11, 12)),
    (+1, (13, 15, 14)), (-1, (0, 4, 13)), (-1, (1, 5, 7)), (-1, (2, 8, 14, 10)),
    (-1, (3, 11, 9)), (-1, (6, 15, 12)),
]


def l152() -> ToricQuiver:
    """The explicit L^{1,5,2} quiver (a non-Y member of the L^{a,b,c} family)."""
    diagram = labc_polygon(1, 5, 2)
    return ToricQuiver("L(1,5,2)", "L^{a,b,c}",
                       "Cone over L^{1,5,2} (non-Y member).",
                       [0, 1, 2, 3, 4, 5], dict(_L152_ARROWS), list(_L152_W),
                       diagram)


def del_pezzo(n: int) -> ToricQuiver:
    """del Pezzo cones realised as toric quivers:
        dP0 = C^3/Z_3(1,1,1) = C(P^2),  dP1 = Y^{2,1}."""
    if n == 0:
        # dP0: 3 nodes, 3 arrows each i->i+1 (9 fields), cubic toric W.
        nodes = [0, 1, 2]
        arrows = {}
        for a in range(3):
            arrows[f"X{a}"] = (0, 1)
            arrows[f"Y{a}"] = (1, 2)
            arrows[f"Z{a}"] = (2, 0)
        # toric (two-term) cubic superpotential, antisymmetric in the flavour
        # index: + X_a Y_b Z_c - X_a Y_c Z_b  over the 3 even/odd flavour triples.
        W = [
            (+1, ("X0", "Y1", "Z2")), (+1, ("X1", "Y2", "Z0")), (+1, ("X2", "Y0", "Z1")),
            (-1, ("X0", "Y2", "Z1")), (-1, ("X1", "Y0", "Z2")), (-1, ("X2", "Y1", "Z0")),
        ]
        # local-P^2 triangle: B=3 boundary points, single interior point -> dim 2.
        diagram = [(1, 0), (0, 1), (-1, -1)]
        return ToricQuiver("dP0", "del Pezzo",
                           "del Pezzo 0 = C^3/Z_3(1,1,1) = C(P^2).",
                           nodes, arrows, W, diagram)
    if n == 1:
        q = ypq(2, 1)
        q.label, q.family = "dP1", "del Pezzo"
        q.description = "del Pezzo 1 = Y^{2,1}."
        return q
    raise ValueError("del_pezzo: only n = 0, 1 are provided as explicit quivers")


# ===========================================================================
# Toric geometries specified by their toric diagram ALONE (no hand-built tiling)
# ===========================================================================
# There are infinitely many toric Calabi-Yau three-folds -- one per convex
# lattice polygon.  Hand-building a consistent brane tiling for each is a project
# in itself, but the two invariants this package reports follow from the toric
# diagram directly:
#
#     dim_C M_conf      = B - 1                (boundary lattice points - 1)
#     # gauge groups    = 2 * area             (normalized area of the polygon)
#
# `ToricDiagram` exposes exactly those, so ANY toric geometry can be catalogued
# by its polygon -- the del Pezzo / Hirzebruch cones, the C^3/(Z_n x Z_m)
# orbifolds, the full L^{a,b,c} family, the 16 reflexive polygons, etc.
# ---------------------------------------------------------------------------
@dataclass
class ToricDiagram:
    """A toric CY3 given by its toric diagram (convex lattice polygon) alone.

    Reports the conformal-manifold dimension (B-1) and the gauge-group count
    (2*area) without an explicit quiver.  Same `dim_C M_conf` convention as
    `ToricQuiver` (toric / beta-deformation sector)."""
    label: str
    family: str
    description: str
    vertices: list
    note: str = ""

    def hull(self):
        return convex_hull(self.vertices)

    @property
    def num_gauge_groups(self) -> int:
        return normalized_area(self.hull())

    def boundary_points(self) -> int:
        return len(boundary_lattice_points(self.hull()))

    def interior_points(self) -> int:
        return interior_lattice_points_count(self.hull())

    def dim_conf(self) -> int:
        """dim_C M_conf = B - 1."""
        return dim_from_polygon(self.hull())

    def signature(self):
        return polygon_signature(self.hull())

    def __str__(self) -> str:
        h = self.hull()
        area2, B, I, edges = polygon_signature(h)
        lines = [
            f"Toric CY3  {self.label}   [{self.family}]  (diagram only)",
            f"  {self.description}",
            f"  toric diagram: vertices {h}",
            f"    signature (2*area, B_boundary, I_interior, edges) "
            f"= ({area2}, {B}, {I}, {list(edges)})",
            f"  # gauge groups = 2*area = {area2}",
            f"  dim_C M_conf = B - 1 = {self.dim_conf()}",
        ]
        if self.note:
            lines.append(f"  note: {self.note}")
        return "\n".join(lines)


def from_diagram(label, vertices, family="custom", description="", note=""):
    """Build a ToricDiagram from any lattice-polygon vertex list."""
    return ToricDiagram(label, family, description or f"Toric CY3 {label}.",
                        list(vertices), note)


# del Pezzo / Hirzebruch cones (reflexive polygons, one interior point) --------
_DEL_PEZZO_DIAGRAMS = {
    0: [(1, 0), (0, 1), (-1, -1)],                                  # P^2
    1: [(1, 0), (0, 1), (-1, -1), (0, -1)],                         # dP1 = F1 = Y^{2,1}
    2: [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)],                # dP2
    3: [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],        # dP3
}


def del_pezzo_diagram(n: int) -> ToricDiagram:
    """Toric diagram of the cone over del Pezzo n (n=0..3): n+3 gauge groups,
    dim_C M_conf = n+2."""
    if n not in _DEL_PEZZO_DIAGRAMS:
        raise ValueError("del_pezzo_diagram: n must be 0, 1, 2 or 3")
    return ToricDiagram(f"dP{n}", "del Pezzo",
                        f"Cone over the del Pezzo surface dP{n}.",
                        _DEL_PEZZO_DIAGRAMS[n])


def f0_diagram() -> ToricDiagram:
    """Cone over P^1 x P^1 (Hirzebruch F0): 4 gauge groups, dim 3."""
    return ToricDiagram("F0", "Hirzebruch", "Cone over P^1 x P^1 (F0).",
                        [(1, 0), (-1, 0), (0, 1), (0, -1)])


def labc_diagram(a: int, b: int, c: int) -> ToricDiagram:
    """Toric diagram of the L^{a,b,c} cone (needs gcd(b,c)=1)."""
    return ToricDiagram(f"L({a},{b},{c})", "L^{a,b,c}",
                        f"Cone over L^{{{a},{b},{c}}}.",
                        labc_polygon(a, b, c))


def orbifold_znm_diagram(n: int, m: int) -> ToricDiagram:
    """C^3/(Z_n x Z_m) abelian orbifold: the n*m-area triangle with vertices
    (0,0),(n,0),(0,m).  # gauge groups = n*m, dim_C M_conf = n+m+gcd(n,m)-1."""
    return ToricDiagram(f"Z({n},{m})", "C^3/(Z_n x Z_m)",
                        f"C^3/(Z_{n} x Z_{m}) abelian orbifold.",
                        [(0, 0), (n, 0), (0, m)])


def default_toric_diagram_library():
    """A broad spread of named toric CY3 catalogued by diagram (no hand tiling).

    Demonstrates the reach beyond Y^{p,q}/L^{a,b,c}: the del Pezzo and Hirzebruch
    cones, the C^3/(Z_n x Z_m) orbifold family, and a slice of general L^{a,b,c}.
    """
    out = [del_pezzo_diagram(n) for n in range(4)]
    out.append(f0_diagram())
    out += [orbifold_znm_diagram(n, m)
            for (n, m) in [(2, 2), (2, 3), (3, 3), (2, 4), (3, 4)]]
    out += [labc_diagram(a, b, c)
            for (a, b, c) in [(1, 5, 2), (1, 7, 3), (2, 3, 1), (1, 9, 4),
                              (2, 5, 3), (1, 11, 5)]]
    return out


# ---------------------------------------------------------------------------
# registry + dispatch
# ---------------------------------------------------------------------------
_PRESETS = {
    "C3": c3,
    "conifold": conifold,
    "dP0": lambda: del_pezzo(0),
    "dP1": lambda: del_pezzo(1),
    "L(1,5,2)": l152,
}

_DIAGRAM_PRESETS = {
    "dP2": lambda: del_pezzo_diagram(2),
    "dP3": lambda: del_pezzo_diagram(3),
    "F0": f0_diagram,
}


def make_toric(label: str):
    """Build a toric geometry from a label string.

    Returns a `ToricQuiver` (explicit quiver) where one is built in, otherwise a
    `ToricDiagram` (toric diagram + dimension only).  Both expose `dim_conf()`.

    Explicit quiver:  C3, conifold, dP0, dP1, L(1,5,2), Y(p,q)  e.g. 'Y(3,1)'.
    Diagram only:     dP2, dP3, F0, Z(n,m)  e.g. 'Z(2,3)',  L(a,b,c) general,
                      or any polygon via `from_diagram`.
    """
    s = label.strip()
    if s in _PRESETS:
        return _PRESETS[s]()
    if s in _DIAGRAM_PRESETS:
        return _DIAGRAM_PRESETS[s]()
    if s.startswith("Y(") and s.endswith(")"):
        p, q = (int(x) for x in s[2:-1].split(","))
        return ypq(p, q)
    if s.startswith("Z(") and s.endswith(")"):
        n, m = (int(x) for x in s[2:-1].split(","))
        return orbifold_znm_diagram(n, m)
    if s.startswith("L(") and s.endswith(")"):
        a, b, c = (int(x) for x in s[2:-1].split(","))
        if (a, b, c) == (1, 5, 2):
            return l152()           # explicit quiver available
        return labc_diagram(a, b, c)  # general L: diagram + dimension only
    raise ValueError(f"unknown toric label: {label!r}")


def default_toric_library():
    """The named toric CY3 with an explicit, validated quiver (LS == B-1)."""
    out = [c3(), conifold(), del_pezzo(0), del_pezzo(1), l152()]
    out += [ypq(p, q) for (p, q) in [(2, 1), (3, 1), (3, 2), (4, 1), (4, 3),
                                     (5, 2), (6, 5), (7, 3)]]
    return out


def list_toric():
    """Labels of the built-in toric library (explicit quivers + diagram-only)."""
    return ([t.label for t in default_toric_library()]
            + [t.label for t in default_toric_diagram_library()])
