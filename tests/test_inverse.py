"""Inverse algorithm: toric diagram -> brane tiling -> quiver gauge theory."""

import pytest

from conformalmanifold.inverse import (
    inverse_quiver, inverse_quiver_json, kasteleyn_newton_polygon,
    forward_extract, canonical_adjacency, phase_invariant, solve_homology,
    urban_renewal, enumerate_toric_phases, integrate_masses, DimerGraph,
    square_gauge_faces, face_polygons, dualize_path, dualize_path_json,
    seiberg_path, quiver_seiberg, _strand_polygon, _trace_faces, _FACE_HAND,
    _integer_kernel, _normalize_tiling, inverse_phases_json)
from conformalmanifold import toric as T
from conformalmanifold.chartable import build_character_table
from conformalmanifold.groups import cyclic
from conformalmanifold.quiver import build_quiver


# (label, vertices, expected #gauge, expected #fields)
CASES = [
    ("C3",       [(0, 0), (1, 0), (0, 1)],                              1, 3),
    ("conifold", [(0, 0), (1, 0), (1, 1), (0, 1)],                      2, 4),
    ("dP0",      [(1, 0), (0, 1), (-1, -1)],                            3, 9),
    ("F0",       [(-1, 0), (1, 0), (0, 1), (0, -1)],                    4, 8),
    ("dP1",      [(1, 0), (0, 1), (-1, -1), (0, -1)],                   4, 10),
    ("dP2",      [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)],          5, 11),
    ("dP3",      [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],  6, 12),
    ("Z2xZ2",    [(0, 0), (2, 0), (0, 2)],                              4, 12),
    ("SPP",      [(0, 0), (2, 0), (1, 1), (0, 1)],                      3, 7),
]


@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_inverse_counts_and_consistency(label, verts, ngauge, nfields):
    t = inverse_quiver(verts)
    assert t.num_gauge == ngauge, label
    assert t.num_fields == nfields, label
    # a consistent brane tiling: counts, anomaly, two-term W, Euler on T^2
    c = t.checks
    assert c["gauge_eq_2area"] and c["fields_eq_sum_det"]
    assert c["white_eq_black"] and c["anomaly_free"] and c["toric_superpotential"]
    assert c["euler_V_minus_E_plus_F"] == 0
    assert t.num_white == t.num_black


@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_inverse_superpotential_is_two_term(label, verts, ngauge, nfields):
    t = inverse_quiver(verts)
    # every field appears in exactly one positive and one negative W term
    pos, neg = {}, {}
    for term in t.superpotential:
        d = pos if term["sign"] > 0 else neg
        for f in term["fields"]:
            d[f] = d.get(f, 0) + 1
    for fld in (f["label"] for f in t.fields):
        assert pos.get(fld) == 1 and neg.get(fld) == 1, (label, fld)
    # total field-slots = 2 * num_fields (each field in two terms)
    assert sum(len(term["fields"]) for term in t.superpotential) == 2 * t.num_fields


def test_inverse_adjacency_anomaly_free():
    t = inverse_quiver([(1, 0), (0, 1), (-1, -1)])    # dP0
    A = t.adjacency_int()
    n = len(A)
    for i in range(n):
        assert sum(A[i]) == sum(A[j][i] for j in range(n))   # in == out


def test_inverse_matches_mckay_for_orbifold():
    # C^3/Z3(1,1,1) = dP0: the reconstructed quiver must match the McKay quiver
    # (3 nodes, 9 fields, cyclic 3-arrow structure) up to relabelling of nodes.
    t = inverse_quiver([(1, 0), (0, 1), (-1, -1)])
    g = cyclic(3, (1, 1, 1))
    mckay = build_quiver(g, build_character_table(g))
    assert t.num_gauge == mckay.num_nodes == 3
    assert t.num_fields == mckay.num_arrows == 9
    # both are the cyclic quiver: each node has 3 outgoing + 3 incoming arrows
    A = t.adjacency_int()
    assert sorted(sum(row) for row in A) == [3, 3, 3]


def test_inverse_z2xz2_matches_mckay():
    # C^3/(Z2 x Z2): 4 nodes, 12 bifundamentals -- matches the McKay quiver.
    t = inverse_quiver([(0, 0), (2, 0), (0, 2)])
    assert t.num_gauge == 4 and t.num_fields == 12
    A = t.adjacency_int()
    assert sorted(sum(row) for row in A) == [3, 3, 3, 3]


def test_inverse_json_shape_and_error():
    out = inverse_quiver_json([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert out["available"] and out["num_gauge"] == 2
    assert len(out["fields"]) == 4
    assert set(out["checks"]) >= {"anomaly_free", "toric_superpotential"}
    # too-small / degenerate input -> graceful error payload
    bad = inverse_quiver_json([(0, 0), (1, 1)])
    assert bad["available"] is False and "error" in bad


def test_inverse_determinism():
    # the placement search is seeded -> identical results across runs
    a = inverse_quiver([(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)])   # dP2
    b = inverse_quiver([(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)])
    assert a.adjacency_int() == b.adjacency_int()


@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_kasteleyn_newton_polygon_recovers_diagram(label, verts, ngauge, nfields):
    # *Independent* certificate: the Newton polygon of the Kasteleyn determinant
    # of the reconstructed tiling must equal the input toric diagram up to
    # GL(2,Z) + translation (the toric diagram is only defined up to SL(2,Z)).
    t = inverse_quiver(verts)
    newt = kasteleyn_newton_polygon(t)
    want = T.convex_hull(verts)
    assert T.normalized_area(newt) == T.normalized_area(want), label
    assert T.gl2z_equiv(newt, want), label


# ---------------------------------------------------------------------------
# forward extractor (DimerGraph -> BraneTiling) round-trip gate
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_forward_extract_roundtrip(label, verts, ngauge, nfields):
    """forward_extract(t.to_dimer()) must reproduce the Gulotta tiling: same
    counts, consistent checks, identical quiver up to relabelling, and the same
    Kasteleyn Newton polygon.  This gates every Seiberg-duality feature."""
    t = inverse_quiver(verts)
    r = forward_extract(t.to_dimer(), verts)
    assert r.num_gauge == t.num_gauge == ngauge, label
    assert r.num_fields == t.num_fields == nfields, label
    assert r.num_white == t.num_white and r.num_black == t.num_black, label
    c = r.checks
    assert c["gauge_eq_2area"] and c["white_eq_black"], label
    assert c["anomaly_free"] and c["toric_superpotential"], label
    assert c["euler_V_minus_E_plus_F"] == 0, label
    # same quiver up to node relabelling
    assert canonical_adjacency(r.adjacency_int()) == \
        canonical_adjacency(t.adjacency_int()), label
    # same toric diagram (homology carried through the combinatorial map)
    assert T.gl2z_equiv(kasteleyn_newton_polygon(r),
                        kasteleyn_newton_polygon(t)), label


@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_solve_homology_from_rotation_system(label, verts, ngauge, nfields):
    """`solve_homology` must reconstruct a certifying edge cochain from the
    rotation system ALONE -- the situation after a Seiberg/urban-renewal move,
    where the stored homology is meaningless.  We wipe the homology to (0,0) and
    require the solver to recover a cochain whose Kasteleyn Newton polygon equals
    the toric diagram up to GL(2,Z) + translation.  This (with the saturated
    integer kernel) is what lets a mutated dual tiling re-certify."""
    d = inverse_quiver(verts).to_dimer()
    for ed in d.edges:                       # erase the geometric homology
        ed["h"] = [0, 0]
    want = T.convex_hull(verts)
    sol = solve_homology(d, want)
    assert sol is not None and len(sol) == len(d.edges), label
    for i, ed in enumerate(d.edges):         # install the solved cochain & verify
        ed["h"] = list(sol[i])
    newt = kasteleyn_newton_polygon(forward_extract(d, verts))
    assert T.gl2z_equiv(newt, want), label
    assert T.normalized_area(newt) == T.normalized_area(want), label


def test_solve_homology_rejects_wrong_diagram():
    """A genuine cochain for F0 must NOT certify against an incompatible diagram
    (a non-equivalent Newton polygon), guarding against false positives."""
    d = inverse_quiver([(-1, 0), (1, 0), (0, 1), (0, -1)]).to_dimer()  # F0
    for ed in d.edges:
        ed["h"] = [0, 0]
    assert solve_homology(d, T.convex_hull([(0, 0), (1, 0), (0, 1)])) is None  # C3


def test_urban_renewal_f0_gives_phase_two():
    """Urban renewal on a square gauge face of F0 phase I (4 nodes / 8 fields)
    must produce phase II (4 nodes / 12 fields), anomaly-free, toric, and with a
    homology cochain that still Newton-certifies to the F0 square."""
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    t = inverse_quiver(F0)
    dimer = t.to_dimer()
    faces, _ = _trace_faces(dimer, _FACE_HAND)
    square = next(orb for orb in faces if len(orb) == 4)
    dual = urban_renewal(dimer, square, F0)
    assert dual is not None
    assert dual.num_gauge == 4 and dual.num_fields == 12
    assert dual.checks["anomaly_free"] and dual.checks["toric_superpotential"]
    assert dual.checks["euler_V_minus_E_plus_F"] == 0
    assert T.gl2z_equiv(kasteleyn_newton_polygon(dual), T.convex_hull(F0))
    # it is a genuinely different phase from the seed
    assert phase_invariant(dual) != phase_invariant(t)


@pytest.mark.parametrize("label,verts,nphases,fields", [
    ("C3",       [(0, 0), (1, 0), (0, 1)],                          1, [3]),
    ("conifold", [(0, 0), (1, 0), (1, 1), (0, 1)],                  1, [4]),
    ("dP0",      [(1, 0), (0, 1), (-1, -1)],                        1, [9]),
    ("dP1",      [(1, 0), (0, 1), (-1, -1), (0, -1)],               1, [10]),
    ("F0",       [(-1, 0), (1, 0), (0, 1), (0, -1)],                2, [8, 12]),
    ("dP2",      [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)],      2, [11, 13]),
    ("dP3",      [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],
                                                             4, [12, 14, 14, 18]),
])
def test_enumerate_toric_phases(label, verts, nphases, fields):
    """Phase enumeration finds exactly the literature's distinct toric phases
    (Feng-Hanany-He hep-th/0104259; Feng-Franco-Hanany-He hep-th/0205144:
    dP0/dP1 1, F0/dP2 2, dP3 4 -- the dP2/dP3 extras need mass integration),
    all certifying to the same Newton polygon (a Seiberg-duality invariant).
    `fields` is the field-count MULTISET: dP3's two 14-field phases (Models II
    and III) are genuinely distinct quivers."""
    phases = enumerate_toric_phases(verts)
    assert len(phases) == nphases, label
    assert sorted(p.num_fields for p in phases) == fields, label
    keys = {phase_invariant(p) for p in phases}
    assert len(keys) == nphases, label                     # all distinct
    want = T.convex_hull(verts)
    for p in phases:
        assert T.gl2z_equiv(kasteleyn_newton_polygon(p), want), label
        c = p.checks
        assert c["anomaly_free"] and c["toric_superpotential"], label
        assert c["euler_V_minus_E_plus_F"] == 0 and c["gauge_eq_2area"], label


def _adjacency_from_arrows(n, arrows):
    A = [[0] * n for _ in range(n)]
    for i, j, m in arrows:
        A[i - 1][j - 1] = m
    return A


def _conj_canonical(A):
    """Canonical key up to node relabelling AND charge conjugation."""
    n = len(A)
    At = [[A[j][i] for j in range(n)] for i in range(n)]
    return min(canonical_adjacency(A), canonical_adjacency(At))


# literature quivers: dP2 from the FHH incidence matrices (hep-th/0104259
# section 5.2); dP3 Models I-IV from the FFHH superpotentials (hep-th/0205144
# section 4.1), fields X_ij = arrow i->j.
DP2_LIT = {
    11: _adjacency_from_arrows(5, [(1, 2, 1), (1, 3, 1), (2, 5, 2), (3, 2, 1),
                                   (3, 5, 1), (4, 1, 1), (4, 3, 1), (5, 1, 1),
                                   (5, 4, 2)]),
    13: _adjacency_from_arrows(5, [(1, 2, 1), (1, 4, 1), (1, 5, 1), (2, 3, 2),
                                   (3, 1, 3), (3, 4, 1), (4, 2, 1), (4, 5, 1),
                                   (5, 3, 2)]),
}
DP3_LIT = {
    "I":   _adjacency_from_arrows(6, [(1, 2, 1), (2, 3, 1), (3, 4, 1),
                                      (4, 5, 1), (5, 6, 1), (6, 1, 1),
                                      (1, 3, 1), (3, 5, 1), (5, 1, 1),
                                      (2, 4, 1), (4, 6, 1), (6, 2, 1)]),
    "II":  _adjacency_from_arrows(6, [(1, 2, 1), (2, 6, 1), (6, 1, 1),
                                      (2, 5, 1), (5, 1, 1), (3, 6, 1),
                                      (6, 4, 1), (4, 3, 1), (3, 5, 1),
                                      (5, 4, 1), (1, 3, 2), (4, 1, 1),
                                      (3, 2, 1)]),
    "III": _adjacency_from_arrows(6, [(4, 1, 1), (1, 5, 2), (5, 4, 1),
                                      (4, 3, 1), (3, 5, 2), (5, 2, 1),
                                      (2, 3, 1), (2, 1, 1), (5, 6, 2),
                                      (6, 4, 1), (6, 2, 1)]),
    "IV":  _adjacency_from_arrows(6, [(4, 1, 1), (1, 6, 2), (6, 4, 3),
                                      (4, 3, 1), (3, 6, 2), (4, 2, 1),
                                      (2, 6, 2), (5, 1, 1), (6, 5, 3),
                                      (5, 3, 1), (5, 2, 1)]),
}


def test_dp2_phases_match_literature_quivers():
    """Both dP2 phases must BE the Feng-Hanany-He quivers (up to node
    relabelling + charge conjugation), not merely have the right field count."""
    phases = enumerate_toric_phases([(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)])
    got = {_conj_canonical(p.adjacency_int()) for p in phases}
    want = {_conj_canonical(A) for A in DP2_LIT.values()}
    assert got == want


def test_dp3_phases_match_literature_models():
    """The four dP3 phases must be exactly Models I-IV of hep-th/0205144, and
    the two 14-field phases (II: one double arrow; III: three double arrows)
    must be the two DIFFERENT literature quivers."""
    phases = enumerate_toric_phases(
        [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)])
    got = {_conj_canonical(p.adjacency_int()) for p in phases}
    want = {_conj_canonical(A) for A in DP3_LIT.values()}
    assert got == want
    assert len(want) == 4                        # II and III really differ


def test_phase_invariant_identifies_charge_conjugate():
    """A quiver and its transpose (charge conjugation = reversing every arrow)
    are the same physical theory and the same toric phase in the standard
    counting (FFHH hep-th/0205144 treats arrow reversal as a symmetry)."""
    class _Fake:
        def __init__(self, A):
            self._A = A
            self.num_gauge = len(A)
            self.num_fields = sum(sum(r) for r in A)

        def adjacency_int(self):
            return [list(r) for r in self._A]

    A = DP3_LIT["III"]
    At = [[A[j][i] for j in range(len(A))] for i in range(len(A))]
    assert phase_invariant(_Fake(A)) == phase_invariant(_Fake(At))
    # but genuinely different quivers stay distinct
    assert phase_invariant(_Fake(A)) != phase_invariant(_Fake(DP3_LIT["II"]))


def test_integrate_masses_roundtrip():
    """Un-integrating an edge (splitting it into a 3-edge path through a new
    2-valent black + 2-valent white vertex = an explicit mass pair) and then
    `integrate_masses` must give back the original theory: same quiver, same
    checks, and the seed's Newton polygon after re-solving the homology."""
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    t = inverse_quiver(F0)
    d = t.to_dimer()
    E, nW, nB = len(d.edges), d.nW, d.nB
    # split edge 0 = (w0, b0):  w0 -e0- b_new -eA- w_new -eB- b0
    w0, b0 = d.edges[0]["w"], d.edges[0]["b"]
    eA, eB, w_new, b_new = E, E + 1, nW, nB
    d.edges[0] = {"w": w0, "b": b_new, "h": [0, 0]}       # e0 re-hung on b_new
    d.edges.append({"w": w_new, "b": b_new, "h": [0, 0]})  # eA
    d.edges.append({"w": w_new, "b": b0, "h": [0, 0]})     # eB
    d.rot_w.append([eA, eB])                               # 2-valent white
    d.rot_b.append([0, eA])                                # 2-valent black
    d.rot_b[b0] = [eB if e == 0 else e for e in d.rot_b[b0]]
    d = DimerGraph(nW + 1, nB + 1, d.edges, d.rot_w, d.rot_b)
    r = integrate_masses(d)
    assert r is not None and r is not d                    # something contracted
    assert r.nW == nW and r.nB == nB and len(r.edges) == E
    back = forward_extract(r, F0)
    assert phase_invariant(back) == phase_invariant(t)
    assert back.checks["toric_superpotential"] and back.checks["anomaly_free"]
    assert back.checks["euler_V_minus_E_plus_F"] == 0
    sol = solve_homology(r, T.convex_hull(F0))
    assert sol is not None                                 # still certifies


def test_integrate_masses_noop_returns_input():
    """A dimer with no 2-valent vertex must come back unchanged (same object),
    preserving its homology cochain."""
    d = inverse_quiver([(0, 0), (1, 0), (1, 1), (0, 1)]).to_dimer()
    assert integrate_masses(d) is d


def test_urban_renewal_spp_square_face_self_dual():
    """SPP's square face visits the same black vertex twice (it carries all 4
    face edges), so corner detection must work at the dart-transition level --
    vertex-level detection returns None for every orientation (regression for
    a confirmed completeness bug).  The SPP square move is self-dual: the dual
    is again SPP (3 nodes, 7 fields, same quiver up to relabelling)."""
    SPP = [(0, 0), (2, 0), (1, 1), (0, 1)]
    t = inverse_quiver(SPP)
    dimer = t.to_dimer()
    faces, _ = _trace_faces(dimer, _FACE_HAND)
    squares = [orb for orb in faces if len(orb) == 4]
    assert len(squares) == 2
    for orb in squares:
        dual = urban_renewal(dimer, orb, SPP)
        assert dual is not None
        assert dual.num_gauge == 3 and dual.num_fields == 7
        assert phase_invariant(dual) == phase_invariant(t)
        assert T.gl2z_equiv(kasteleyn_newton_polygon(dual), T.convex_hull(SPP))
    assert len(enumerate_toric_phases(SPP)) == 1


@pytest.mark.parametrize("label,verts,ngauge,nfields", [
    ("Z2xZ3", [(0, 0), (2, 0), (0, 3)],  6, 18),
    ("Z3xZ3", [(0, 0), (3, 0), (0, 3)],  9, 27),
    ("skew",  [(0, 0), (2, 0), (1, 3)],  6, 18),
    ("Z4xZ4", [(0, 0), (4, 0), (0, 4)], 16, 48),
])
def test_orbifold_honeycomb_triangles(label, verts, ngauge, nfields):
    """Triangle diagrams (= abelian orbifolds C^3/Gamma) are built EXACTLY via
    the quotient-honeycomb construction -- including Z3xZ3 / Z4xZ4, where the
    random Gulotta placement search fails.  Checks, the zig-zag strand-polygon
    certificate, and (where affordable) the independent Kasteleyn certificate
    must all hold; every node of the McKay quiver has 3 in + 3 out arrows."""
    t = inverse_quiver(verts)
    assert t.num_gauge == ngauge and t.num_fields == nfields, label
    c = t.checks
    assert c["gauge_eq_2area"] and c["anomaly_free"], label
    assert c["toric_superpotential"] and c["euler_V_minus_E_plus_F"] == 0, label
    assert "honeycomb" in t.note, label
    want = T.convex_hull(verts)
    assert T.gl2z_equiv(_strand_polygon(t.to_dimer()), want), label
    if t.num_white <= 9:
        assert T.gl2z_equiv(kasteleyn_newton_polygon(t), want), label
    A = t.adjacency_int()
    assert sorted(sum(r) for r in A) == [3] * ngauge, label
    # abelian orbifolds are honeycombs: no square face, so no UR move (their
    # toric phase is unique)
    assert square_gauge_faces(t) == [], label
    # the drawn embedding is a LINEAR image of the planar honeycomb: every
    # true universal-cover edge (white -> black - h, the renderer convention)
    # is one of exactly THREE direction vectors -- no crisscrossing.
    dirs = set()
    for e, f in enumerate(t.fields):
        w = t.white_glob[f["white"]]
        b = t.black_glob[f["black"]]
        h = f["homology"]
        dirs.add((round(b[0] - h[0] - w[0], 6), round(b[1] - h[1] - w[1], 6)))
    assert len(dirs) == 3, (label, dirs)


def test_orbifold_honeycomb_field_R_available():
    """The honeycomb seed must carry the zig-zag leg pair per field, so the
    superconformal per-field R-charges still compute (marginal W + NSVZ)."""
    from conformalmanifold import api
    d = api.summarize_toric_web([(0, 0), (3, 0), (0, 3)])
    fr = d["scft"].get("field_R")
    assert fr and fr["W_marginal"] and fr["nsvz_ok"]


def test_dualize_path_interactive_urban_renewal():
    """EVERY square gauge face of the displayed tiling is an available
    Seiberg-duality (urban renewal) move: F0's seed offers all four nodes;
    dualizing any takes 8 -> 12 fields, and dualizing an offered square face
    of the dual goes back to 8 (involutive up to the phase).  dP2 offers four
    nodes; two are self-dual (11 fields) and two reach the 13-field phase."""
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    seed = inverse_quiver_json(F0)
    assert seed["square_faces"] == [0, 1, 2, 3]     # all four nodes N_f = 2N_c
    for f in seed["square_faces"]:
        d = dualize_path_json(F0, [f])
        assert d["available"] and d["num_fields"] == 12
        assert d["dual_path"] == [f] and d["square_faces"]
    # walk back down on a square face of the dual
    d = dualize_path_json(F0, [0])
    back = [dualize_path_json(F0, [0, f2])["num_fields"]
            for f2 in d["square_faces"]]
    assert 8 in back

    dP2 = [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)]
    seed = inverse_quiver_json(dP2)
    got = {f: dualize_path_json(dP2, [f])["num_fields"]
           for f in seed["square_faces"]}
    assert sorted(got.values()) == [11, 11, 13, 13]

    # illegal moves fail loudly
    with pytest.raises(ValueError):
        dualize_path(F0, [99])
    bad = dualize_path_json(F0, [99])
    assert bad["available"] is False and "no gauge node" in bad["error"]


def _seiberg_rule(A, k):
    """Labelled field-theory Seiberg duality on node k: reverse k's flavors,
    add mesons between its in/out neighbours, integrate out massive
    vector-like pairs."""
    n = len(A)
    B = [row[:] for row in A]
    for j in range(n):
        B[k][j], B[j][k] = A[j][k], A[k][j]
    for i in range(n):
        for j in range(n):
            if i != k and j != k and i != j:
                B[i][j] += A[i][k] * A[k][j]
    for i in range(n):
        for j in range(n):
            m = min(B[i][j], B[j][i])
            B[i][j] -= m
            B[j][i] -= m
    return B


@pytest.mark.parametrize("label,verts", [
    ("F0",  [(-1, 0), (1, 0), (0, 1), (0, -1)]),
    ("dP1", [(1, 0), (0, 1), (-1, -1), (0, -1)]),
    ("dP2", [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)]),
    ("dP3", [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]),
])
def test_dualize_preserves_node_identity(label, verts):
    """Dualizing node k must return the LABELLED quiver of field-theory
    Seiberg duality at node k -- spectators keep their labels, node k keeps
    its label with reversed flavors, mesons connect its neighbours (regression:
    the trace-order renumbering used to scramble labels, so dualizing node 2
    of F0 displayed the same labelled quiver as node 0)."""
    seed = inverse_quiver_json(verts)
    A0 = seed["adjacency"]
    for f in seed["square_faces"]:
        t = dualize_path(verts, [f])
        assert t.adjacency_int() == _seiberg_rule(A0, f), (label, f)
    # composition: dualizing the same node twice is the identity (labelled)
    f = seed["square_faces"][0]
    back = dualize_path(verts, [f, f])
    assert back.adjacency_int() == A0, label


def test_general_seiberg_beyond_dimer_regime():
    """Seiberg duality is allowed on ANY node.  On an N_f != 2N_c node the
    ranks become unequal (N_c -> N_f - N_c): a genuine non-toric phase with no
    brane tiling -- the quiver (with ranks) is still tracked and reported with
    `dimer_available: False`."""
    dP0 = [(1, 0), (0, 1), (-1, -1)]
    # dP0 has no square face; dualizing node 0 (N_f = 3N): rank 1 -> 2, and
    # the quiver is the exceptional-collection dual (3, 3, 6 arrows).
    st = seiberg_path(dP0, [0])
    assert st.tiling is None and st.ranks == [2, 1, 1]
    assert st.adjacency == [[0, 0, 3], [3, 0, 0], [0, 6, 0]]
    assert "non-toric" in st.reason
    # rank-weighted anomaly freedom holds
    n = 3
    for k in range(n):
        assert sum(st.adjacency[i][k] * st.ranks[i] for i in range(n)) == \
            sum(st.adjacency[k][j] * st.ranks[j] for j in range(n))
    # involution: the same node twice returns the seed quiver and ranks
    st2 = seiberg_path(dP0, [0, 0])
    assert st2.ranks == [1, 1, 1]
    assert st2.adjacency == inverse_quiver_json(dP0)["adjacency"]

    # F0: a square-node move keeps the dimer; a non-square follow-up leaves
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    assert seiberg_path(F0, [0]).tiling is not None
    st = seiberg_path(F0, [0, 1])                # node 1 has N_f = 4N there
    assert st.tiling is None and st.ranks == [1, 3, 1, 1]
    d = dualize_path_json(F0, [0, 1])
    assert d["available"] and d["dimer_available"] is False
    assert d["checks"]["anomaly_free"] and d["ranks"] == [1, 3, 1, 1]
    assert d["square_faces"] == []
    # a toric-path json still reports a dimer
    d = dualize_path_json(F0, [0])
    assert d["dimer_available"] is True and d["ranks"] == [1, 1, 1, 1]

    # ill-defined duality (N_f < 2N_c would give rank <= 0) is rejected:
    # conifold node has N_f = 2N (fine), C3's single node has adjoints only
    with pytest.raises(ValueError):
        quiver_seiberg([[0]], [1], 0)            # N_f = 0 -> rank -1


def test_superpotential_tracked_through_general_seiberg():
    """The superpotential is tracked through general (non-toric) Seiberg
    dualities by DWZ mutation: mesons M[ab], dual quarks a*/b*, the Delta
    coupling, and F-term integration of every mass term.  Counts: the dP0
    dual has 9 cubic terms; the double dual (involution) recovers the seed's
    6; F0's [0,1] non-toric phase has 16 and [0,1,1] returns to phase II's 8."""
    dP0 = [(1, 0), (0, 1), (-1, -1)]
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    st = seiberg_path(dP0, [0])
    assert st.W is not None and len(st.W) == 9
    assert all(len(w) == 3 for _, w in st.W)     # all cubic
    assert any(l.startswith("M[") for _, w in st.W for l in w)
    st = seiberg_path(dP0, [0, 0])
    assert len(st.W) == 6                        # the seed's term count
    st = seiberg_path(F0, [0, 1])
    assert st.W is not None and len(st.W) == 16
    st = seiberg_path(F0, [0, 1, 1])
    assert st.ranks == [1, 1, 1, 1] and len(st.W) == 8   # F0 phase II count
    # the json payload carries it
    d = dualize_path_json(dP0, [0])
    assert d["dimer_available"] is False
    assert len(d["superpotential_w"]) == 9
    assert all(t["coeff"] in ("1", "-1") for t in d["superpotential_w"])


def test_dualized_tiling_has_harmonic_torus_layout():
    """A Seiberg-dualized tiling is drawn with the HARMONIC flat-torus
    embedding (each vertex at the centroid of its true universal-cover
    neighbour images, using the solved homology) -- not the schematic
    spanning-tree layout, which tangled the picture."""
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    d = dualize_path(F0, [0])
    acc = {}
    for f in d.fields:
        w, b, h = f["white"], f["black"], f["homology"]
        wp, bp = d.white_glob[w], d.black_glob[b]
        acc.setdefault(("w", w), []).append((bp[0] - h[0], bp[1] - h[1]))
        acc.setdefault(("b", b), []).append((wp[0] + h[0], wp[1] + h[1]))
    for (kind, i), nbrs in acc.items():
        p = d.white_glob[i] if kind == "w" else d.black_glob[i]
        cx = sum(q[0] for q in nbrs) / len(nbrs)
        cy = sum(q[1] for q in nbrs) / len(nbrs)
        assert abs(p[0] - cx) < 1e-6 and abs(p[1] - cy) < 1e-6, (kind, i)


def test_integer_kernel_is_saturated():
    """`_integer_kernel` must return a SATURATED basis (the full integer lattice
    of the rational nullspace, not a finite-index sublattice) -- otherwise a
    correct homology cochain could be unreachable.  Regression for the classic
    failure of clearing each rational nullspace vector independently: the row
    [2,1,1] has an integer solution (0,1,-1) that an index-2 sublattice basis
    cannot reach."""
    import itertools
    basis = _integer_kernel([[2, 1, 1]], 3)
    assert len(basis) == 2
    # every small integer solution of 2x+y+z=0 must lie in the basis' Z-span
    def in_span(target):
        for ca, cb in itertools.product(range(-4, 5), repeat=2):
            v = [ca * basis[0][i] + cb * basis[1][i] for i in range(3)]
            if v == list(target):
                return True
        return False
    sols = [x for x in itertools.product(range(-3, 4), repeat=3)
            if 2 * x[0] + x[1] + x[2] == 0]
    assert all(in_span(x) for x in sols)        # index 1 == fully saturated


def test_canonical_adjacency_permutation_stable():
    """A node relabelling must not change the canonical adjacency key."""
    import itertools
    t = inverse_quiver([(-1, 0), (1, 0), (0, 1), (0, -1)])   # F0
    A = t.adjacency_int()
    n = len(A)
    key = canonical_adjacency(A)
    for p in list(itertools.permutations(range(n)))[:24]:
        B = [[A[p[i]][p[j]] for j in range(n)] for i in range(n)]
        assert canonical_adjacency(B) == key


@pytest.mark.parametrize("label,verts,ngauge,nfields", CASES)
def test_face_polygons_cover_close_and_flag_squares(label, verts, ngauge, nfields):
    """The clickable-dimer face data: one boundary polygon per gauge node, in
    the drawing's universal-cover coordinates."""
    t = _normalize_tiling(inverse_quiver(verts), verts)
    fps = face_polygons(t)
    # one polygon per gauge node, in displayed-node order
    assert [f["node"] for f in fps] == list(range(t.num_gauge)), label
    # square flags reproduce the interactive urban-renewal move list
    assert sorted(f["node"] for f in fps if f["square"]) == square_gauge_faces(t)
    # every dart contributes one polygon corner: sizes sum to 2E, even, >= 4
    sizes = sorted(len(f["poly"]) for f in fps)
    assert sum(sizes) == 2 * t.num_fields, label
    assert all(k % 2 == 0 and k >= 4 for k in sizes), label
    # walk closure: the integer translate accumulated around each face orbit
    # (including the wraparound transition) returns to zero -- the geometric
    # face-cocycle condition that makes the drawn polygon closed
    dimer = t.to_dimer()
    faces, _ = _trace_faces(dimer, _FACE_HAND)
    eh = [f["homology"] for f in t.fields]
    for orb in faces:
        tx = ty = 0
        for k in range(len(orb)):
            (e, _s), (pe, ps) = orb[k], orb[k - 1]
            if ps == 1:
                tx += eh[e][0] - eh[pe][0]
                ty += eh[e][1] - eh[pe][1]
        assert (tx, ty) == (0, 0), (label, orb)


def test_face_polygons_on_mutated_phase():
    """F0 phase II (urban-renewal product, schematic layout): the face data is
    still emitted, sizes {4,8,8,4}, squares matching the move list."""
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    ph = inverse_phases_json(F0)
    assert ph["available"] and ph["num_phases"] == 2
    for p in ph["phases"]:
        fs = p["tiling"]["faces"]
        assert [f["node"] for f in fs] == list(range(p["num_gauge"]))
        assert sorted(f["node"] for f in fs if f["square"]) == p["square_faces"]
        assert sum(len(f["poly"]) for f in fs) == 2 * p["num_fields"]
    sizes = sorted(len(f["poly"]) for f in ph["phases"][1]["tiling"]["faces"])
    assert sizes == [4, 4, 8, 8]
