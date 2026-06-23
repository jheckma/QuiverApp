"""Inverse algorithm: toric diagram -> brane tiling -> quiver gauge theory."""

import pytest

from conformalmanifold.inverse import (
    inverse_quiver, inverse_quiver_json, kasteleyn_newton_polygon,
    forward_extract, canonical_adjacency, phase_invariant, solve_homology,
    urban_renewal, enumerate_toric_phases, _trace_faces, _FACE_HAND,
    _integer_kernel)
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


@pytest.mark.parametrize("label,verts,nphases,fieldset", [
    ("C3",       [(0, 0), (1, 0), (0, 1)],                 1, {3}),
    ("conifold", [(0, 0), (1, 0), (1, 1), (0, 1)],         1, {4}),
    ("F0",       [(-1, 0), (1, 0), (0, 1), (0, -1)],       2, {8, 12}),
])
def test_enumerate_toric_phases(label, verts, nphases, fieldset):
    """Phase enumeration finds exactly the expected distinct toric phases, all
    certifying to the same Newton polygon (a Seiberg-duality invariant)."""
    phases = enumerate_toric_phases(verts)
    assert len(phases) == nphases, label
    assert {p.num_fields for p in phases} == fieldset, label
    keys = {phase_invariant(p) for p in phases}
    assert len(keys) == nphases, label                     # all distinct
    want = T.convex_hull(verts)
    for p in phases:
        assert T.gl2z_equiv(kasteleyn_newton_polygon(p), want), label


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
