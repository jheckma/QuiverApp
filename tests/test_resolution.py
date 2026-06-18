"""Tests for toric resolutions: triangulations, dual (p,q) webs, and flops."""

import pytest

from conformalmanifold import resolution as R
from conformalmanifold import toric as T
from conformalmanifold.toric import convex_hull, normalized_area, polygon_signature


CASES = {
    "conifold": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "C3": [(0, 0), (1, 0), (0, 1)],
    "dP0": [(1, 0), (0, 1), (-1, -1)],
    "dP3": [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],
    "Z2xZ2": [(0, 0), (2, 0), (0, 2)],
    "Z3xZ3": [(0, 0), (3, 0), (0, 3)],
    "Y41": T.labc_polygon(3, 5, 4),
    "L194": T.labc_polygon(1, 9, 4),
    "big": [(0, 0), (4, 0), (4, 3), (1, 4), (0, 3)],
}


def _tri_area2(pts, t):
    a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


@pytest.mark.parametrize("verts", CASES.values(), ids=CASES.keys())
def test_triangulation_is_unimodular_and_complete(verts):
    hull = convex_hull(verts)
    pts, tris = R.triangulate(hull)
    area2 = normalized_area(hull)
    # number of triangles == 2*area == number of gauge groups
    assert len(tris) == area2
    # every triangle has lattice area 1/2 (normalized area 1) -> unimodular
    assert all(_tri_area2(pts, t) == 1 for t in tris)
    # the triangles tile the whole polygon (areas add up)
    assert sum(_tri_area2(pts, t) for t in tris) == area2
    # every lattice point is used as a triangulation vertex
    used = set()
    for t in tris:
        used.update(t)
    assert used == set(range(len(pts)))
    # every triangle is oriented CCW
    for t in tris:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        assert R._orient(a, b, c) > 0


@pytest.mark.parametrize("verts", CASES.values(), ids=CASES.keys())
def test_edge_map_manifold(verts):
    pts, tris = R.triangulate(convex_hull(verts))
    for edge, owners in R.edge_map(tris).items():
        assert len(owners) in (1, 2)         # boundary or internal, never >2


@pytest.mark.parametrize("verts", CASES.values(), ids=CASES.keys())
def test_dual_web(verts):
    hull = convex_hull(verts)
    pts, tris = R.triangulate(hull)
    web = R.dual_web(pts, tris, hull)
    B = polygon_signature(hull)[1]
    # one junction per triangle
    assert len(web["junctions"]) == len(tris)
    # one external leg per boundary segment (= B), charges conserved
    legs = web["external_legs"]
    assert len(legs) == B
    assert (sum(l["pq"][0] for l in legs), sum(l["pq"][1] for l in legs)) == (0, 0)
    # internal web edges == internal triangulation edges
    n_internal = sum(1 for o in R.edge_map(tris).values() if len(o) == 2)
    assert len(web["internal_edges"]) == n_internal
    # each external leg perpendicular to its boundary segment is encoded as the
    # outward normal; the junction index is a real triangle
    for l in legs:
        assert 0 <= l["junction"] < len(tris)


def test_flop_validity_and_involution():
    # conifold square has exactly one flop; flopping twice returns the start
    hull = convex_hull(CASES["conifold"])
    pts, tris = R.triangulate(hull)
    fe = R.flippable_edges(pts, tris)
    assert len(fe) == 1
    e = fe[0]
    t1 = R.flop(pts, tris, e)
    assert len(t1) == len(tris)
    assert all(_tri_area2(pts, t) == 1 for t in t1)          # still unimodular
    assert {frozenset(t) for t in t1} != {frozenset(t) for t in tris}  # changed
    # the new triangulation has its own single flop; taking it returns the start
    fe1 = R.flippable_edges(pts, t1)
    assert len(fe1) == 1
    t2 = R.flop(pts, t1, fe1[0])
    assert {frozenset(t) for t in t2} == {frozenset(t) for t in tris}  # involutive


def test_flop_rejects_non_flippable():
    hull = convex_hull(CASES["dP0"])
    pts, tris = R.triangulate(hull)
    # dP0's star triangulation has no flops (unique fine triangulation)
    assert R.flippable_edges(pts, tris) == []
    with pytest.raises(ValueError):
        R.flop(pts, tris, [0, 1])


def test_interior_points():
    # dP0 = P^2 has exactly one interior lattice point (reflexive)
    assert len(R.interior_lattice_points(CASES["dP0"])) == 1
    # the conifold square has none
    assert R.interior_lattice_points(CASES["conifold"]) == []
    # Z3xZ3 triangle has interior points = total - boundary
    hull = convex_hull(CASES["Z3xZ3"])
    allp = R.lattice_points(hull)
    inter = R.interior_lattice_points(hull)
    assert len(inter) == len(allp) - polygon_signature(hull)[1]


def test_is_valid_triangulation_guard():
    hull = convex_hull([(0, 0), (2, 0), (0, 2)])
    pts, tris = R.triangulate(hull)
    assert R.is_valid_triangulation(pts, tris, hull)
    # a bogus "triangulation" that covers all indices and has the right count but
    # contains an oversized (non-unimodular) overlapping triangle must be rejected
    bogus = [(0, 5, 2), (0, 3, 1), (1, 4, 2), (3, 5, 4)]
    assert not R.is_valid_triangulation(pts, bogus, hull)
    # the API falls back to the default when handed a bogus triangulation
    from conformalmanifold import api
    r = api.summarize_toric_web([(0, 0), (2, 0), (0, 2)], bogus)
    assert r["resolution"]["num_triangles"] == 4   # default, not the bogus one


def test_circumcenter_equidistant():
    cc = R.circumcenter((0, 0), (2, 0), (0, 2))
    import math
    d = [math.hypot(cc[0] - x, cc[1] - y) for (x, y) in [(0, 0), (2, 0), (0, 2)]]
    assert max(d) - min(d) < 1e-9
