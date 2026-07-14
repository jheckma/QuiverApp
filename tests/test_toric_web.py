"""Tests for the toric (p,q)-web builder: geometry helpers + the web API."""

import math

import pytest

from conformalmanifold import api
from conformalmanifold import toric as T


# --------------------------------------------------------------------------
# pq_web: external legs perpendicular to edges, charges sum to zero
# --------------------------------------------------------------------------
@pytest.mark.parametrize("verts", [
    [(0, 0), (1, 0), (1, 1), (0, 1)],                 # conifold
    [(0, 0), (1, 0), (0, 1)],                          # C3
    [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],  # dP3
    [(0, 0), (2, 0), (0, 2)],                          # Z2 x Z2 (multiplicity-2 edges)
])
def test_web_charge_conservation(verts):
    legs = T.pq_web(T.convex_hull(verts))
    sx = sum(l["pq"][0] for l in legs)
    sy = sum(l["pq"][1] for l in legs)
    assert (sx, sy) == (0, 0)                          # 5-brane charge conservation
    # one leg per boundary lattice segment == B
    assert len(legs) == len(T.boundary_lattice_points(T.convex_hull(verts)))


def test_web_legs_perpendicular_to_edges():
    # each leg (p,q) is orthogonal to its edge direction
    hull = T.convex_hull([(0, 0), (1, 0), (1, 1), (0, 1)])
    n = len(hull)
    for leg in T.pq_web(hull):
        a = hull[leg["edge"]]
        b = hull[(leg["edge"] + 1) % n]
        e = (b[0] - a[0], b[1] - a[1])
        p, q = leg["pq"]
        assert e[0] * p + e[1] * q == 0               # leg perpendicular to edge


# --------------------------------------------------------------------------
# gl2z_equiv: invariance under unimodular maps + translation
# --------------------------------------------------------------------------
def test_gl2z_equiv_under_shear_and_translation():
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    # shear M=[[1,1],[0,1]] then translate by (3,-2)
    sheared = [(x + y + 3, y - 2) for (x, y) in sq]
    assert T.gl2z_equiv(sq, sheared)
    # a genuinely different polygon (triangle) is not equivalent
    assert not T.gl2z_equiv(sq, [(0, 0), (1, 0), (0, 1)])


# --------------------------------------------------------------------------
# identify_toric: recognize named geometries up to GL(2,Z)
# --------------------------------------------------------------------------
def test_identify_conifold_sheared():
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    sheared = [(x + 2 * y - 1, y + 4) for (x, y) in sq]
    g = T.identify_toric(sheared)
    assert isinstance(g, T.ToricQuiver) and g.label == "conifold"


def test_identify_unknown_returns_none():
    # a large generic polygon not in the library
    assert T.identify_toric([(0, 0), (3, 0), (3, 2), (1, 3), (0, 2)]) is None


# --------------------------------------------------------------------------
# summarize_toric_web: the API the web endpoint calls
# --------------------------------------------------------------------------
def test_summarize_conifold():
    r = api.summarize_toric_web([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert r["conformal"]["dim_conf"] == 3
    assert r["conformal"]["num_gauge_groups"] == 2          # 2 * area
    assert r["identified"]["label"] == "conifold"
    assert r["identified"]["has_quiver"] is True
    assert r["web"]["charge_sum"] == [0, 0]
    q = r["quiver"]
    assert q["num_nodes"] == 2 and q["num_w_terms"] == 2
    assert q["dim_conf_ls"] == 3 and q["valid"] is True
    # adjacency is a 2x2 matrix of arrow counts
    assert q["adjacency"] == [[0, 2], [2, 0]]


def test_summarize_dp3_diagram_only():
    r = api.summarize_toric_web([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)])
    assert r["conformal"]["dim_conf"] == 5
    assert r["conformal"]["num_gauge_groups"] == 6
    assert r["identified"]["label"] == "dP3"
    assert r["identified"]["has_quiver"] is False
    assert "quiver" not in r


def test_summarize_generic_unmatched():
    r = api.summarize_toric_web([(0, 0), (3, 0), (3, 2), (1, 3), (0, 2)])
    assert r["identified"]["matched"] is False
    assert r["identified"]["has_quiver"] is False
    # invariants still reported from the polygon
    assert r["conformal"]["dim_conf"] == r["diagram"]["boundary_points"] - 1


def test_summarize_rejects_collinear():
    with pytest.raises(ValueError):
        api.summarize_toric_web([(0, 0), (1, 0), (2, 0)])


# --------------------------------------------------------------------------
# dual (p,q) web: tropical/Legendre junctions -> perpendicular, balanced web
# --------------------------------------------------------------------------
from conformalmanifold import resolution as Rz   # noqa: E402


_WEB_DIAGRAMS = [
    [(0, 0), (1, 0), (1, 1), (0, 1)],                       # conifold
    [(0, 0), (2, 0), (0, 2)],                                # Z2 x Z2
    [(0, 0), (3, 0), (0, 1)],                                # asymmetric
    [(0, 0), (2, 0), (2, 1), (0, 1)],                        # 2x1 rectangle
    [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],    # dP3
]


def _toric_edge_of(tris, t1, t2):
    """The shared lattice edge (i,j) of triangles t1,t2."""
    return tuple(sorted(set(tris[t1]) & set(tris[t2])))


@pytest.mark.parametrize("verts", _WEB_DIAGRAMS)
def test_dual_web_internal_edges_perpendicular(verts):
    hull = T.convex_hull(verts)
    pts, tris = Rz.triangulate(hull)
    web = Rz.dual_web(pts, tris, hull)
    assert web["junction_kind"] == "tropical"          # convex lift was found
    J = web["junctions"]
    for e in web["internal_edges"]:
        t1, t2 = e["tris"]
        i, j = _toric_edge_of(tris, t1, t2)
        d = (pts[j][0] - pts[i][0], pts[j][1] - pts[i][1])     # toric edge dir
        w = (J[t2][0] - J[t1][0], J[t2][1] - J[t1][1])          # web edge dir
        assert abs(d[0] * w[0] + d[1] * w[1]) < 1e-9            # perpendicular
        assert math.hypot(*w) > 1e-6                            # non-degenerate
        # the (p,q) label is parallel to the drawn finite leg
        p, q = e["pq"]
        assert abs(p * w[1] - q * w[0]) < 1e-9
        assert math.gcd(abs(p), abs(q)) == 1                    # primitive charge


@pytest.mark.parametrize("verts", _WEB_DIAGRAMS)
def test_dual_web_local_charge_conservation(verts):
    """At every junction the legs (internal + external) balance to zero."""
    hull = T.convex_hull(verts)
    pts, tris = Rz.triangulate(hull)
    web = Rz.dual_web(pts, tris, hull)
    bal = [[0, 0] for _ in tris]
    for e in web["internal_edges"]:                    # pq points t1 -> t2
        t1, t2 = e["tris"]
        bal[t1][0] += e["pq"][0]; bal[t1][1] += e["pq"][1]
        bal[t2][0] -= e["pq"][0]; bal[t2][1] -= e["pq"][1]
    for leg in web["external_legs"]:                   # outward
        t = leg["junction"]
        bal[t][0] += leg["pq"][0]; bal[t][1] += leg["pq"][1]
    assert all(b == [0, 0] for b in bal)


def test_dual_web_flop_changes_web_visibly():
    """The conifold flop rotates a finite internal edge by 90 degrees (it does
    not collapse to a point as the circumcenter web would)."""
    hull = T.convex_hull([(0, 0), (1, 0), (1, 1), (0, 1)])
    pts, tris = Rz.triangulate(hull)
    e0 = Rz.flippable_edges(pts, tris)[0]
    web1 = Rz.dual_web(pts, tris, hull)
    web2 = Rz.dual_web(pts, Rz.flop(pts, tris, e0), hull)

    def internal_dir(web):
        J, e = web["junctions"], web["internal_edges"][0]
        t1, t2 = e["tris"]
        return (J[t2][0] - J[t1][0], J[t2][1] - J[t1][1])

    d1, d2 = internal_dir(web1), internal_dir(web2)
    assert math.hypot(*d1) > 0.5 and math.hypot(*d2) > 0.5    # both visible
    assert abs(d1[0] * d2[0] + d1[1] * d2[1]) < 1e-9          # rotated 90 degrees


def test_adjacency_matrix_consistency():
    # the adjacency total equals the number of bifundamental fields
    c = T.conifold()
    A = c.adjacency_matrix()
    assert sum(sum(row) for row in A) == c.num_fields
