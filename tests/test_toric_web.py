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
# identify_toric: recognise named geometries up to GL(2,Z)
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


def test_adjacency_matrix_consistency():
    # the adjacency total equals the number of bifundamental fields
    c = T.conifold()
    A = c.adjacency_matrix()
    assert sum(sum(row) for row in A) == c.num_fields
