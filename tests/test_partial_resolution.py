"""Partial resolution of the singular cone: blow up / blow down exceptional
divisors (lattice points) at FIXED polygon.

Blowing every divisor down leaves the singular cone (one cell, legs with
multiplicity); blowing them all up is the full crepant resolution (the default,
unimodular phase).  Cells of 2*area > 1 in between are the residual orbifold
singularities.
"""

import pytest

from conformalmanifold import api
from conformalmanifold import resolution as R
from conformalmanifold.toric import convex_hull, normalized_area

Z33 = [(0, 0), (3, 0), (0, 3)]
Z22 = [(0, 0), (2, 0), (0, 2)]


def test_triangulate_active_subset():
    hull = convex_hull(Z33)
    pts, tris = R.triangulate(hull, active=[(0, 0), (3, 0), (0, 3)])
    assert len(tris) == 1                         # the singular cone: one cell
    assert R.is_valid_subdivision(pts, tris, hull, set(convex_hull(Z33)))
    cells = R.residual_cells(pts, tris)
    assert len(cells) == 1 and cells[0]["area2"] == 9
    # corners are forced even if omitted from `active`
    pts2, tris2 = R.triangulate(hull, active=[(1, 1)])
    assert sum(R._orient(pts2[a], pts2[b], pts2[c]) for (a, b, c) in tris2) == 9
    assert len(tris2) == 3                        # corners + the interior point


def test_full_active_equals_default():
    hull = convex_hull(Z22)
    pts, tris = R.triangulate(hull)
    pts2, tris2 = R.triangulate(hull, active=pts)
    assert pts == pts2 and sorted(tris) == sorted(tris2)
    assert R.residual_cells(pts, tris) == []      # fully resolved: no residues


def test_partial_resolution_api_singular_cone():
    d = api.summarize_toric_web(Z33, active=[[0, 0], [3, 0], [0, 3]],
                                include_inverse=False)
    r = d["resolution"]
    assert not r["fully_resolved"] and r["num_active"] == 3
    assert r["num_triangles"] == 1
    sing = r["residual_singularities"]
    assert len(sing) == 1 and sing[0]["area2"] == 9
    assert sing[0]["label"]                        # identified (Z(3,3)) or generic
    # the singular cone's web: 3 external legs, each of multiplicity 3
    legs = r["web"]["external_legs"]
    assert len(legs) == 3
    assert sorted(l["mult"] for l in legs) == [3, 3, 3]
    # global charge conservation with multiplicity
    assert [sum(l["pq"][0] * l["mult"] for l in legs),
            sum(l["pq"][1] * l["mult"] for l in legs)] == [0, 0]
    # the polygon-level data (quiver of the singular theory) is unchanged
    assert d["conformal"]["num_gauge_groups"] == 9


def test_partial_resolution_api_one_divisor():
    # blow up just the interior divisor of C^3/Z3xZ3
    d = api.summarize_toric_web(Z33, active=[[0, 0], [3, 0], [0, 3], [1, 1]],
                                include_inverse=False)
    r = d["resolution"]
    assert r["num_triangles"] == 3 and r["num_active"] == 4
    assert len(r["residual_singularities"]) == 3
    assert all(c["area2"] == 3 for c in r["residual_singularities"])


def test_partial_resolution_full_active_is_canonical():
    # passing ALL lattice points is the same as passing nothing
    d0 = api.summarize_toric_web(Z22, include_inverse=False)
    lat = d0["resolution"]["lattice_points"]
    d1 = api.summarize_toric_web(Z22, active=lat, include_inverse=False)
    assert d1["resolution"]["fully_resolved"]
    assert d1["resolution"]["triangulation"] == d0["resolution"]["triangulation"]
    assert d1["resolution"]["residual_singularities"] == []


def test_partial_resolution_identifies_named_residue():
    # corners-only Z2xZ2: the single cell is the named Z(2,2) orbifold
    d = api.summarize_toric_web(Z22, active=[[0, 0], [2, 0], [0, 2]],
                                include_inverse=False)
    sing = d["resolution"]["residual_singularities"]
    assert len(sing) == 1 and sing[0]["area2"] == 4
    assert "2" in sing[0]["label"]                 # Z(2,2) from the library


def test_partial_resolution_flops_still_work():
    # F0 corners-only: two area-2 cells sharing a flippable diagonal
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    d = api.summarize_toric_web(F0, active=[[-1, 0], [1, 0], [0, 1], [0, -1]],
                                include_inverse=False)
    r = d["resolution"]
    assert r["num_triangles"] == 2
    assert all(c["area2"] == 2 for c in r["residual_singularities"])
    assert len(r["flippable_edges"]) == 1
    d2 = api.summarize_toric_web(F0, triangulation=r["triangulation"],
                                 flop_edge=r["flippable_edges"][0],
                                 active=[[-1, 0], [1, 0], [0, 1], [0, -1]],
                                 include_inverse=False)
    assert d2["resolution"]["num_triangles"] == 2
    assert d2["resolution"]["triangulation"] != r["triangulation"]


def test_surface_blowup_resets_resolution_state():
    # a base-surface blow-up changes the polygon: the new diagram starts fully
    # resolved regardless of any prior partial-resolution state
    dp0 = [(1, 0), (0, 1), (-1, -1)]
    d = api.summarize_toric_web(dp0, surface_blowup=[1, 1],
                                active=[[1, 0], [0, 1], [-1, -1]],
                                include_inverse=False)
    assert d["resolution"]["fully_resolved"]
