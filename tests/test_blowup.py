"""Tests for the corner blow-up (chamfer) of a toric diagram."""

import pytest

from conformalmanifold import api
from conformalmanifold import resolution as R
from conformalmanifold.toric import (convex_hull, normalized_area,
                                      polygon_signature)


def _sig(points):
    return polygon_signature(convex_hull(points))


def test_blowup_chamfers_a_singular_corner():
    # Z3xZ3 = the C^3/(Z3xZ3) triangle; corner (0,0) is a C^2/Z3 x C singularity.
    tri = [(0, 0), (3, 0), (0, 3)]
    out = R.blowup_corner(tri, (0, 0))
    hull = convex_hull(out)
    # the corner vertex is gone, replaced by the two one-step chamfer points
    assert (0, 0) not in hull
    assert (1, 0) in out and (0, 1) in out
    # still a valid 2-dimensional diagram, now a quadrilateral
    assert len(hull) == 4
    # the chamfer cut off exactly one unimodular triangle: 2*area drops by 1
    assert normalized_area(hull) == normalized_area(convex_hull(tri)) - 1


def test_blowup_recomputes_the_pipeline():
    before = api.summarize_toric_web([(0, 0), (3, 0), (0, 3)])
    after = api.summarize_toric_web([(0, 0), (3, 0), (0, 3)], blowup=[0, 0])
    # the returned diagram is the rewritten (chamfered) one
    pts = {tuple(p) for p in after["diagram"]["input_points"]}
    assert (0, 0) not in pts and (1, 0) in pts and (0, 1) in pts
    # one fewer gauge group (2*area dropped by one unimodular triangle)
    assert after["conformal"]["num_gauge_groups"] == \
        before["conformal"]["num_gauge_groups"] - 1
    # the resolution is recomputed for the new diagram (#triangles == 2*area)
    assert after["resolution"]["num_triangles"] == \
        after["conformal"]["num_gauge_groups"]


def test_smooth_corner_cannot_be_blown_up():
    # C^3 (unit triangle) and the conifold square are smooth at every corner.
    for verts in ([(0, 0), (1, 0), (0, 1)],
                  [(0, 0), (1, 0), (1, 1), (0, 1)]):
        for c in verts:
            with pytest.raises(ValueError):
                R.blowup_corner(verts, c)


def test_non_corner_is_rejected():
    tri = [(0, 0), (3, 0), (0, 3)]
    # (1, 1) is an interior lattice point, not a convex-hull vertex
    with pytest.raises(ValueError):
        R.blowup_corner(tri, (1, 1))
    # a point not in the diagram at all
    with pytest.raises(ValueError):
        R.blowup_corner(tri, (5, 5))


def test_half_singular_corner_uses_only_the_singular_side():
    # a corner with one unit edge and one length-2 edge: chamfer only the long side.
    verts = [(0, 0), (2, 0), (0, 1)]            # corner (0,0): edges len 2 and 1
    out = R.blowup_corner(verts, (0, 0))
    assert (1, 0) in out                        # one step along the length-2 edge
    assert (0, 0) not in convex_hull(out)


def test_blowup_then_blowup_again():
    # blowing up one corner, then another, stays valid and keeps shrinking the area.
    out1 = R.blowup_corner([(0, 0), (3, 0), (0, 3)], (0, 0))
    a1 = normalized_area(convex_hull(out1))
    out2 = R.blowup_corner(out1, (3, 0))        # (3,0) is still a singular corner
    a2 = normalized_area(convex_hull(out2))
    assert a2 == a1 - 1
    assert len(convex_hull(out2)) >= 3
