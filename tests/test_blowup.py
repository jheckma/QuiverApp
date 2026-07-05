"""Tests for the toric blow-up (add an exceptional divisor / star subdivision).

A blow-up of a point adds a new exceptional ray to the toric diagram -- the
diagram GROWS, the gauge-group count (= 2*area) goes UP by one, and the (p,q)-web
gains an external leg.  The minimal blow-up over a boundary edge A--B with
interior apex C adds the point W = A + B - C.  Iterating the del Pezzo ladder
dP0 (P^2) -> dP1 -> dP2 -> dP3 is exactly this operation.
"""

import pytest

from conformalmanifold import api
from conformalmanifold import resolution as R
from conformalmanifold.toric import convex_hull, normalized_area


def test_blowup_adds_an_exceptional_divisor_and_grows_area():
    # dP0 = P^2: blowing up edge {(1,0),(0,1)} (interior apex (0,0)) adds (1,1).
    dp0 = [(1, 0), (0, 1), (-1, -1)]
    cands = {tuple(c["new_point"]) for c in R.surface_blowup_candidates(dp0)}
    assert (1, 1) in cands
    out = R.surface_blowup(dp0, (1, 1))
    hull = convex_hull(out)
    assert (1, 1) in out
    # exactly one unimodular triangle added: 2*area goes up by one.
    assert normalized_area(hull) == normalized_area(convex_hull(dp0)) + 1
    assert len(hull) == 4                       # triangle -> quadrilateral (dP1)


def test_del_pezzo_ladder():
    # dP0 -> dP1 -> dP2 -> dP3 by repeated blow-ups; nodes 3 -> 4 -> 5 -> 6.
    pts = [(1, 0), (0, 1), (-1, -1)]
    adds = [(1, 1), (-1, 0), (0, -1)]
    expected_nodes = [3, 4, 5, 6]
    for step, want_nodes in enumerate(expected_nodes):
        d = api.summarize_toric_web(pts)
        assert d["conformal"]["num_gauge_groups"] == want_nodes
        assert d["inverse_quiver"]["num_gauge"] == want_nodes
        if step < len(adds):
            cands = {tuple(c["new_point"]) for c in R.surface_blowup_candidates(pts)}
            assert adds[step] in cands
            pts = R.surface_blowup(pts, adds[step])
    # dP3 is the maximal hexagon (6 boundary points, one interior).
    assert len(convex_hull(pts)) == 6


def test_blowup_candidates_lie_outside_and_are_minimal():
    # blow-up of the BASE SURFACE: defined only for cones over a toric surface
    # (exactly one interior lattice point).  C^3 and the conifold have none --
    # there is no compact divisor to blow up (that was the source of the bogus
    # C^3 -> ... -> C^3/Z_MxZ_N "blow-up" chain).
    assert R.surface_blowup_candidates([(0, 0), (1, 0), (0, 1)]) == []          # C^3
    assert R.surface_blowup_candidates([(0, 0), (1, 0), (1, 1), (0, 1)]) == []  # conifold
    for verts in ([(1, 0), (0, 1), (-1, -1)],               # dP0
                  [(-1, 0), (1, 0), (0, 1), (0, -1)]):      # F0
        hull = convex_hull(verts)
        a0 = normalized_area(hull)
        cands = R.surface_blowup_candidates(verts)
        assert cands
        for c in cands:
            W = tuple(c["new_point"])
            assert not R.point_in_convex(W, hull)           # strictly outside
            out = R.surface_blowup(verts, W)
            assert normalized_area(convex_hull(out)) == a0 + 1


def test_blowup_recomputes_the_pipeline_and_reads_the_new_quiver():
    before = api.summarize_toric_web([(1, 0), (0, 1), (-1, -1)])    # dP0, 3 nodes
    after = api.summarize_toric_web([(1, 0), (0, 1), (-1, -1)], surface_blowup=[1, 1])
    pts = {tuple(p) for p in after["diagram"]["input_points"]}
    assert (1, 1) in pts                                    # diagram grew
    assert after["conformal"]["num_gauge_groups"] == \
        before["conformal"]["num_gauge_groups"] + 1
    # the resolution + reconstructed quiver are recomputed for the NEW diagram
    assert after["resolution"]["num_triangles"] == \
        after["conformal"]["num_gauge_groups"]
    assert after["inverse_quiver"]["num_gauge"] == \
        after["conformal"]["num_gauge_groups"]
    # the API advertises the new blow-up sites of the grown diagram
    assert after["diagram"]["surface_blowup_candidates"]


def test_illegal_blowup_sites_are_rejected():
    dp0 = [(1, 0), (0, 1), (-1, -1)]
    # an interior point is not a blow-up site
    with pytest.raises(ValueError):
        R.surface_blowup(dp0, (0, 0))
    # a far-away point is not a minimal star subdivision
    with pytest.raises(ValueError):
        R.surface_blowup(dp0, (5, 5))
    # a point on an existing edge would not grow the area by exactly one
    with pytest.raises(ValueError):
        R.surface_blowup(dp0, (1, 0))


def test_surface_blowup_counts_and_minimal_surfaces():
    # dP0 = P^2: three torus-fixed points -> three blow-up sites (all -> dP1);
    # P^2 is minimal: no -1-curve to blow down.
    dp0 = [(1, 0), (0, 1), (-1, -1)]
    assert {tuple(c["new_point"]) for c in R.surface_blowup_candidates(dp0)} == \
        {(1, 1), (-1, 0), (0, -1)}
    assert R.surface_blowdown_candidates(dp0) == []
    for c in R.surface_blowup_candidates(dp0):
        d = api.summarize_toric_web(dp0, surface_blowup=c["new_point"],
                                    include_inverse=False)
        assert d["conformal"]["num_gauge_groups"] == 4          # dP1
    # F0 = P^1 x P^1: four fixed points, minimal (no -1-curves).
    F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]
    assert len(R.surface_blowup_candidates(F0)) == 4
    assert R.surface_blowdown_candidates(F0) == []
    # the -1-curve counts of the del Pezzos: dP1 -> 1, dP2 -> 3, dP3 -> 6.
    assert len(R.surface_blowdown_candidates([(1, 0), (0, 1), (-1, -1), (0, -1)])) == 1
    assert len(R.surface_blowdown_candidates(
        [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)])) == 3
    assert len(R.surface_blowdown_candidates(
        [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)])) == 6


# ---------------------------------------------------------------------------
# blow-down (the inverse: contract a smooth corner)
# ---------------------------------------------------------------------------
def test_blowdown_is_the_inverse_of_blowup():
    dp0 = [(1, 0), (0, 1), (-1, -1)]
    dp1 = R.surface_blowup(dp0, (1, 1))                             # dP0 -> dP1
    # (1,1) is now a smooth corner and the only blow-down site
    bd = {tuple(c["corner"]) for c in R.surface_blowdown_candidates(dp1)}
    assert bd == {(1, 1)}
    back = R.surface_blowdown(dp1, (1, 1))                          # dP1 -> dP0
    # same POLYGON as dP0 (blowdown returns all its lattice points, so compare
    # hulls, not raw point sets)
    assert convex_hull(back) == convex_hull(dp0)
    assert normalized_area(convex_hull(back)) == \
        normalized_area(convex_hull(dp1)) - 1


def test_blowdown_long_edge_corners():
    """A corner between two LONG edges is contractible when its ear against
    the adjacent boundary LATTICE points is unimodular -- regression for the
    adjacent-hull-corner ear test, which reported NO shrink sites on the
    C^3/Z_3xZ_3 triangle (one could grow into the orbifold but never shrink
    back, breaking grow/shrink inverseness)."""
    Z33 = [(0, 0), (3, 0), (0, 3)]
    cd = R.surface_blowdown_candidates(Z33)
    assert {tuple(c["corner"]) for c in cd} == {(0, 0), (3, 0), (0, 3)}
    for c in cd:
        W = c["corner"]
        for A in c["neighbors"]:                # one primitive step along edge
            assert max(abs(A[0] - W[0]), abs(A[1] - W[1])) == 1
        out = R.surface_blowdown(Z33, c["corner"])
        h = convex_hull(out)
        assert normalized_area(h) == 8 and len(h) == 4
        # exact inverse: growing back must offer the removed corner
        grows = {tuple(g["new_point"]) for g in R.surface_blowup_candidates(out)}
        assert tuple(c["corner"]) in grows


def test_blowdown_candidates_and_minimal_diagrams():
    # C^3 and the conifold are not cones over a compact surface (no interior
    # lattice point): nothing to blow down.  In particular conifold -> C^3 is
    # NOT a divisorial blow-down (only a small/flop transition).
    assert R.surface_blowdown_candidates([(0, 0), (1, 0), (0, 1)]) == []
    assert R.surface_blowdown_candidates([(0, 0), (1, 0), (1, 1), (0, 1)]) == []


def test_blowdown_via_api_recomputes_and_rejects_illegal():
    dp3 = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]  # hexagon, 6 nodes
    before = api.summarize_toric_web(dp3)
    corner = api.summarize_toric_web(dp3)["diagram"]["surface_blowdown_candidates"][0]["corner"]
    after = api.summarize_toric_web(dp3, surface_blowdown=corner)
    assert after["conformal"]["num_gauge_groups"] == \
        before["conformal"]["num_gauge_groups"] - 1
    assert tuple(corner) not in {tuple(p) for p in after["diagram"]["input_points"]}
    # an interior / non-smooth point cannot be blown down
    with pytest.raises(ValueError):
        R.surface_blowdown(dp3, (0, 0))


# ---------------------------------------------------------------------------
# the inverse algorithm is decoupled / gated so the geometry stays responsive
# ---------------------------------------------------------------------------
def test_include_inverse_false_skips_reconstruction_but_keeps_geometry():
    pts = [(0, 0), (3, 0), (0, 3)]                          # Z3xZ3: slow/failing inverse
    d = api.summarize_toric_web(pts, include_inverse=False)
    # geometry is present...
    assert d["conformal"]["num_gauge_groups"] == 9
    assert d["diagram"]["surface_blowup_candidates"] is not None
    assert d["web"]["num_external_legs"] > 0
    # ...but the expensive reconstruction is deferred
    assert d["inverse_quiver"].get("deferred") is True
    assert d["inverse_phases"].get("deferred") is True


def test_large_diagram_gates_the_inverse():
    big = [(0, 0), (6, 0), (6, 6), (0, 6)]                  # 2*area = 72 >> gate
    d = api.summarize_toric_web(big)
    assert d["inverse_quiver"]["available"] is False        # skipped, not hung
    assert "reason" in d["inverse_quiver"]
    assert d["conformal"]["num_gauge_groups"] == 72         # geometry still exact
