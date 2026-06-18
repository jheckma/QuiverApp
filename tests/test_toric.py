"""Tests for the toric Calabi-Yau quiver module.

The central invariant: for every shipped toric quiver the field-theory dimension
(Leigh-Strassler / NSVZ incidence counting) equals the geometric closed form
B - 1 from the toric diagram, and the quiver passes every toric consistency
check (closed superpotential loops, the toric two-term condition, anomaly
cancellation).
"""

import json

import pytest

from conformalmanifold import toric as T
from conformalmanifold.toric import (boundary_lattice_points, convex_hull,
                                     dim_from_polygon, normalized_area,
                                     polygon_signature)


# --------------------------------------------------------------------------
# the whole shipped library: LS dim == B-1, and every quiver is consistent
# --------------------------------------------------------------------------
@pytest.mark.parametrize("q", T.default_toric_library(), ids=lambda q: q.label)
def test_library_dim_agrees_and_valid(q):
    ls = q.dim_conf_ls()
    geo = q.dim_conf_geometric()
    assert geo is not None, f"{q.label} has no toric diagram"
    assert ls == geo, f"{q.label}: LS dim {ls} != B-1 {geo}"
    assert q.validate() == [], f"{q.label}: consistency errors {q.validate()}"
    assert ls >= 1


# --------------------------------------------------------------------------
# pinned known dimensions
# --------------------------------------------------------------------------
@pytest.mark.parametrize("label,expected", [
    ("C3", 2),
    ("conifold", 3),
    ("dP0", 2),
    ("dP1", 3),
    ("L(1,5,2)", 3),
    ("Y(2,1)", 3),
    ("Y(3,1)", 3),
    ("Y(3,2)", 3),
    ("Y(7,3)", 3),
])
def test_known_dimensions(label, expected):
    assert T.make_toric(label).dim_conf_ls() == expected


# --------------------------------------------------------------------------
# the Y^{p,q} family is always dim 3, for a wide grid
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p", range(2, 13))
def test_ypq_family_all_dim_three(p):
    for q in range(1, p):
        from math import gcd
        quiver = T.ypq(p, q)
        assert quiver.dim_conf_ls() == 3, f"Y({p},{q}) LS != 3"
        # geometric cross-check only where the L^{a,b,c} polygon is smooth
        if gcd(p + q, p) == 1 and gcd(p - q, p + q) == 1:
            assert quiver.dim_conf_geometric() == 3, f"Y({p},{q}) B-1 != 3"
        assert quiver.validate() == []


# --------------------------------------------------------------------------
# dP1 is literally Y^{2,1}: same quiver size, same toric signature
# --------------------------------------------------------------------------
def test_dp1_is_y21():
    dp1 = T.make_toric("dP1")
    y21 = T.ypq(2, 1)
    assert dp1.signature() == y21.signature()
    assert (dp1.num_nodes, dp1.num_fields, dp1.num_w_terms) == \
           (y21.num_nodes, y21.num_fields, y21.num_w_terms)


# --------------------------------------------------------------------------
# the brane-tiling Euler identity  E = G + W  (faces - edges + vertices = 0)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("q", T.default_toric_library(), ids=lambda q: q.label)
def test_tiling_euler_relation(q):
    # number of gauge nodes == 2 * area of the toric diagram (minimal phase),
    # and fields == gauge + W terms (Euler characteristic of T^2 tiling = 0).
    assert q.num_fields == q.num_nodes + q.num_w_terms
    if q.diagram:
        assert normalized_area(convex_hull(q.diagram)) == q.num_nodes


# --------------------------------------------------------------------------
# orbifold overlap: for the C^3/Z_K triangle, B - 1 == sum gcd - 1
# --------------------------------------------------------------------------
@pytest.mark.parametrize("K,a,b,c", [
    (3, 1, 1, 1), (6, 1, 2, 3), (7, 1, 2, 4), (9, 1, 2, 6), (12, 1, 4, 7),
])
def test_orbifold_polygon_matches_gcd_formula(K, a, b, c):
    from math import gcd
    # build the C^3/Z_K(a,b,c) junior-simplex triangle the same way orbifold.py
    # does, and check B = sum gcd(K, .) so that B-1 = the character formula.
    pts = [(K, 0), (0, K), (0, 0)]
    for t in range(1, K):
        ra, rb, rc = (t * a) % K, (t * b) % K, (t * c) % K
        if ra + rb + rc == K:
            pts.append((ra, rb))

    def egcd(m, n):
        if n == 0:
            return (m, 1, 0)
        g, x, y = egcd(n, m % n)
        return (g, y, x - (m // n) * y)
    g, u, x3 = egcd(K, a)
    yy = (x3 * b) % K
    mapped = []
    for (px, py) in pts:
        num = py - yy * px
        assert num % K == 0
        mapped.append((px, num // K))
    hull = convex_hull(list(set(mapped)))
    B = len(boundary_lattice_points(hull))
    assert B == gcd(K, a) + gcd(K, b) + gcd(K, c)
    assert dim_from_polygon(hull) == gcd(K, a) + gcd(K, b) + gcd(K, c) - 1


# --------------------------------------------------------------------------
# geometry helpers: Pick's theorem (2*area = 2I + B - 2)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("q", T.default_toric_library(), ids=lambda q: q.label)
def test_picks_theorem(q):
    if not q.diagram:
        return
    hull = convex_hull(q.diagram)
    area2, B, I, _ = polygon_signature(hull)
    assert area2 == 2 * I + B - 2


# --------------------------------------------------------------------------
# the explicit L^{1,5,2} reconstruction has the documented size
# --------------------------------------------------------------------------
def test_l152_size():
    q = T.l152()
    assert (q.num_nodes, q.num_fields, q.num_w_terms) == (6, 16, 10)
    assert q.dim_conf_ls() == 3
    assert q.validate() == []


# --------------------------------------------------------------------------
# make_toric round-trips and rejects nonsense
# --------------------------------------------------------------------------
def test_make_toric_dispatch():
    assert T.make_toric("Y(4,3)").label == "Y(4,3)"
    assert T.make_toric("conifold").label == "conifold"
    with pytest.raises(ValueError):
        T.make_toric("not-a-geometry")
    with pytest.raises(ValueError):
        T.ypq(3, 3)          # need q < p
