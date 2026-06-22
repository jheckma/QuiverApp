"""Inverse algorithm: toric diagram -> brane tiling -> quiver gauge theory."""

import pytest

from conformalmanifold.inverse import inverse_quiver, inverse_quiver_json
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
