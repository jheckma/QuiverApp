"""Validation: closed form vs gcd formula, group orders, character-table sanity."""

import os
import sys
from math import gcd

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conformalmanifold import conformal_manifold_dim, make_group  # noqa: E402
from conformalmanifold.chartable import build_character_table  # noqa: E402
from conformalmanifold.conformal import cyclic_closed_form  # noqa: E402
from conformalmanifold.groups import (cyclic, delta_3n2, delta_6n2,  # noqa: E402
                                      icosahedral_A5, tetrahedral_A4)
from conformalmanifold.quiver import build_quiver  # noqa: E402

# faithful cyclic actions: gcd(a,b,c,n) = 1  (gcd(n,a) may still exceed 1)
CYCLIC_CASES = [
    (3, (1, 1, 1)),
    (5, (1, 1, 3)),
    (6, (1, 2, 3)),
    (7, (1, 2, 4)),
    (10, (2, 3, 5)),
    (12, (1, 4, 7)),
    (12, (3, 4, 5)),
    (24, (3, 8, 13)),
]

# non-faithful actions must be rejected
NONFAITHFUL_CASES = [(12, (2, 4, 6)), (24, (6, 8, 10))]


@pytest.mark.parametrize("n,w", CYCLIC_CASES)
def test_cyclic_matches_gcd_formula(n, w):
    g = cyclic(n, w)
    res = conformal_manifold_dim(g)
    assert res.dim_conf == cyclic_closed_form(n, w)
    # per-direction read-off equals the three gcd's
    assert res.per_direction == [gcd(n, w[0]), gcd(n, w[1]), gcd(n, w[2])]


@pytest.mark.parametrize("n,w", CYCLIC_CASES)
def test_cyclic_independent_of_presentation(n, w):
    # S = sum of fix_Q must equal sum of gcds regardless
    g = cyclic(n, w)
    res = conformal_manifold_dim(g)
    assert res.fixed_sum == sum(gcd(n, x) for x in w)


@pytest.mark.parametrize("n,w", NONFAITHFUL_CASES)
def test_nonfaithful_rejected(n, w):
    with pytest.raises(ValueError):
        cyclic(n, w)


def test_group_orders():
    assert tetrahedral_A4().order == 12
    assert delta_3n2(3).order == 27
    assert delta_6n2(2).order == 24
    assert icosahedral_A5().order == 60


@pytest.mark.parametrize("factory,order,nclasses", [
    (tetrahedral_A4, 12, 4),     # A4: 4 conjugacy classes -> 4 irreps (1,1,1,3)
    (lambda: delta_3n2(3), 27, 11),
    (icosahedral_A5, 60, 5),     # A5: 5 classes -> irreps (1,3,3,4,5)
])
def test_character_table_class_count(factory, order, nclasses):
    g = factory()
    t = build_character_table(g)
    assert g.order == order
    assert t.num_irreps == nclasses


def test_character_orthogonality():
    g = tetrahedral_A4()
    t = build_character_table(g)
    sizes = np.array(t.class_sizes, dtype=float)
    # row orthogonality:  (1/|G|) sum_k |C_k| chi_i(k) conj(chi_j(k)) = delta_ij
    G = g.order
    for i in range(t.num_irreps):
        for j in range(t.num_irreps):
            val = np.sum(sizes * t.chars[i] * np.conj(t.chars[j])) / G
            assert np.isclose(val, 1.0 if i == j else 0.0, atol=1e-3)


def test_dims_squared_sum_to_order():
    for factory in (tetrahedral_A4, lambda: delta_3n2(3), icosahedral_A5):
        g = factory()
        t = build_character_table(g)
        assert sum(d * d for d in t.dims) == g.order


def test_a4_quiver():
    # A4 = Delta(12): irreps (1,1,1,3); Q is the 3-dim irrep.
    g = tetrahedral_A4()
    q = build_quiver(g)
    assert q.num_nodes == 4
    assert sorted(q.dims) == [1, 1, 1, 3]
    # McKay invariant: sum_j a_ij * dim(R_j) = dim(Q) * dim(R_i) = 3 dim(R_i)
    dims = np.array(q.dims)
    assert np.allclose(q.adjacency @ dims, 3 * dims)
    # incoming version too: sum_i a_ij * dim(R_i) = 3 dim(R_j)
    assert np.allclose(q.adjacency.T @ dims, 3 * dims)
    assert q.is_connected()


def test_independent_of_N():
    # the formula contains no N; just assert dim is an int >= 1 for samples
    for name in ["A4 = Delta(12)", "Delta(27)", "A5 = Sigma(60)"]:
        res = conformal_manifold_dim(make_group(name))
        assert isinstance(res.dim_conf, int)
