"""Validation of the superconformal index."""

import pytest

sympy = pytest.importorskip("sympy")

from sqcdkit import index_series, index_pq, index_symbolic, SuperconformalIndexError

CONIFOLD = {
    "name": "conifold",
    "node_labels": ["n0", "n1"],
    "ranks": [2, 2],
    "arrows": [
        {"label": "A1", "source": 0, "target": 1, "r_charge": "1/2"},
        {"label": "A2", "source": 0, "target": 1, "r_charge": "1/2"},
        {"label": "B1", "source": 1, "target": 0, "r_charge": "1/2"},
        {"label": "B2", "source": 1, "target": 0, "r_charge": "1/2"},
    ],
    "superpotential": [
        {"coefficient": "1", "factors": ["A1", "B1", "A2", "B2"]},
        {"coefficient": "-1", "factors": ["A1", "B2", "A2", "B1"]},
    ],
}

SQCD_SU2_NF3 = {
    "ranks": [2],
    "arrows": [
        {"label": "Q", "source": 0, "target": 1, "r_charge": "1/3"},
        {"label": "Qt", "source": 0, "target": 2, "r_charge": "1/3"},
    ],
    "superpotential": [],
}


def test_conifold_mesons_and_baryons():
    # R=1 (u^2): 4 mesons Tr(A_i B_j) + 6 SU(2) baryons = 10.
    series = index_series(CONIFOLD, 4)
    assert series[0] == 1
    assert series[2] == 10


def test_su2_nf3_s_confinement():
    # SU(2) with three flavors s-confines to 15 free chirals (R=2/3): the
    # gauge integral equals C(6,2)=15 at the meson order, flavor fugacity -> 1.
    series = index_series(SQCD_SU2_NF3, 4, flavor_ranks=[3, 3])
    subs = {
        s: 1
        for c in series.values()
        if hasattr(c, "free_symbols")
        for s in c.free_symbols
    }
    at_unit = {
        k: sympy.expand(c.subs(subs)) if hasattr(c, "subs") else c
        for k, c in series.items()
    }
    assert at_unit[0] == 1
    assert at_unit[2] == 15


def test_full_pq_index_collapses_to_unrefined():
    pq = index_pq(CONIFOLD, 4)
    collapsed = {}
    for (i, j), c in pq.items():
        collapsed[i + j] = sympy.expand(collapsed.get(i + j, 0) + c)
    assert collapsed[0] == 1
    assert collapsed[4] == 10  # the p=q slice reproduces the unrefined u^2 = 10


def test_irrational_r_out_of_scope():
    bad = {
        "ranks": [2, 2],
        "arrows": [{"label": "X", "source": 0, "target": 1, "r_charge": "sqrt(2)/2"}],
        "superpotential": [],
    }
    with pytest.raises(SuperconformalIndexError):
        index_series(bad, 2)


def test_index_symbolic_rational_matches_series():
    # The conifold's a-max R = 1/2 is rational: index_symbolic (a tau series)
    # must reproduce index_series term by term (u^k = tau^(k/2)).
    tau = sympy.Symbol("tau", positive=True)
    expr = index_symbolic(CONIFOLD, 4, derive_r="amax")
    ref = sum(
        sympy.Integer(c) * tau ** sympy.Rational(k, 2)
        for k, c in index_series(CONIFOLD, 4, derive_r="amax").items()
    )
    assert sympy.simplify(expr - ref) == 0


def test_index_symbolic_irrational_R():
    # An injected sqrt(13) target R (a dP1-like mixing) gives a genuine index
    # with irrational tau exponents; at tau=1 it counts 1 + 6 + 3 + 1 = 11.
    tau = sympy.Symbol("tau", positive=True)
    a = (sympy.sqrt(13) - 3) / 12
    Rirr = {
        "A1": sympy.Rational(1, 2) - 3 * a,
        "A2": sympy.Rational(1, 2) + a,
        "B1": sympy.Rational(1, 2) + a,
        "B2": sympy.Rational(1, 2) + a,
    }
    expr = index_symbolic(CONIFOLD, 2, r_charges=Rirr)
    assert expr.has(sympy.sqrt(13))
    assert sympy.simplify(expr.subs(tau, 1) - 11) == 0


def test_index_symbolic_rejects_flavor_nodes():
    # SQCD has arrows to flavor nodes (outside ranks): out of scope here.
    with pytest.raises(SuperconformalIndexError):
        index_symbolic(SQCD_SU2_NF3, 2)
