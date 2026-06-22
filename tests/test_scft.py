"""Superconformal-data tests: orbifold anomalies + toric volume minimisation."""

from fractions import Fraction

import pytest

from conformalmanifold.chartable import build_character_table
from conformalmanifold.groups import cyclic, make_group
from conformalmanifold.quiver import build_quiver
from conformalmanifold.scft import (
    central_charges,
    minimize_volume,
    orbifold_scft,
    toric_scft_json,
)
from conformalmanifold import toric as T


# ---------------------------------------------------------------------------
# 't Hooft-anomaly central-charge helper
# ---------------------------------------------------------------------------
def test_central_charges_N4():
    # N=4 SYM, U(N): gaugino R=1 (N^2) + three adjoint chiral fermions R=-1/3
    # (3 N^2) -> Tr R = 0, Tr R^3 = (8/9) N^2 -> a = c = N^2/4.
    P = Fraction(1)                      # |Gamma| = 1 (trivial)
    Qd = Fraction(3)                     # 3 chirals
    tr_R = P + Qd * Fraction(-1, 3)
    tr_R3 = P + Qd * Fraction(-1, 3) ** 3
    a, c = central_charges(tr_R, tr_R3)
    assert tr_R == 0
    assert a == c == Fraction(1, 4)


# ---------------------------------------------------------------------------
# Orbifold SCFT (exact)
# ---------------------------------------------------------------------------
def _scft_for(group):
    quiver = build_quiver(group, build_character_table(group))
    return orbifold_scft(group, quiver)


@pytest.mark.parametrize("n,a,b,c", [(3, 1, 1, 1), (5, 1, 1, 3), (7, 1, 2, 4)])
def test_orbifold_cyclic_a_equals_c_and_value(n, a, b, c):
    g = cyclic(n, (a, b, c))
    s = _scft_for(g)
    # Tr R = 0  ->  a = c = |Gamma|/4 * N^2
    assert s.tr_R == 0
    assert s.a_eq_c
    assert s.a_coeff == Fraction(g.order, 4)
    assert s.c_coeff == Fraction(g.order, 4)
    assert s.a_over_aN4 == g.order          # |Gamma| times the N=4 value
    assert s.R_chiral == Fraction(2, 3)
    assert s.R_superpotential == 2


def test_orbifold_nonabelian_delta27():
    g = make_group("Delta(27)")
    s = _scft_for(g)
    assert s.tr_R == 0 and s.a_eq_c
    assert s.a_coeff == Fraction(g.order, 4)     # 27/4
    assert "×" in s.gauge_group


# ---------------------------------------------------------------------------
# Toric SCFT via MSY volume minimisation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("corners,expect", [
    ([(0, 0), (1, 0), (0, 1)], 1 / 4),                 # C^3 / S^5
    ([(0, 0), (1, 0), (1, 1), (0, 1)], 27 / 64),       # conifold / T^{1,1}
    ([(1, 0), (0, 1), (-1, -1)], 3 / 4),               # dP0 = C^3/Z3 = C(P^2)
])
def test_toric_central_charge(corners, expect):
    out = toric_scft_json(corners)
    assert out["converged"]
    assert out["a_eq_c"]
    assert out["a"]["val"] == pytest.approx(expect, rel=1e-6)
    assert out["c"]["val"] == pytest.approx(expect, rel=1e-6)
    assert out["R_sum"] == pytest.approx(2.0, abs=1e-6)
    for cr in out["corner_R"]:
        assert cr["R"] > 0


def test_toric_orbifold_crosscheck_dP0():
    # dP0 = C^3/Z3(1,1,1): the toric (MSY) route must match the orbifold route.
    orb = _scft_for(cyclic(3, (1, 1, 1)))
    tor = toric_scft_json([(1, 0), (0, 1), (-1, -1)])
    assert float(orb.a_coeff) == pytest.approx(tor["a"]["val"], rel=1e-6)


def test_toric_conifold_R_charges_half():
    out = toric_scft_json([(0, 0), (1, 0), (1, 1), (0, 1)])
    for cr in out["corner_R"]:
        assert cr["R"] == pytest.approx(0.5, abs=1e-5)


def test_toric_ypq_three_dim_and_converges():
    # Y^{2,1} = dP1; volume minimisation must converge and give a finite a=c.
    g = T.ypq(2, 1)
    hull = T.convex_hull(g.diagram)
    out = toric_scft_json(hull)
    assert out["converged"] and out["a_eq_c"]
    assert 0 < out["a"]["val"] < 1
    assert out["R_sum"] == pytest.approx(2.0, abs=1e-5)


def test_minimize_volume_returns_reeb_with_b1_three():
    g, b, t, ok = minimize_volume([(0, 0), (1, 0), (0, 1)])
    assert b[0] == pytest.approx(3.0)
    assert g == pytest.approx(3.0, rel=1e-6)        # C^3 -> g = 3
