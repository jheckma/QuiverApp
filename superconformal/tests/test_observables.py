"""Validation of the SCFT observables against textbook results.

Self-contained: SQCD and the conifold are built inline; dP0 (SU(3)^3) and
dP1 (SU(2)^4) are loaded from tests/fixtures/. No external project needed.
"""

import json
from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")

from sqcdkit import (
    scft_observables,
    superconformal_central_charges,
    one_loop_beta_coefficients,
    flavor_thooft_anomalies,
    abelian_flavor_anomalies,
)

FIX = Path(__file__).parent / "fixtures"


def _fixture(name):
    return json.loads((FIX / name).read_text())


def _sqcd(Nc, Nf):
    # SU(Nc) gauge (node 0) + two SU(Nf) global flavor nodes (1, 2).
    return {
        "name": f"sqcd_{Nc}_{Nf}",
        "node_labels": ["g0"],
        "ranks": [Nc],
        "arrows": [
            {"label": "Q", "source": 0, "target": 1, "r_charge": "1/2"},
            {"label": "Qb", "source": 2, "target": 0, "r_charge": "1/2"},
        ],
        "superpotential": [],
    }


# ----------------------------------------------------------------------
# SQCD (flavored): the canonical anomaly-fixed results.
# ----------------------------------------------------------------------


def test_sqcd_superconformal_R():
    for Nc, Nf in [(2, 3), (3, 6), (3, 4), (4, 6)]:
        R = superconformal_central_charges(_sqcd(Nc, Nf), flavor_ranks=[Nf, Nf]).r_charges
        expected = 1 - sympy.Rational(Nc, Nf)
        assert sympy.simplify(R["Q"] - expected) == 0
        assert sympy.simplify(R["Qb"] - expected) == 0


def test_sqcd_one_loop_beta():
    for Nc, Nf in [(2, 3), (3, 6), (3, 4)]:
        b0 = one_loop_beta_coefficients(_sqcd(Nc, Nf), flavor_ranks=[Nf, Nf])
        assert b0[0] == 3 * Nc - Nf


def test_sqcd_flavor_anomalies():
    for Nc, Nf in [(2, 3), (3, 6)]:
        R = superconformal_central_charges(_sqcd(Nc, Nf), flavor_ranks=[Nf, Nf]).r_charges
        anom = flavor_thooft_anomalies(_sqcd(Nc, Nf), R, flavor_ranks=[Nf, Nf])
        for d in anom.values():
            assert abs(d["SU3"]) == Nc  # cubic flavor anomaly, matched by Seiberg
            assert sympy.simplify(d["SU2_R"] - sympy.Rational(-(Nc ** 2), 2 * Nf)) == 0


def test_sqcd_baryonic_anomalies():
    R = superconformal_central_charges(_sqcd(2, 3), flavor_ranks=[3, 3]).r_charges
    ab = abelian_flavor_anomalies(_sqcd(2, 3), R, flavor_ranks=[3, 3])
    assert ab["n_u1"] == 1                 # U(1)_B only
    assert ab["grav2"][0] == 0             # Tr B = 0
    assert ab["cubic"][(0, 0, 0)] == 0     # B^3 = 0 (charge conjugation)


def test_hofman_maldacena_gate():
    # SU(3) Nf=4 (= Nc+1, s-confining, below the window): a/c < 1/2, flagged.
    assert not scft_observables(_sqcd(3, 4), flavor_ranks=[4, 4])["hofman_maldacena"]["ok"]
    # SU(2) Nf=3 saturates a/c = 1/2 (its IR is free).
    obs = scft_observables(_sqcd(2, 3), flavor_ranks=[3, 3])
    assert obs["hofman_maldacena"]["ok"]
    assert abs(obs["hofman_maldacena"]["a_over_c"] - 0.5) < 1e-9


# ----------------------------------------------------------------------
# Kutasov: adjoint matter with a higher-order single-trace superpotential.
# ----------------------------------------------------------------------


def _kutasov(k, Nf):
    # SU(2) + adjoint Phi + Nf flavors, W = Tr Phi^(k+1).
    return {
        "name": f"kutasov_k{k}",
        "node_labels": ["g0"],
        "ranks": [2],
        "arrows": [
            {"label": "Phi", "source": 0, "target": 0, "r_charge": "1/2"},
            {"label": "Q", "source": 0, "target": 1, "r_charge": "1/2"},
            {"label": "Qb", "source": 2, "target": 0, "r_charge": "1/2"},
        ],
        "superpotential": [{"coefficient": "1", "factors": ["Phi"] * (k + 1)}],
    }


def test_kutasov_adjoint_r_charge():
    # W = Tr Phi^(k+1) pins the adjoint at R(Phi) = 2/(k+1) (Kutasov 1995).
    for k in (2, 3, 4, 5):
        R = superconformal_central_charges(_kutasov(k, 4), flavor_ranks=[4, 4]).r_charges
        assert R["Phi"] == sympy.Rational(2, k + 1)


# ----------------------------------------------------------------------
# Toric quivers: exact rational and irrational central charges.
# ----------------------------------------------------------------------


def test_dp0_central_charges():
    # dP0 = C^3/Z_3 at SU(3)^3: a = 99/16, c = 51/8 (Tr R = -3, Tr R^3 = 21).
    obs = scft_observables(_fixture("dp0.json"))
    assert obs["tr_R"] == "-3" and obs["tr_R3"] == "21"
    assert all(v == "0" for v in obs["one_loop_b0"].values())
    assert obs["hofman_maldacena"]["ok"]


def test_dp1_exact_irrational_central_charge():
    # dP1 = Y^{2,1} at SU(2)^4: a = -739/4 + 52 sqrt(13), c = -369/2 + 52 sqrt(13).
    res = superconformal_central_charges(_fixture("dp1.json"))
    a_expected = sympy.Rational(-739, 4) + 52 * sympy.sqrt(13)
    c_expected = sympy.Rational(-369, 2) + 52 * sympy.sqrt(13)
    assert sympy.simplify(res.a - a_expected) == 0
    assert sympy.simplify(res.c - c_expected) == 0
    assert res.exact


def test_spp_irrational_central_charge():
    # SPP (suspended pinch point), a non-orbifold toric quiver with an adjoint:
    # a-max gives irrational R (sqrt(97)); a, c exact, HM in window, no decoupling.
    res = superconformal_central_charges(_fixture("spp.json"))
    assert sympy.simplify(res.r_charges["X22"] - (sympy.Rational(5, 2) - sympy.sqrt(97) / 6)) == 0
    assert sympy.simplify(res.a - (sympy.Rational(-189, 64) + 97 * sympy.sqrt(97) / 192)) == 0
    assert sympy.simplify(res.c - (sympy.Rational(-171, 64) + 95 * sympy.sqrt(97) / 192)) == 0
    obs = scft_observables(_fixture("spp.json"))
    assert obs["hofman_maldacena"]["ok"]
    assert obs["mesonic_below_bound"] == ()


def test_c3_a4_nonabelian_orbifold():
    # C^3/A4 (tetrahedral) McKay quiver at N=2 (ranks [6,2,2,2]): every field at
    # R=2/3, a=3N^2-19/24=269/24, c=3N^2-7/12=137/12, matching the published
    # values of Fang-Feng-Xie (arXiv:2310.15792).
    res = superconformal_central_charges(_fixture("c3_a4.json"))
    assert all(sympy.simplify(v - sympy.Rational(2, 3)) == 0 for v in res.r_charges.values())
    assert sympy.simplify(res.a - sympy.Rational(269, 24)) == 0
    assert sympy.simplify(res.c - sympy.Rational(137, 12)) == 0
    assert scft_observables(_fixture("c3_a4.json"))["hofman_maldacena"]["ok"]
