"""Magnetic-quiver Coulomb-branch physics: the monopole formula and the
partition -> quiver construction for nilpotent orbit closures.

All Hilbert series are in the 2Delta convention (a scalar contributes t^1), so
coefficients are integers and the closed forms below are the ones reproduced by
``monopole_hilbert_series``.
"""

import pytest

from conformalmanifold import magnetic as M
from conformalmanifold import partitions as P
from conformalmanifold import api
from conformalmanifold.magnetic import _series_geom, _series_mul


ORDER = 12


def _closed_c2zn(N, order):
    # C^2/Z_N: (1 - t^{2N}) / [(1 - t^2)(1 - t^N)^2]
    inv = _series_mul(_series_geom(2, order),
                      _series_mul(_series_geom(N, order), _series_geom(N, order), order),
                      order)
    res = inv[:]
    if 2 * N <= order:
        for k in range(2 * N, order + 1):
            res[k] -= inv[k - 2 * N]
    return res


def _closed_dn(N, order):
    # SU(2) with N flavors: 1/(1 - t^4) + t^{2(N-2)} / [(1 - t^2)(1 - t^{2(N-2)})]
    a = _series_geom(4, order)
    tmp = _series_mul(_series_geom(2, order), _series_geom(2 * (N - 2), order), order)
    b = [0] * (order + 1)
    sh = 2 * (N - 2)
    for k in range(sh, order + 1):
        b[k] += tmp[k - sh]
    return [a[k] + b[k] for k in range(order + 1)]


# --- Fixture 1: U(1) with N flavors -> C^2/Z_N ----------------------------

@pytest.mark.parametrize("N", [1, 2, 3, 4, 5])
def test_u1_with_N_flavors_is_c2_zn(N):
    q = M.Quiver(gauge=[("U", 1)], flavors=[N], edges={})
    coeffs, info = M.monopole_hilbert_series(q, ORDER)
    assert coeffs == _closed_c2zn(N, ORDER)
    assert info["converged"]


def test_c2z2_coefficients_are_odd_numbers():
    # N=2 gives (1+t^2)/(1-t^2)^2 = 1 + 3 t^2 + 5 t^4 + ...
    q = M.Quiver(gauge=[("U", 1)], flavors=[2], edges={})
    coeffs, _ = M.monopole_hilbert_series(q, ORDER)
    assert [coeffs[2 * k] for k in range(ORDER // 2 + 1)] == [2 * k + 1 for k in range(ORDER // 2 + 1)]
    assert all(coeffs[2 * k + 1] == 0 for k in range((ORDER - 1) // 2 + 1))


# --- Fixture 2: SU(2) with N flavors -> D_N -------------------------------

@pytest.mark.parametrize("N", [3, 4, 5, 6])
def test_su2_with_N_flavors_is_dn(N):
    q = M.Quiver(gauge=[("SU", 2)], flavors=[N], edges={})
    coeffs, info = M.monopole_hilbert_series(q, ORDER)
    assert coeffs == _closed_dn(N, ORDER)


def test_su2_classification():
    # SU(2) with 4 flavors is a good theory (min monopole dim > 1/2).
    q = M.Quiver(gauge=[("SU", 2)], flavors=[4], edges={})
    _c, info = M.monopole_hilbert_series(q, ORDER)
    assert info["classification"] == "good"


# --- Fixture 3: minimal orbit of sl(N), [1]-(1)-...-(1)-[1] ----------------

def test_minimal_orbit_sl3_dimensions():
    # coefficient of t^{2k} = dim of the SU(3) irrep [k,k].
    q = M.Quiver(gauge=[("U", 1), ("U", 1)], flavors=[1, 1], edges={(0, 1): 1})
    coeffs, info = M.monopole_hilbert_series(q, 10)
    def dim_kk(k):
        return (k + 1) * (k + 1) * (2 * k + 2) // 2
    assert [coeffs[2 * k] for k in range(6)] == [dim_kk(k) for k in range(6)]
    assert q.coulomb_dim() == 2
    assert info["classification"] == "good"


def test_minimal_orbit_from_partition_matches_hand_quiver():
    q = M.magnetic_quiver_from_partition("su", (2, 1))
    assert [r for _k, r in q.gauge] == [1, 1]
    assert q.flavors == [1, 1]
    assert q.coulomb_dim() == P.orbit_dim("su", (2, 1)) // 2


# --- Fixture 4: T[SU(N)] and the balanced-node SU(N) symmetry --------------

def test_tsun_ranks_and_symmetry():
    # M((N)) = T[SU(N)]: ranks (N-1, N-2, ..., 1), N flavors on the big node,
    # every gauge node balanced -> SU(N) Coulomb-branch symmetry.
    for N in range(2, 6):
        q = M.magnetic_quiver_from_partition("su", (N,))
        assert sorted(r for _k, r in q.gauge) == list(range(1, N))
        assert q.coulomb_dim() == N * (N - 1) // 2
        assert set(M.balanced_nodes(q)) == set(range(q.n))
        sym = M.coulomb_global_symmetry(q)
        assert sym["display"] == f"SU({N})"


def test_tsu3_hilbert_series_regression():
    # Coulomb branch of T[SU(3)] = nilcone of sl(3); freeze the first coefficients.
    q = M.magnetic_quiver_from_partition("su", (3,))
    coeffs, info = M.monopole_hilbert_series(q, 8)
    assert coeffs == [1, 0, 8, 0, 35, 0, 111, 0, 286]
    assert info["classification"] == "good"


# --- partition -> quiver builder edge cases --------------------------------

def test_a_type_ranks_are_integers_for_all_orbits_up_to_6():
    for N in range(2, 7):
        for part in P.orbits("su", N):
            q = M.magnetic_quiver_from_partition("su", part)
            assert all(isinstance(r, int) and r > 0 for _k, r in q.gauge)
            assert q.coulomb_dim() == P.orbit_dim("su", part) // 2
            # every node of a nilpotent-orbit quiver is balanced
            assert set(M.balanced_nodes(q)) == set(range(q.n))


# --- JSON entry point / API layer -----------------------------------------

def test_magnetic_from_partition_json_shape():
    d = M.magnetic_from_partition_json("su", 3, "2,1", order=8)
    assert d["available"] and d["mode"] == "partition"
    assert d["quiver"]["ranks"] == [1, 1]
    assert d["orbit"]["partition"] == [2, 1]
    assert d["orbit"]["quaternionic_dim"] == 2
    assert d["physics"]["hilbert_series"]["coeffs"][0] == 1
    assert "t^2" in d["physics"]["hilbert_series"]["string"]


def test_api_summarize_magnetic_quiver_partition():
    d = api.summarize_magnetic_quiver(mode="partition", kind="su", N=3, partition="3")
    assert d["available"]
    assert d["physics"]["global_symmetry"]["display"] == "SU(3)"


def test_large_quiver_degrades_gracefully():
    # T[SU(5)] (ranks 4,3,2,1) is too large to sum in-browser: the series is omitted
    # but the quiver, symmetry and orbit data must still come back, fast.
    q = M.magnetic_quiver_from_partition("su", (5,))
    coeffs, info = M.monopole_hilbert_series(q, 8)
    assert coeffs is None and info["too_expensive"] is True
    assert info["classification"] == "good"          # still classified (unit-flux sample)
    d = M.magnetic_from_partition_json("su", 5, "5", order=8)
    assert d["available"]
    assert d["quiver"]["ranks"] == [4, 3, 2, 1]
    assert d["physics"]["global_symmetry"]["display"] == "SU(5)"
    assert d["physics"]["hilbert_series"]["available"] is False
    assert "reason" in d["physics"]["hilbert_series"]


def test_minimal_orbit_large_N_still_computes():
    # minimal orbit of sl(6) is all-U(1): cheap even though N is large.
    d = M.magnetic_from_partition_json("su", 6, "2,1,1,1,1", order=8)
    assert d["physics"]["hilbert_series"]["available"] is True


def test_parse_quiver_linear():
    q = M.parse_quiver("1-2-[3]")
    assert [r for _k, r in q.gauge] == [1, 2]
    assert q.flavors == [0, 3]
    assert q.edges == {(0, 1): 1}
    q2 = M.parse_quiver("[1]-1-1-[1]")
    assert [r for _k, r in q2.gauge] == [1, 1]
    assert q2.flavors == [1, 1]


def test_custom_monopole_series_matches_tsun():
    # "1-2-[3]" is T[SU(3)]: same Coulomb HS as M((3)).
    d = M.monopole_series_json("1-2-[3]", order=8)
    assert d["mode"] == "custom" and d["orbit"] is None
    assert d["physics"]["global_symmetry"]["display"] == "SU(3)"
    assert d["physics"]["hilbert_series"]["coeffs"] == [1, 0, 8, 0, 35, 0, 111, 0, 286]


def test_custom_monopole_matches_partition_builder():
    # a hand-typed min-orbit quiver reproduces the partition-mode series.
    custom = M.monopole_series_json("[1]-1-1-[1]", order=8)
    part = M.magnetic_from_partition_json("su", 3, "2,1", order=8)
    assert custom["physics"]["hilbert_series"]["coeffs"] == part["physics"]["hilbert_series"]["coeffs"]


def test_api_custom_mode():
    d = api.summarize_magnetic_quiver(mode="custom", quiver="1-2-3-[4]", order=6)
    assert d["available"] and d["mode"] == "custom"


def test_identify_theory():
    assert "T[SU(3)]" in M.identify_theory(M.parse_quiver("1-2-[3]"))["name"]
    assert "minimal nilpotent orbit" in M.identify_theory(M.parse_quiver("[1]-1-1-[1]"))["name"]
    # single U(1) with k flavors = C^2/Z_k Kleinian
    q = M.Quiver(gauge=[("U", 1)], flavors=[4], edges={})
    assert "C^2/Z_4" in M.identify_theory(q)["name"]
    q1 = M.Quiver(gauge=[("U", 1)], flavors=[1], edges={})
    assert "free" in M.identify_theory(q1)["name"].lower()


def test_identified_in_json_payloads():
    d = M.monopole_series_json("1-2-3-[4]", order=6)
    assert "T[SU(4)]" in d["physics"]["identified"]["name"]
    d2 = M.magnetic_from_partition_json("su", 4, "2,1,1", order=8)
    assert "minimal" in d2["physics"]["identified"]["name"]


def test_parse_quiver_rejects_garbage():
    with pytest.raises(ValueError):
        M.parse_quiver("1-x-3")
    with pytest.raises(ValueError):
        M.parse_quiver("")


# --- orthosymplectic (gated) ----------------------------------------------

@pytest.mark.parametrize("N", [3, 4, 5, 6])
def test_usp2_equals_su2(N):
    # USp(2) == SU(2): the whole USp path (roots 2e_i, vector +-e_i, C_p dressing)
    # must reproduce the SU(2) series exactly.
    qs = M.Quiver(gauge=[("SU", 2)], flavors=[N], edges={})
    qu = M.Quiver(gauge=[("USp", 2)], flavors=[N], edges={})
    assert M.monopole_hilbert_series(qs, ORDER)[0] == M.monopole_hilbert_series(qu, ORDER)[0]


def test_so_odd_single_node_is_sane():
    q = M.Quiver(gauge=[("SO", 3)], flavors=[4], edges={})
    coeffs, info = M.monopole_hilbert_series(q, ORDER)
    assert coeffs[0] == 1 and all(c >= 0 for c in coeffs)
    assert info["classification"] in ("good", "ugly")


def test_coulomb_dim_kind_aware():
    assert M.Quiver(gauge=[("USp", 2)], flavors=[4], edges={}).coulomb_dim() == 1
    assert M.Quiver(gauge=[("SO", 5)], flavors=[1], edges={}).coulomb_dim() == 2
    assert M.Quiver(gauge=[("SU", 2)], flavors=[4], edges={}).coulomb_dim() == 1


def _osp_coeffs(gauge, flavors, order=12):
    c, _info = M.monopole_hilbert_series(M.Quiver(gauge=gauge, flavors=flavors, edges={}), order)
    return [c[2 * k] for k in range(order // 2 + 1)]


def test_usp2_c2z4_fixture():
    # arXiv:2505.03875: USp(2) with 3 fundamental flavors -> C^2/Z4 (A_3 Kleinian).
    assert _osp_coeffs([("USp", 2)], [3]) == [1, 1, 3, 3, 5, 5, 7]


def test_so4_d_type_fixture():
    # arXiv:2505.03875 eq. 3.11/5.11: SO(4) with 3 vector hypers -- exercises the D_2
    # Pfaffian dressing and the m_2<0 lattice.  Locked against the paper.
    assert _osp_coeffs([("SO", 4)], [3]) == [1, 1, 5, 6, 15, 19, 35]


def test_d_type_is_validated_not_experimental():
    # SO(2n) (D-type) is now validated -> it carries the "validated" note, not "experimental"
    _c, info = M.monopole_hilbert_series(M.Quiver(gauge=[("SO", 4)], flavors=[3], edges={}), 6)
    assert any("validated" in w for w in info["warnings"])
    assert not any("EXPERIMENTAL" in w for w in info["warnings"])


def test_O_node_is_experimental():
    q = M.Quiver(gauge=[("O", 5)], flavors=[1], edges={})
    coeffs, info = M.monopole_hilbert_series(q, 6)
    assert coeffs is not None
    assert any("EXPERIMENTAL" in w for w in info["warnings"])


def test_osp_bifundamental_matter_matches_flavors():
    # a half-hyper bifundamental to an ungauged SO(2N) reproduces N USp fundamental
    # flavors -- validates the bifundamental weight sum and the half-hyper factor.
    from conformalmanifold.magnetic import _bifund_weight_sum, _fund_weight_sum
    for m in (1, 2, 3):
        for N in (1, 2, 3):
            fund = N * _fund_weight_sum("USp", (m,))
            bif = _bifund_weight_sum("USp", 2, (m,), "SO", 2 * N, (0,) * N, half=True)
            assert fund == bif
    # U-U: w flavors == full bifundamental to an ungauged U(w)
    assert _bifund_weight_sum("U", 1, (2,), "U", 3, (0, 0, 0), half=False) == 3 * 2


def test_osp_chain_computes_with_warning():
    q = M.parse_quiver("USp2-SO3")
    assert q.half_edges == {(0, 1)}                 # OSp-OSp defaults to half-hyper
    coeffs, info = M.monopole_hilbert_series(q, 8)
    assert coeffs is not None
    assert any("half-hyper" in w for w in info["warnings"])


def test_so_partition_mode_returns_orbit_data_only():
    d = M.magnetic_from_partition_json("so", 5, "3,1,1", order=8)
    assert d["available"] and d["quiver"] is None and d["physics"] is None
    assert d["orbit"]["quaternionic_dim"] == 3
    assert d["warnings"]                        # experimental notice


def test_parse_partition_validation():
    assert M.parse_partition("[3,2,1]") == (3, 2, 1)
    assert M.parse_partition("2 1 1", N=4) == (2, 1, 1)
    with pytest.raises(ValueError):
        M.parse_partition("2,1", N=4)          # wrong sum
    with pytest.raises(ValueError):
        M.parse_partition("")                   # empty
