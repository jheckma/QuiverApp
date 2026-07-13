"""Finite-chamber BPS quiver tests."""

import pytest

from conformalmanifold import api
from conformalmanifold.bps import bps_quiver_json, parse_arrows, parse_matrix


A2_EXCHANGE = "0 1; -1 0"


def test_a2_maximal_green_sequence_stable_charges():
    out = bps_quiver_json(A2_EXCHANGE, "1,2,1")

    assert out["chamber"]["green_sequence"] is True
    assert out["chamber"]["maximal_green"] is True
    assert out["chamber"]["stable_count"] == 3
    assert out["chamber"]["warnings"] == []
    assert [obj["dimension_vector"] for obj in out["stable_objects"]] == [
        [1, 0],
        [1, 1],
        [0, 1],
    ]


def test_partial_green_sequence_warns_when_not_maximal():
    out = bps_quiver_json(A2_EXCHANGE, "1,2")

    assert out["chamber"]["green_sequence"] is True
    assert out["chamber"]["maximal_green"] is False
    assert out["chamber"]["stable_count"] == 2
    assert len(out["chamber"]["warnings"]) == 1
    assert "not maximal green" in out["chamber"]["warnings"][0]
    assert "partial finite chamber" in out["chamber"]["warnings"][0]


def test_non_green_sequence_warns_and_reports_green_stable_objects_only():
    out = bps_quiver_json(A2_EXCHANGE, "1,1")

    assert out["chamber"]["green_sequence"] is False
    assert out["chamber"]["maximal_green"] is False
    assert [m["status"] for m in out["mutations"]] == ["green", "red"]
    assert [obj["dimension_vector"] for obj in out["stable_objects"]] == [[1, 0]]
    assert len(out["chamber"]["warnings"]) == 1
    assert "not green throughout" in out["chamber"]["warnings"][0]


def test_arrow_parser_accepts_multiplicity_suffix():
    parsed = parse_arrows("1->2 x2")

    assert parsed.labels_are_zero_based is False
    assert parsed.matrix == [[0, 2], [-2, 0]]


def test_invalid_exchange_matrix_raises_value_error():
    with pytest.raises(ValueError, match="skew-symmetric"):
        parse_matrix("0 1; 0 0")

def test_zero_based_arrow_input_round_trips_node_labels():
    out = bps_quiver_json("0->1", "0 1 0", kind="arrows")

    assert out["quiver"]["label_base"] == 0
    assert out["quiver"]["node_labels"] == ["gamma_0", "gamma_1"]
    assert out["chamber"]["mutation_sequence"] == [0, 1, 0]
    assert [obj["node"] for obj in out["stable_objects"]] == [0, 1, 0]
    assert [obj["charge"] for obj in out["stable_objects"]] == [
        "gamma_0",
        "gamma_0 + gamma_1",
        "gamma_1",
    ]


def test_arrow_parser_rejects_two_cycles():
    with pytest.raises(ValueError, match="2-cycle"):
        parse_arrows("1->2; 2->1")


def test_masses_must_be_positive():
    with pytest.raises(ValueError, match="positive"):
        bps_quiver_json(A2_EXCHANGE, "1,2,1", phases_text="0.6 0.2", masses_text="1 0")



def test_bps_toric_geometry_uses_existing_inverse_quiver():
    out = api.summarize_bps_toric_quiver([(1, 0), (0, 1), (-1, -1)], "0 1 2")

    assert out["source"]["kind"] == "toric_geometry"
    assert out["source"]["quiver_source"] == "inverse brane tiling"
    assert [0, 0] in out["source"]["lattice_points"]
    assert out["source"]["sw_curve"].startswith("H(x,y) = ")
    assert "c[0,0]" in out["source"]["sw_curve"]
    assert out["quiver"]["adjacency"] == [[0, 3, 0], [0, 0, 3], [3, 0, 0]]
    assert out["quiver"]["exchange_matrix"] == [[0, 3, -3], [-3, 0, 3], [3, -3, 0]]
    assert out["chamber"]["green_sequence"] is True
    assert [obj["dimension_vector"] for obj in out["stable_objects"]] == [
        [1, 0, 0],
        [3, 1, 0],
        [0, 0, 1],
    ]


def test_bps_toric_geometry_reports_reciprocal_netting():
    out = api.summarize_bps_toric_quiver([(0, 0), (1, 0), (1, 1), (0, 1)], "0 1")

    assert out["quiver"]["adjacency"] == [[0, 2], [2, 0]]
    assert out["quiver"]["exchange_matrix"] == [[0, 0], [0, 0]]
    assert any("reciprocal arrows" in w for w in out["chamber"]["warnings"])


# --- canonical orbifold chamber + invariants + spectral network (5d SCFTs) ---

_FIVED = {
    "dP0":  ([(1, 0), (0, 1), (-1, -1)], 3, 0),
    "F0":   ([(-1, 0), (1, 0), (0, 1), (0, -1)], 4, 1),
    "dP1":  ([(1, 0), (0, 1), (-1, -1), (0, -1)], 4, 1),
    "dP2":  ([(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)], 5, 2),
    "dP3":  ([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)], 6, 3),
}


@pytest.mark.parametrize("label", sorted(_FIVED))
def test_orbifold_chamber_and_invariants(label):
    """Every toric 5d SCFT must produce its canonical finite chamber -- the n
    fractional branes as hypermultiplets (Closset-Del Zotto) -- with the
    KK/D0 charge delta = (1,...,1) a null vector of the exchange matrix and the
    correct E_n flavor rank.  An empty mutation sequence selects this chamber."""
    verts, n_states, flavor = _FIVED[label]
    out = api.summarize_bps_toric_quiver(verts, "")     # empty -> orbifold chamber
    assert out["chamber"]["stable_count"] == n_states
    assert out["chamber"]["kind"].startswith("orbifold")
    assert len(out["stable_objects"]) == n_states
    # charges are the standard basis, phases strictly decreasing
    dims = [s["dimension_vector"] for s in out["stable_objects"]]
    assert sorted(dims) == sorted([1 if i == j else 0 for j in range(n_states)]
                                  for i in range(n_states))
    ph = [s["central_charge"]["phase"] for s in out["stable_objects"]]
    assert all(ph[i] >= ph[i + 1] for i in range(len(ph) - 1))
    inv = out["invariants"]
    assert inv["kk_charge"] == [1] * n_states
    assert inv["kk_charge_is_null"] is True
    assert inv["flavor_rank"] == flavor


def test_local_p2_has_no_maximal_green_sequence():
    """Local P^2 (the (3,3,3) Markov quiver) is mutation-infinite: no maximal
    green sequence exists, so the orbifold 3-hyper chamber is THE canonical
    finite chamber (Closset-Del Zotto)."""
    from conformalmanifold.bps import find_maximal_green_sequence
    B = [[0, 3, -3], [-3, 0, 3], [3, -3, 0]]
    assert find_maximal_green_sequence(B) is None
    out = api.summarize_bps_toric_quiver([(1, 0), (0, 1), (-1, -1)], "")
    assert out["chamber"]["maximal_green_sequence"] is None
    assert out["chamber"]["stable_count"] == 3


def test_maximal_green_sequence_finder_a2():
    from conformalmanifold.bps import find_maximal_green_sequence
    assert find_maximal_green_sequence([[0, 1], [-1, 0]]) == [0, 1, 0]


def test_spectral_network_description_from_dimer():
    """The GMN spectral-network reading: genus = Coulomb rank, the asymptotic
    (p,q) legs = the toric polygon edges = the dimer zig-zags."""
    out = api.summarize_bps_toric_quiver([(1, 0), (0, 1), (-1, -1)], "")   # dP0
    sn = out["spectral_network"]
    assert sn["coulomb_rank"] == 1                 # local P^2 has rank 1
    assert sn["num_asymptotic_legs"] == 3          # triangle: 3 legs
    assert sn["flavor_rank"] == 0
    assert "Gaiotto" in sn["framework"]
    legs = {tuple(l["pq_leg"]) for l in sn["asymptotic_legs"]}
    assert len(legs) == 3
    # 5-brane charge conservation: sum of (length-weighted) (p,q) legs vanishes
    def leg_sum(spec):
        sx = sum(l["length"] * l["pq_leg"][0] for l in spec["asymptotic_legs"])
        sy = sum(l["length"] * l["pq_leg"][1] for l in spec["asymptotic_legs"])
        return (sx, sy)
    assert leg_sum(sn) == (0, 0)
    # F0 = P1xP1: 4 legs, rank 1, flavor 1, legs still balance
    snf = api.summarize_bps_toric_quiver([(-1, 0), (1, 0), (0, 1), (0, -1)], "")["spectral_network"]
    assert snf["num_asymptotic_legs"] == 4 and snf["coulomb_rank"] == 1
    assert snf["flavor_rank"] == 1 and leg_sum(snf) == (0, 0)
