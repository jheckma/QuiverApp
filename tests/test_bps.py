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
