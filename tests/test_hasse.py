"""Higgs-branch Hasse diagrams from the nilpotent-orbit dominance poset.

Cross-checked against the Kraft-Procesi structure of the sl(N) nilcone: the chain of
transverse slices from the regular orbit down to zero telescopes to the nilcone's
quaternionic dimension.
"""

import pytest

from conformalmanifold import hasse as H
from conformalmanifold import partitions as P
from conformalmanifold import api


def test_sl3_hasse_is_a_chain_with_A2_then_a2():
    d = H.nilpotent_hasse("su", 3)
    # three orbits: zero, minimal (2,1), regular (3).
    labels = {nd["id"]: tuple(nd["partition"]) for nd in d["nodes"]}
    assert sorted(labels.values(), key=len) == [(3,), (2, 1), (1, 1, 1)]
    # two covering edges, labeled A_2 (top) then a_2.
    edge_by_parts = {}
    for e in d["edges"]:
        edge_by_parts[(labels[e["from"]], labels[e["to"]])] = e
    top = edge_by_parts[((3,), (2, 1))]
    bot = edge_by_parts[((2, 1), (1, 1, 1))]
    assert (top["letter"], top["k"], top["dim_H"]) == ("A", 2, 1)
    assert (bot["letter"], bot["k"], bot["dim_H"]) == ("a", 2, 2)


def test_layers_ordered_by_dimension():
    d = H.nilpotent_hasse("su", 4)
    # layers ascend in quaternionic dimension; bottom layer is the zero orbit (dim 0).
    dims_per_layer = []
    id2dim = {nd["id"]: nd["dim_H"] for nd in d["nodes"]}
    for layer in d["layers"]:
        layer_dims = {id2dim[i] for i in layer}
        assert len(layer_dims) == 1               # each layer is one dimension
        dims_per_layer.append(layer_dims.pop())
    assert dims_per_layer == sorted(dims_per_layer)
    assert dims_per_layer[0] == 0


def test_edge_dims_telescope_along_a_chain():
    for N in range(2, 7):
        d = H.nilpotent_hasse("su", N)
        id2part = {nd["id"]: tuple(nd["partition"]) for nd in d["nodes"]}
        children = {}
        for e in d["edges"]:
            children.setdefault(e["from"], []).append(e)
        # start at the regular orbit (max dim_H), walk down.
        top = max(d["nodes"], key=lambda nd: nd["dim_H"])["id"]
        node, total = top, 0
        while node in children:
            e = children[node][0]
            total += e["dim_H"]
            node = e["to"]
        assert id2part[node] == (1,) * N
        assert total == P.orbit_dim("su", (N,)) // 2


def test_su6_hasse_branches():
    d = H.nilpotent_hasse("su", 6)
    out_degree = {}
    for e in d["edges"]:
        out_degree[e["from"]] = out_degree.get(e["from"], 0) + 1
    assert max(out_degree.values()) >= 2          # the diagram forks


def test_orthosymplectic_reports_dims_with_warning():
    d = H.nilpotent_hasse("so", 5)
    assert d["warnings"]                          # slice-type not classified
    for e in d["edges"]:
        assert e["letter"] == "?" and e["dim_H"] >= 1


def test_restricted_to_orbit_closure():
    # The Higgs-branch Hasse of the minimal orbit (2,1,1) of sl(4) is just the
    # two-leaf chain (2,1,1) > (1^4); the larger orbits (4), (3,1), (2,2) drop out.
    d = H.nilpotent_hasse("su", 4, top=(2, 1, 1))
    parts = {tuple(nd["partition"]) for nd in d["nodes"]}
    assert parts == {(2, 1, 1), (1, 1, 1, 1)}
    assert len(d["edges"]) == 1
    assert d["edges"][0]["slice"].startswith("a")     # minimal slice a_3


def test_hasse_json_partition_restricts():
    d = H.hasse_json(kind="su", N=4, partition="2,1,1")
    assert d["top"] == [2, 1, 1]
    assert len(d["nodes"]) == 2


def test_api_summarize_magnetic_hasse():
    d = api.summarize_magnetic_hasse(kind="su", N=3, method="partition")
    assert d["available"] and d["method"] == "partition"
    assert len(d["nodes"]) == 3
    assert len(d["edges"]) == 2


def test_subtraction_recognizes_orbit_quiver():
    # "1-2-[3]" = T[SU(3)] = M((3)); subtraction returns the sl(3) nilcone Hasse.
    d = H.quiver_subtraction_hasse("1-2-[3]")
    assert d["method"] == "subtraction"
    assert d["recognized_as"]["partition"] == [3]
    assert len(d["nodes"]) == 3
    assert [e["slice"][0] for e in d["edges"]] == ["A", "a"]


def test_subtraction_rejects_unrecognized_quiver():
    with pytest.raises(ValueError):
        H.quiver_subtraction_hasse("2-2-2")


def test_api_subtraction_via_recognition():
    # [1]-1-1-[1] is M((2,1)), the minimal orbit of sl(3).
    d = api.summarize_magnetic_hasse(method="subtraction", quiver="[1]-1-1-[1]")
    assert d["available"] and d["recognized_as"]["partition"] == [2, 1]


def test_subtraction_recognizes_sl2_nilcone():
    # a single U(1) with two flavors is M((2)) of sl(2).
    d = H.quiver_subtraction_hasse("[1]-1-[1]")
    assert d["recognized_as"]["partition"] == [2]
