"""Partition combinatorics: transpose, orbit dimension, dominance order,
covering relations, and Kraft-Procesi transverse-slice types (type A).

The Kraft-Procesi fixtures are cross-checked against Brieskorn-Slodowy: the
transverse slice to the subregular orbit of sl(N) is the Kleinian singularity
C^2/Z_N = A_{N-1}.
"""

import pytest

from conformalmanifold import partitions as P


def test_transpose_involution_and_values():
    assert P.transpose((3,)) == (1, 1, 1)
    assert P.transpose((2, 1)) == (2, 1)
    assert P.transpose((2, 2, 1)) == (3, 2)
    for part in [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1), (4, 2, 1)]:
        assert P.transpose(P.transpose(part)) == P.normalize(part)


def test_orbit_dimension_su():
    # sl(3): regular (nilcone) 6, subregular/minimal 4, zero 0.
    assert P.orbit_dim("su", (3,)) == 6
    assert P.orbit_dim("su", (2, 1)) == 4
    assert P.orbit_dim("su", (1, 1, 1)) == 0
    # sl(N) nilcone (regular orbit) has complex dimension N(N-1).
    for N in range(2, 8):
        assert P.orbit_dim("su", (N,)) == N * (N - 1)
    # minimal orbit of sl(N): complex dim 2(N-1).
    for N in range(2, 8):
        minimal = (2,) + (1,) * (N - 2)
        assert P.orbit_dim("su", minimal) == 2 * (N - 1)


def test_orbit_dimension_classical_parity():
    # so/sp orbit dimensions are even (quaternionic dim is an integer).
    for part in P.orbits("so", 6):
        assert P.orbit_dim("so", part) % 2 == 0
    for part in P.orbits("sp", 6):
        assert P.orbit_dim("sp", part) % 2 == 0


def test_valid_partition_constraints():
    # so: even parts need even multiplicity.
    assert P.valid_partition("so", (3, 3, 1)) is True
    assert P.valid_partition("so", (2, 1)) is False        # single 2 (even, odd mult)
    assert P.valid_partition("so", (2, 2, 1)) is True
    # sp: odd parts need even multiplicity.
    assert P.valid_partition("sp", (2, 1, 1)) is True
    assert P.valid_partition("sp", (3, 1)) is False         # 3 and 1 each odd mult
    assert P.valid_partition("sp", (2, 2)) is True
    # su: anything goes.
    assert P.valid_partition("su", (3, 1)) is True


def test_dominance_order():
    assert P.dominates((3,), (2, 1))
    assert P.dominates((2, 1), (1, 1, 1))
    assert P.strictly_dominates((3,), (1, 1, 1))
    assert not P.dominates((2, 1), (3,))
    # incomparable pair in partitions of 6.
    assert not P.dominates((3, 3), (4, 1, 1))
    assert not P.dominates((4, 1, 1), (3, 3))
    with pytest.raises(ValueError):
        P.dominates((3,), (2, 2))       # different totals


def test_covering_relations_su3_and_su4_are_chains():
    assert P.covering_relations("su", 3) == [((3,), (2, 1)), ((2, 1), (1, 1, 1))]
    # partitions of 4 form a chain in dominance order.
    chain = [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
    rels = P.covering_relations("su", 4)
    assert rels == list(zip(chain, chain[1:]))


def test_dominance_branches_at_n6():
    # (3,3) and (4,1,1) are incomparable -> the dominance poset is not a chain.
    assert not P.dominates((3, 3), (4, 1, 1))
    assert not P.dominates((4, 1, 1), (3, 3))
    rels = P.covering_relations("su", 6)
    out_degree = {}
    for a, _b in rels:
        out_degree[a] = out_degree.get(a, 0) + 1
    # some orbit degenerates in more than one minimal way -> the Hasse diagram forks.
    assert max(out_degree.values()) >= 2


def test_kp_slice_types_su3():
    top = P.kp_slice_type((3,), (2, 1))
    assert (top["letter"], top["k"], top["dim_H"]) == ("A", 2, 1)   # C^2/Z_3
    bot = P.kp_slice_type((2, 1), (1, 1, 1))
    assert (bot["letter"], bot["k"], bot["dim_H"]) == ("a", 2, 2)   # min orbit sl(3)


def test_kp_slice_subregular_is_kleinian_A_Nminus1():
    # Brieskorn-Slodowy: sl(N) subregular transverse slice = C^2/Z_N = A_{N-1}.
    for N in range(2, 8):
        subreg = (N - 1, 1) if N > 1 else (1,)
        s = P.kp_slice_type((N,), subreg)
        assert s["letter"] == "A" and s["k"] == N - 1 and s["dim_H"] == 1


def test_kp_dim_matches_orbit_dim_drop():
    # Each covering's transverse dim equals the drop in quaternionic orbit dim;
    # in particular every Kleinian ('A') slice must genuinely drop dim_H by 1.
    for N in range(2, 8):
        for a, b in P.covering_relations("su", N):
            s = P.kp_slice_type(a, b)
            drop = (P.orbit_dim("su", a) - P.orbit_dim("su", b)) // 2
            assert s["dim_H"] == drop
            if s["letter"] == "A":
                assert drop == 1


def test_kp_chain_dims_telescope_along_a_maximal_chain():
    # Walk one maximal chain from the nilcone to zero; transverse dims sum to dim_H.
    for N in range(2, 8):
        rels = P.covering_relations("su", N)
        children = {}
        for a, b in rels:
            children.setdefault(a, []).append(b)
        node, total = (N,), 0
        while node in children:
            nxt = children[node][0]
            total += P.kp_slice_type(node, nxt)["dim_H"]
            node = nxt
        assert node == (1,) * N
        assert total == P.orbit_dim("su", (N,)) // 2


def test_a1_equals_A1():
    # sl(2): (2) -> (1,1) is C^2/Z_2, labeled A_1 (== a_1).
    s = P.kp_slice_type((2,), (1, 1))
    assert s["k"] == 1 and s["dim_H"] == 1
