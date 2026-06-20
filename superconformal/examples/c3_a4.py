"""C^3/A4: the tetrahedral non-abelian orbifold (McKay quiver).

    python examples/c3_a4.py

The tetrahedral group A4 (order 12) has irreducible representations of
dimension 1, 1, 1, 3, so the McKay quiver has a central 3N node (the 3-dim
irrep) carrying two adjoints u, v, and three outer N nodes (the 1-dim irreps),
with a cubic superpotential
    W = u(Aa + Bb + Cc) + v(Aa + Bb + Cc) + u^3 + v^3.
a-maximization puts every field at R = 2/3; the central charges come out as
a = 3N^2 - 19/24 and c = 3N^2 - 7/12, matching the published values in
arXiv:2310.15792 (Fang-Feng-Xie).
"""

from sqcdkit import scft_observables


def c3_a4(N):
    """McKay quiver of A4 with N regular branes: ranks (3N, N, N, N)."""
    pairs = [("A", "a", 1), ("B", "b", 2), ("C", "c", 3)]
    arrows = [
        {"label": "u", "source": 0, "target": 0, "r_charge": "2/3"},
        {"label": "v", "source": 0, "target": 0, "r_charge": "2/3"},
    ]
    W = [["u", "u", "u"], ["v", "v", "v"]]
    for up, lo, node in pairs:
        arrows.append({"label": up, "source": node, "target": 0, "r_charge": "2/3"})
        arrows.append({"label": lo, "source": 0, "target": node, "r_charge": "2/3"})
        W += [["u", up, lo], ["v", up, lo]]
    return {
        "name": f"C3_A4_N{N}",
        "node_labels": ["n0", "n1", "n2", "n3"],
        "ranks": [3 * N, N, N, N],
        "arrows": arrows,
        "superpotential": [{"coefficient": "1", "factors": f} for f in W],
    }


if __name__ == "__main__":
    print("C^3/A4 (tetrahedral non-abelian orbifold)")
    for N in (1, 2):
        obs = scft_observables(c3_a4(N))
        print(f"  N={N}, ranks={c3_a4(N)['ranks']}:  R(all)=2/3,  "
              f"a={obs['a']},  c={obs['c']},  a/c={obs['hofman_maldacena']['a_over_c']:.4f}")
    print("  (a = 3N^2 - 19/24,  c = 3N^2 - 7/12;  all fields at R = 2/3)")
