"""The conifold: Klebanov-Witten SU(2) x SU(2) with the quartic superpotential.

    python examples/conifold.py

The superconformal R is the symmetric R(A)=R(B)=1/2; the lowest-order index
counts the four mesons Tr(A_i B_j) and the six SU(2) baryons.
"""

from sqcdkit import scft_observables, index_series, format_index

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


if __name__ == "__main__":
    obs = scft_observables(CONIFOLD)
    print("Conifold  SU(2) x SU(2)  (Klebanov-Witten)")
    print(f"  R(A) = {obs['r_charges']['A1']}")
    print(f"  a, c = {obs['a']}, {obs['c']}    a/c = {obs['hofman_maldacena']['a_over_c']:.4f}")
    I = index_series(CONIFOLD, 4)
    print(f"  index = {format_index(I)}    (u^2 = 10: 4 mesons + 6 baryons)")
