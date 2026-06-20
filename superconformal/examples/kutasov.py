"""Kutasov-like theory: SU(2) with an adjoint Phi and Nf flavors, W = Tr Phi^3.

    python examples/kutasov.py

The superpotential W = Tr Phi^(k+1) (here k=2) fixes R(Phi) = 2/(k+1) = 2/3;
a-maximization over the baryonic U(1) then gives R(Q). All R-charges are
rational, so the index is the ordinary index_series.
"""

from sqcdkit import scft_observables

KUTASOV = {
    "name": "kutasov",
    "node_labels": ["g0"],
    "ranks": [2],
    "arrows": [
        {"label": "Phi", "source": 0, "target": 0, "r_charge": "1/2"},
        {"label": "Q", "source": 0, "target": 1, "r_charge": "1/2"},
        {"label": "Qb", "source": 2, "target": 0, "r_charge": "1/2"},
    ],
    "superpotential": [{"coefficient": "1", "factors": ["Phi", "Phi", "Phi"]}],
}


if __name__ == "__main__":
    obs = scft_observables(KUTASOV, flavor_ranks=[3, 3])
    print("Kutasov  SU(2) + adjoint + Nf=3,   W = Tr Phi^3")
    print(f"  R(Phi) = {obs['r_charges']['Phi']}    (= 2/(k+1) = 2/3, fixed by W)")
    print(f"  R(Q)   = {obs['r_charges']['Q']}")
    print(f"  a, c   = {obs['a']}, {obs['c']}    a/c = {obs['hofman_maldacena']['a_over_c']:.4f}")
