"""SPP: the suspended pinch point, a non-orbifold toric quiver with an adjoint.

    python examples/spp.py

Three SU(2) nodes, an adjoint X22 at node 2, and a mixed cubic/quartic
superpotential
    W = X12 X21 X22 - X22 X23 X32 + X13 X23 X31 X32 - X12 X13 X21 X31.
a-maximization gives irrational R-charges (involving sqrt(97)); a, c come back
as exact closed forms and no operator decouples.
"""

from sqcdkit import scft_observables

SPP = {
    "name": "SPP",
    "node_labels": ["n1", "n2", "n3"],
    "ranks": [2, 2, 2],
    "arrows": [
        {"label": "X22", "source": 1, "target": 1, "r_charge": "1/2"},  # adjoint at node 2
        {"label": "X12", "source": 0, "target": 1, "r_charge": "1/2"},
        {"label": "X21", "source": 1, "target": 0, "r_charge": "1/2"},
        {"label": "X23", "source": 1, "target": 2, "r_charge": "1/2"},
        {"label": "X32", "source": 2, "target": 1, "r_charge": "1/2"},
        {"label": "X13", "source": 0, "target": 2, "r_charge": "1/2"},
        {"label": "X31", "source": 2, "target": 0, "r_charge": "1/2"},
    ],
    "superpotential": [
        {"coefficient": "1", "factors": ["X12", "X21", "X22"]},
        {"coefficient": "-1", "factors": ["X22", "X23", "X32"]},
        {"coefficient": "1", "factors": ["X13", "X23", "X31", "X32"]},
        {"coefficient": "-1", "factors": ["X12", "X13", "X21", "X31"]},
    ],
}


if __name__ == "__main__":
    obs = scft_observables(SPP)
    print("SPP (suspended pinch point = L^{1,2,1}), SU(2)^3")
    print(f"  R(X22) = {obs['r_charges']['X22']}   (adjoint at node 2)")
    print(f"  R(X12) = {obs['r_charges']['X12']}    R(X13) = {obs['r_charges']['X13']}")
    print(f"  a = {obs['a']}")
    print(f"  c = {obs['c']}")
    print(f"  a/c = {obs['hofman_maldacena']['a_over_c']:.5f}   HM ok: {obs['hofman_maldacena']['ok']}")
    print(f"  mesonic_below_bound: {obs['mesonic_below_bound']}  (no operator decoupling)")
