"""Irrational-R index: the toric quivers whose a-maximized R is irrational.

    python examples/irrational_index.py

index_series stores integer u-powers and so needs rational R; index_symbolic
returns the index at the irrational a-max R as a symbolic series in
tau = (pq)^(1/2), with irrational exponents.
"""

import sympy as sp

from sqcdkit import index_symbolic

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
    "superpotential": [],
}


if __name__ == "__main__":
    # Consistency: at the conifold's RATIONAL a-max R = 1/2, the symbolic index
    # is an ordinary tau series (equal to index_series).
    print("conifold, rational R = 1/2:")
    print("   ", index_symbolic(CONIFOLD, 4, derive_r="amax"))

    # An injected sqrt(13) target R (as for dP1) gives genuine irrational
    # exponents -- the form index_series cannot hold.
    a = (sp.sqrt(13) - 3) / 12
    Rirr = {
        "A1": sp.Rational(1, 2) - 3 * a,
        "A2": sp.Rational(1, 2) + a,
        "B1": sp.Rational(1, 2) + a,
        "B2": sp.Rational(1, 2) + a,
    }
    print("conifold, irrational R (sqrt(13)):")
    print("   ", index_symbolic(CONIFOLD, 2, r_charges=Rirr))
    print()
    print("A genuine toric quiver (dP1: R ~ sqrt(13)) works the same way via")
    print("index_symbolic(theory, order); its operators sit at high order, so")
    print("the fugacity index is slower there.")
