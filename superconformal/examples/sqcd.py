"""SQCD: SU(Nc) with Nf flavors (no superpotential).

    python examples/sqcd.py

Reproduces R_Q = 1 - Nc/Nf, the central charges and Hofman-Maldacena ratio,
the flavor SU(Nf)^3 anomaly, the one-loop b0, and the (s-confining)
superconformal index of SU(2) with three flavors.
"""

import sympy as sp

from sqcdkit import scft_observables, index_series, format_index


def sqcd(Nc, Nf):
    """SU(Nc) gauge node (0) + two SU(Nf) global flavor nodes (1, 2)."""
    return {
        "name": f"sqcd_{Nc}_{Nf}",
        "node_labels": ["g0"],
        "ranks": [Nc],
        "arrows": [
            {"label": "Q", "source": 0, "target": 1, "r_charge": "1/2"},
            {"label": "Qb", "source": 2, "target": 0, "r_charge": "1/2"},
        ],
        "superpotential": [],
    }


def _unrefined(series):
    subs = {s: 1 for c in series.values()
            if hasattr(c, "free_symbols") for s in c.free_symbols}
    return {k: (sp.expand(c.subs(subs)) if hasattr(c, "subs") else c)
            for k, c in series.items()}


if __name__ == "__main__":
    Nc, Nf = 2, 3
    th = sqcd(Nc, Nf)
    obs = scft_observables(th, flavor_ranks=[Nf, Nf])
    print(f"SQCD  SU({Nc})  with  Nf={Nf}")
    print(f"  R_Q  = {obs['r_charges']['Q']}    (= 1 - Nc/Nf)")
    print(f"  a, c = {obs['a']}, {obs['c']}    a/c = {obs['hofman_maldacena']['a_over_c']:.4f}")
    print(f"  flavor SU(Nf)^3 = {[v['SU3'] for v in obs['flavor_anomalies'].values()]}    (= -+Nc)")
    print(f"  one-loop b0     = {obs['one_loop_b0']}    (= 3Nc - Nf)")
    I = index_series(th, 4, flavor_ranks=[Nf, Nf], derive_r="amax")
    print(f"  index = {format_index(_unrefined(I))}    (s-confinement: 15 at R=2/3)")
