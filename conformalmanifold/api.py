"""JSON-serialisable summaries of the pipeline, for the web frontend."""

from __future__ import annotations

import numpy as np

from .chartable import build_character_table
from .conformal import conformal_manifold_dim
from .groups import MatrixGroup, cyclic, make_group
from .quiver import build_quiver


def _c(z) -> dict:
    z = complex(z)
    re = 0.0 if abs(z.real) < 1e-4 else round(z.real, 4)
    im = 0.0 if abs(z.imag) < 1e-4 else round(z.imag, 4)
    return {"re": re, "im": im}


def summarize(group: MatrixGroup) -> dict:
    """Run Step 1 -> 2 -> 3 and return a JSON-friendly dict."""
    table = build_character_table(group)
    quiver = build_quiver(group, table)
    result = conformal_manifold_dim(group)

    A = quiver.adjacency
    chiral = bool(np.any(A != A.T))
    q_real = all(abs(complex(z).imag) < 1e-4 for z in table.chi_Q)

    return {
        "group": {
            "name": group.name,
            "description": group.description,
            "order": group.order,
            "abelian": group.is_abelian(),
            "num_generators": len(group.generators),
        },
        "quiver": {
            "num_nodes": quiver.num_nodes,
            "dims": list(quiver.dims),
            "node_labels": [quiver.node_rank_label(i) for i in range(quiver.num_nodes)],
            "adjacency": quiver.adjacency.astype(int).tolist(),
            "num_arrows": quiver.num_arrows,
            "num_cubic": quiver.num_cubic_terms(),
            "connected": quiver.is_connected(),
            "chiral": chiral,
            "Q_real": q_real,
            "net_chiral": (A - A.T).astype(int).tolist(),
            "chi_Q": [_c(z) for z in table.chi_Q],
            "class_sizes": list(table.class_sizes),
            "class_orders": list(table.class_orders),
        },
        "conformal": {
            "fixed_sum": result.fixed_sum,
            "dim_conf": result.dim_conf,
            "per_direction": result.per_direction,
            "note": result.note,
        },
    }


def summarize_named(name: str) -> dict:
    return summarize(make_group(name))


def summarize_cyclic(n: int, a: int, b: int, c: int) -> dict:
    return summarize(cyclic(n, (a, b, c)))


def summarize_freeform(expr: str) -> dict:
    from .freeform import build_from_expr
    return summarize(build_from_expr(expr))
