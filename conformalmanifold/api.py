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


# ===========================================================================
# Toric "web builder": a user-drawn toric diagram -> dual (p,q) web + quiver
# ===========================================================================
def summarize_toric_web(points, triangulation=None, flop_edge=None) -> dict:
    """Given the lattice points of a toric diagram (the dot/grid diagram a
    physicist draws), return its convex-hull toric diagram, a triangulation
    (resolution / toric phase) of it, the dual (p,q) 5-brane web for that
    triangulation, the conformal-manifold dimension, and -- when the geometry is
    recognised -- the explicit quiver gauge theory.

    `points`        : iterable of (x, y) integer lattice points.
    `triangulation` : optional list of index-triples (into the canonical lattice
                      point order) -- the current phase to render; a default
                      triangulation is computed if omitted or inconsistent.
    `flop_edge`     : optional [i, j] internal edge to flop in `triangulation`
                      before rendering."""
    from . import toric as T
    from . import resolution as Rz

    pts = [(int(round(x)), int(round(y))) for (x, y) in points]
    hull = T.convex_hull(pts)
    if len(hull) < 3:
        raise ValueError("need at least 3 non-collinear lattice points to span a "
                         "two-dimensional toric diagram")

    area2, B, I, edges = T.polygon_signature(hull)
    legs = T.pq_web(hull)

    # --- resolution: lattice points + triangulation + dual web ---------------
    lat = Rz.lattice_points(hull)             # canonical index order
    interior = set(Rz.interior_lattice_points(hull))

    tri = [tuple(t) for t in triangulation] if triangulation else None
    if not (tri and Rz.is_valid_triangulation(lat, tri, hull)):
        _, tri = Rz.triangulate(hull)         # (re)compute the default phase
        flop_edge = None                      # a stale flop no longer applies
    if flop_edge is not None:
        tri = Rz.flop(lat, tri, flop_edge)    # raises ValueError if not flippable

    web = Rz.dual_web(lat, tri, hull)
    flippable = Rz.flippable_edges(lat, tri)

    out = {
        "diagram": {
            "input_points": [list(p) for p in pts],
            "hull": [list(p) for p in hull],
            "num_corners": len(hull),
            "boundary_points": B,
            "interior_points": I,
            "norm_area2": area2,
            "edge_lengths": list(edges),
        },
        "resolution": {
            "lattice_points": [list(p) for p in lat],
            "boundary_flags": [p not in interior for p in lat],   # True = boundary
            "triangulation": [list(t) for t in tri],              # index-triples
            "num_triangles": len(tri),
            "web": web,                       # junctions / internal_edges / external_legs
            "flippable_edges": flippable,     # [[i,j], ...]
        },
        "web": {
            "legs": legs,                      # each {"pq":[p,q], "base":[x,y], "edge":i}
            "num_external_legs": len(legs),
            "charge_sum": [sum(l["pq"][0] for l in legs),
                           sum(l["pq"][1] for l in legs)],
        },
        "conformal": {
            "dim_conf": B - 1,
            "formula": "dim_C M_conf = B - 1   (B = boundary lattice points "
                       "= # external (p,q) legs)",
            "num_gauge_groups": area2,         # = 2 * area of the toric diagram
        },
    }

    geom = T.identify_toric(hull)
    if geom is None:
        out["identified"] = {
            "matched": False,
            "label": None,
            "family": "(general toric CY3)",
            "description": "A valid toric Calabi-Yau three-fold not in the named "
                           "library; invariants are read from the toric diagram.",
            "has_quiver": False,
            "note": "An explicit quiver needs the inverse (brane-tiling) "
                    "algorithm; the dimension and gauge-group count above hold "
                    "for any toric phase.",
        }
        return out

    note = getattr(geom, "note", "")
    if isinstance(geom, T.ToricQuiver):
        N = geom.num_nodes
        out["identified"] = {
            "matched": True,
            "label": geom.label,
            "family": geom.family,
            "description": geom.description,
            "has_quiver": True,
            "note": note,
        }
        out["quiver"] = {
            "num_nodes": N,
            "dims": [1] * N,                   # all toric nodes are U(N)
            "node_labels": ["U(N)"] * N,
            "adjacency": geom.adjacency_matrix(),
            "num_arrows": geom.num_fields,
            "num_w_terms": geom.num_w_terms,
            "dim_conf_ls": geom.dim_conf_ls(),
            "valid": geom.validate() == [],
        }
    else:  # ToricDiagram (named, diagram-only)
        out["identified"] = {
            "matched": True,
            "label": geom.label,
            "family": geom.family,
            "description": geom.description,
            "has_quiver": False,
            "note": (note + "  " if note else "") +
                    "Named geometry; explicit quiver not built in (diagram-only).",
        }
    return out
