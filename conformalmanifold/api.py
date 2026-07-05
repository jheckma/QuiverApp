"""JSON-serialisable summaries of the pipeline, for the web frontend."""

from __future__ import annotations

import numpy as np

from .chartable import build_character_table
from .conformal import conformal_manifold_dim
from .groups import MatrixGroup, cyclic, make_group
from .quiver import build_quiver
from .scft import orbifold_scft_json, toric_field_R_charges, toric_scft_json
from .inverse import (inverse_quiver_json, inverse_phases_json,
                      dualize_path_json)


def _c(z) -> dict:
    z = complex(z)
    re = 0.0 if abs(z.real) < 1e-4 else round(z.real, 4)
    im = 0.0 if abs(z.imag) < 1e-4 else round(z.imag, 4)
    return {"re": re, "im": im}


# ---------------------------------------------------------------------------
# Cache for the expensive hull-only toric computations (the inverse brane-tiling
# reconstruction in particular can take many seconds, and it depends ONLY on the
# polygon -- NOT on the triangulation/flop/blow-up the user is interacting with).
# So flopping an edge, toggling blow-up mode, or re-selecting the same diagram
# reuses these instead of recomputing them every click.
# ---------------------------------------------------------------------------
_HULL_CACHE: dict = {}
_HULL_CACHE_MAX = 128

# Largest diagram (2*area = gauge-group count) for which the toric tab runs the
# inverse brane-tiling reconstruction inline.  Above this the geometry still
# renders instantly; only the (expensive) reconstructed-quiver panel is skipped.
INVERSE_MAX_AREA2 = 16


def _hull_cached(hull, tag, fn):
    """Memoise `fn()` keyed by (hull, tag); `fn` must depend only on the hull."""
    key = (tuple(map(tuple, hull)), tag)
    if key not in _HULL_CACHE:
        if len(_HULL_CACHE) >= _HULL_CACHE_MAX:
            _HULL_CACHE.clear()               # crude but bounded; dev tool
        _HULL_CACHE[key] = fn()
    return _HULL_CACHE[key]


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
        "scft": orbifold_scft_json(group, quiver),
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
def summarize_toric_web(points, triangulation=None, flop_edge=None,
                        surface_blowup=None, surface_blowdown=None, active=None,
                        dualize=None, dualize_phase=0,
                        include_inverse=True) -> dict:
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
                      before rendering.
    Blow-ups come in TWO distinct kinds, kept separate throughout:

    `surface_blowup`  : optional [x, y] -- blow up a POINT OF THE BASE SURFACE,
                      adding a NEW exceptional P^1 (the polygon GROWS; the
                      returned `diagram.input_points` is the new point set and
                      the triangulation/flop are recomputed from scratch).
                      dP0 -> dP1 -> dP2 -> dP3.
    `surface_blowdown`: optional [x, y] -- contract a -1-curve of the base
                      surface, the inverse of `surface_blowup`.
    `active`        : RESOLUTION OF THE SINGULARITY (fixed polygon): the list of
                      [x, y] lattice points = exceptional divisors currently
                      blown up in the crepant resolution of the singular cone
                      (hull corners are always active).  None = all lattice
                      points = fully resolved (the default phase).  Cells of
                      2*area > 1 in the resulting subdivision are reported as
                      `residual_singularities`.
    `dualize`       : optional list of gauge-face indices -- an INTERACTIVE
                      Seiberg-duality (urban renewal) path.  Every square face
                      (N_f = 2N_c node) of the displayed tiling is a legal
                      move; the moves are applied in order starting from the
                      seed (or from toric phase `dualize_phase`), and the
                      returned `inverse_quiver` is the dualized tiling (with
                      `dual_path` echoed and fresh `square_faces`).
    `include_inverse`: when False, skip the expensive brane-tiling reconstruction
                      (inverse quiver + Seiberg phases + per-field R-charges) so
                      the geometry returns instantly; the UI fetches those
                      separately.  Everything else is always computed."""
    from . import toric as T
    from . import resolution as Rz
    from . import fived as F

    pts = [(int(round(x)), int(round(y))) for (x, y) in points]
    if surface_blowup is not None:
        pts = Rz.surface_blowup(pts, surface_blowup)   # new P^1; ValueError if illegal
        triangulation = None                  # the diagram changed; start fresh
        flop_edge = None
        active = None                         # new diagram starts fully resolved
    elif surface_blowdown is not None:
        pts = Rz.surface_blowdown(pts, surface_blowdown)  # contract a -1-curve
        triangulation = None
        flop_edge = None
        active = None
    hull = T.convex_hull(pts)
    if len(hull) < 3:
        raise ValueError("need at least 3 non-collinear lattice points to span a "
                         "two-dimensional toric diagram")

    area2, B, I, edges = T.polygon_signature(hull)
    legs = T.pq_web(hull)

    # --- resolution: lattice points + triangulation + dual web ---------------
    lat = Rz.lattice_points(hull)             # canonical index order
    interior = set(Rz.interior_lattice_points(hull))

    # partial resolution: which exceptional divisors are blown up.  Corners are
    # forced; None (or the full set) = fully resolved.
    if active is not None:
        activeset = ({(int(round(x)), int(round(y))) for (x, y) in active}
                     & set(lat)) | set(hull)
        if activeset == set(lat):
            active = None                     # everything on = fully resolved
    if active is None:
        activeset = set(lat)

    tri = [tuple(t) for t in triangulation] if triangulation else None
    if not (tri and Rz.is_valid_subdivision(lat, tri, hull, activeset)):
        _, tri = Rz.triangulate(hull, activeset)   # (re)compute the default phase
        flop_edge = None                      # a stale flop no longer applies
    if flop_edge is not None:
        tri = Rz.flop(lat, tri, flop_edge)    # raises ValueError if not flippable

    web = Rz.dual_web(lat, tri, hull)
    flippable = Rz.flippable_edges(lat, tri)

    # residual orbifold singularities of the partial resolution (cells of
    # 2*area > 1), each identified against the named-geometry library.
    residual = Rz.residual_cells(lat, tri)
    for cell in residual:
        geom = T.identify_toric([tuple(v) for v in cell["vertices"]])
        # every simplicial cell is an abelian orbifold C^3/Gamma, |Gamma|=2*area
        cell["label"] = (geom.label if geom is not None
                         else f"C³/Γ orbifold, |Γ| = {cell['area2']}")

    out = {
        "diagram": {
            "input_points": [list(p) for p in pts],
            "hull": [list(p) for p in hull],
            "num_corners": len(hull),
            "boundary_points": B,
            "interior_points": I,
            "norm_area2": area2,
            "edge_lengths": list(edges),
            # BASE-SURFACE sites (polygon changes): blow up a point of the
            # surface = add a NEW exceptional P^1 (W = A+B-O; dP0 -> dP1), or
            # blow down a -1-curve.  Distinct from resolving the singularity
            # (the `active` divisor set at fixed polygon, under "resolution").
            "surface_blowup_candidates": Rz.surface_blowup_candidates(pts),
            "surface_blowdown_candidates": Rz.surface_blowdown_candidates(pts),
        },
        "resolution": {
            "lattice_points": [list(p) for p in lat],
            "boundary_flags": [p not in interior for p in lat],   # True = boundary
            # partial resolution state: which exceptional divisors are blown up
            "active_flags": [p in activeset for p in lat],
            "fully_resolved": len(activeset) == len(lat),
            "num_active": len(activeset),
            "residual_singularities": residual,   # non-unimodular cells + labels
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
        "scft": _hull_cached(hull, "scft", lambda: toric_scft_json(hull)),
    }

    # AdS6/CFT5 reading of the same diagram (kept separate from the 4d quiver).
    # Cheap (pure polygon combinatorics), so it is computed on BOTH the fast
    # include_inverse=False response and the full one.
    factors = F.one_form_symmetry(hull)
    out["fived"] = {
        "rank": I,                          # Coulomb-branch dim = # interior points
        "flavor_rank": F.flavor_rank(B),    # flavor-symmetry rank = # mass params = B - 3
        "one_form_factors": factors,        # invariant factors, e.g. [] / [3] / [2,2]
        "one_form_label": F.abelian_label(factors),
        "note": F.non_isolated_note(edges), # '' unless the singularity is non-isolated
        # full defect group D = Gamma_e (+) Gamma_m, canonical pairing, # global forms
        "defect_group": F.defect_group(hull),
    }

    # Inverse algorithm: reconstruct a quiver + brane tiling from the diagram
    # alone.  This is the ONLY expensive part of the toric tab -- the diagram,
    # triangulation, (p,q) web, gauge-group count, dimension and blow-up/down
    # candidates above are all cheap and always computed.  The reconstruction
    # (Gulotta dimer + Kasteleyn certificate) grows fast with the area, and on
    # some diagrams the dimer search even fails *slowly*, so the UI fetches it in
    # a SEPARATE request (include_inverse=False) -- the geometry renders instantly
    # and blow-up/down never wait on it.  Also size-gated + cached by hull.
    if not include_inverse:
        out["inverse_quiver"] = {"available": False, "deferred": True,
                                 "reason": "fetched separately"}
        out["inverse_phases"] = {"available": False, "deferred": True}
        return _finish_identify(out, hull, T)

    if dualize:
        # interactive Seiberg-duality path: not cached (path-dependent), but a
        # single urban-renewal move is fast.
        out["inverse_quiver"] = dualize_path_json(hull, list(dualize),
                                                  start_phase=dualize_phase)
    elif area2 <= INVERSE_MAX_AREA2:
        out["inverse_quiver"] = _hull_cached(
            hull, "inv", lambda: inverse_quiver_json(hull, max_gauge=INVERSE_MAX_AREA2))
    else:
        out["inverse_quiver"] = {
            "available": False,
            "reason": f"quiver reconstruction skipped for 2*area = {area2} "
                      f"(> {INVERSE_MAX_AREA2}) to keep the builder responsive; "
                      "the diagram, (p,q) web, gauge-group count and dimension "
                      "are still exact.",
        }

    # Seiberg-dual *toric phases* (cycle them in the UI).  The phase search is
    # combinatorial in the face degrees, so only run it for small diagrams; for
    # larger ones the single seed phase above stands.
    if area2 <= 6:
        out["inverse_phases"] = _hull_cached(
            hull, "phases", lambda: inverse_phases_json(hull))
    else:
        out["inverse_phases"] = {
            "available": False,
            "reason": "phase enumeration limited to small diagrams "
                      "(2*area <= 6); showing the seed phase only",
        }

    # per-field superconformal R-charges: needs both the corner R-charges (scft)
    # and the reconstructed dimer (zig-zag legs + superpotential).
    inv = out["inverse_quiver"]
    if inv.get("available") and out["scft"].get("corner_R") and \
            inv.get("fields") and all("zigzag" in f for f in inv["fields"]):
        # needs the seed's zig-zag annotations; a dualized (urban-renewal)
        # tiling has new fields with no zig-zag data, and a NON-TORIC Seiberg
        # dual has no dimer fields at all -- skip in both cases.
        out["scft"]["field_R"] = toric_field_R_charges(
            hull, inv["fields"], inv["superpotential"], out["scft"]["corner_R"])

    return _finish_identify(out, hull, T)


def _finish_identify(out, hull, T):
    """Attach the named-library identification (+ its quiver, if any) to `out`.
    Cheap (a lookup), so it runs on both the fast and full toric responses."""
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
