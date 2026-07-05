"""Build a SQLite database of C^3/Gamma quiver SCFTs, one row per group.

Reports ALL data the pipeline produces (no filtering / ranking here -- that is
left to downstream sweep code).  Columns cover the group, the McKay quiver, and
the conformal-manifold result; matrices and lists are stored as JSON text.

    from conformalmanifold.database import build_database, default_entries
    build_database(default_entries(), "quivers.db")      # create / upsert

    python -m conformalmanifold.database quivers.db      # CLI: build default sweep

Idempotent: re-running upserts on the `name` primary key, so a sweep can be
extended and rebuilt without duplicating rows.
"""
from __future__ import annotations

import json
import sqlite3
import sys

import numpy as np

from . import fived
from .groups import MatrixGroup, library
from .pipeline import run as _run
from .u2groups import u2_library


SCHEMA = """
CREATE TABLE IF NOT EXISTS quivers (
    name              TEXT PRIMARY KEY,
    family            TEXT,
    description       TEXT,
    grp_order         INTEGER,
    is_abelian        INTEGER,   -- 0/1
    third_coord_fixed INTEGER,   -- 0/1: genuine SU(2) embedding (C^2/Gamma x C)
    num_nodes         INTEGER,   -- # irreps = # gauge nodes
    irrep_dims        TEXT,      -- JSON list[int]
    num_arrows        INTEGER,   -- total bifundamentals
    connected         INTEGER,   -- 0/1
    num_cubic_terms   INTEGER,   -- Tr(A^3)
    fixed_sum         INTEGER,   -- S = sum_g fix_Q(g)
    per_direction     TEXT,      -- JSON list[int] or null (abelian only)
    dim_conf          INTEGER,   -- dim_C M_conf = S - 1
    note              TEXT,
    generators        TEXT       -- JSON list of 3x3 complex matrices [[ [re,im], ...]]
);
"""

_COLS = ["name", "family", "description", "grp_order", "is_abelian",
         "third_coord_fixed", "num_nodes", "irrep_dims", "num_arrows",
         "connected", "num_cubic_terms", "fixed_sum", "per_direction",
         "dim_conf", "note", "generators"]


def _gens_json(gens) -> str:
    out = []
    for m in gens:
        out.append([[[float(z.real), float(z.imag)] for z in row] for row in m])
    return json.dumps(out)


def _third_coord_fixed(group: MatrixGroup, tol=1e-9) -> bool:
    """True iff every element fixes the 3rd coordinate (M[2,2]=1) -- i.e. the
    group sits in the SU(2) block (C^2/Gamma x C), not genuine U(2)."""
    return all(abs(g[2, 2] - 1.0) < tol for g in group.elements)


def record_for(group: MatrixGroup, family: str = "") -> dict:
    """Run the pipeline and assemble the full per-group record."""
    res = _run(group, verbose=False)
    q, r = res["quiver"], res["result"]
    return {
        "name": group.name,
        "family": family,
        "description": group.description,
        "grp_order": group.order,
        "is_abelian": int(group.is_abelian()),
        "third_coord_fixed": int(_third_coord_fixed(group)),
        "num_nodes": q.num_nodes,
        "irrep_dims": json.dumps([int(d) for d in res["table"].dims]),
        "num_arrows": q.num_arrows,
        "connected": int(q.is_connected()),
        "num_cubic_terms": q.num_cubic_terms(),
        "fixed_sum": r.fixed_sum,
        "per_direction": json.dumps(r.per_direction) if r.per_direction is not None else None,
        "dim_conf": r.dim_conf,
        "note": r.note,
        "generators": _gens_json(group.generators),
    }


def build_database(entries, db_path: str):
    """entries: iterable of (MatrixGroup, family) or bare MatrixGroup.
    Creates db_path if absent and upserts one row per group."""
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(SCHEMA)
        placeholders = ",".join("?" for _ in _COLS)
        updates = ",".join(f"{c}=excluded.{c}" for c in _COLS if c != "name")
        sql = (f"INSERT INTO quivers ({','.join(_COLS)}) VALUES ({placeholders}) "
               f"ON CONFLICT(name) DO UPDATE SET {updates}")
        n = 0
        for e in entries:
            group, family = e if isinstance(e, tuple) else (e, "")
            rec = record_for(group, family)
            conn.execute(sql, [rec[c] for c in _COLS])
            n += 1
        conn.commit()
    finally:
        conn.close()
    return n


def default_entries():
    """The built-in C^3/Gamma library + the finite-U(2)-subgroup library."""
    out = []
    for name, fac in library().items():
        g = fac()
        g.name = name
        out.append((g, "su3-library"))
    for name, fac in u2_library().items():
        out.append((fac(), "u2-subgroup"))
    return out


# ===========================================================================
# Toric (non-orbifold) Calabi-Yau quivers -> a parallel `toric_quivers` table
# ===========================================================================
TORIC_SCHEMA = """
CREATE TABLE IF NOT EXISTS toric_quivers (
    label           TEXT PRIMARY KEY,
    family          TEXT,
    description     TEXT,
    num_nodes       INTEGER,   -- gauge nodes
    num_fields      INTEGER,   -- bifundamental + adjoint chirals
    num_w_terms     INTEGER,   -- toric (two-term) superpotential monomials
    dim_conf        INTEGER,   -- dim_C M_conf via LS/NSVZ counting (authoritative)
    dim_conf_geom   INTEGER,   -- B - 1 from the toric diagram (cross-check)
    boundary_points INTEGER,   -- B = # lattice points on the toric-diagram boundary
    interior_points INTEGER,   -- I = # interior lattice points
    norm_area2      INTEGER,   -- 2 * area = # gauge nodes of the minimal tiling
    edge_lengths    TEXT,      -- JSON: sorted primitive edge lengths
    diagram         TEXT,      -- JSON: CCW lattice-polygon vertices
    arrows          TEXT,      -- JSON: {field: [src, tgt]}
    superpotential  TEXT,      -- JSON: [[sign, [field,...]], ...]
    rank_5d         INTEGER,   -- 5d SCFT rank (Coulomb-branch dim) = I
    flavor_rank_5d  INTEGER,   -- 5d flavor-symmetry rank = B - 3
    one_form_5d     TEXT,      -- 5d 1-form symmetry, e.g. 'Z_3' / 'trivial'
    defect_group_5d TEXT,      -- D = Gamma (+) Gamma, e.g. 'Z_3 (+) Z_3'
    pairing_5d      TEXT,      -- canonical Dirac pairing, e.g. '1/3' (isolated)
    n_global_forms_5d INTEGER  -- # polarizations of D (isolated; else NULL)
);
"""

_TORIC_COLS = ["label", "family", "description", "num_nodes", "num_fields",
               "num_w_terms", "dim_conf", "dim_conf_geom", "boundary_points",
               "interior_points", "norm_area2", "edge_lengths", "diagram",
               "arrows", "superpotential", "rank_5d", "flavor_rank_5d",
               "one_form_5d", "defect_group_5d", "pairing_5d",
               "n_global_forms_5d"]

# columns added after the first release: ALTERed into pre-existing db files
_TORIC_5D_COLS = {"rank_5d": "INTEGER", "flavor_rank_5d": "INTEGER",
                  "one_form_5d": "TEXT", "defect_group_5d": "TEXT",
                  "pairing_5d": "TEXT", "n_global_forms_5d": "INTEGER"}


def _fived_defect_cols(hull) -> dict:
    """The defect-group columns shared by both toric tables."""
    if not hull:
        return {"defect_group_5d": None, "pairing_5d": None,
                "n_global_forms_5d": None}
    dg = fived.defect_group(hull)
    return {
        "defect_group_5d": dg["label"],
        "pairing_5d": ", ".join(dg["pairing"]) if dg["pairing"] else None,
        "n_global_forms_5d": dg["num_global_forms"],
    }


def _ensure_columns(conn, table: str, cols: dict):
    """ALTER TABLE ADD COLUMN for any of `cols` missing from an existing db
    (CREATE TABLE IF NOT EXISTS does not touch already-created tables)."""
    have = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    for name, sqltype in cols.items():
        if name not in have:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {name} {sqltype}")


def toric_record(q) -> dict:
    """Full per-geometry record for a `toric.ToricQuiver`."""
    from .toric import convex_hull, polygon_signature
    sig = q.signature()
    area2, B, I, edges = sig if sig else (None, None, None, ())
    hull = convex_hull(q.diagram) if q.diagram else None
    return {
        "rank_5d": I,
        "flavor_rank_5d": fived.flavor_rank(B) if B is not None else None,
        "one_form_5d": (fived.abelian_label(fived.one_form_symmetry(hull))
                        if hull else None),
        **_fived_defect_cols(hull),
        "label": q.label,
        "family": q.family,
        "description": q.description,
        "num_nodes": q.num_nodes,
        "num_fields": q.num_fields,
        "num_w_terms": q.num_w_terms,
        "dim_conf": q.dim_conf_ls(),
        "dim_conf_geom": q.dim_conf_geometric(),
        "boundary_points": B,
        "interior_points": I,
        "norm_area2": area2,
        "edge_lengths": json.dumps(list(edges)),
        "diagram": json.dumps([list(v) for v in hull] if hull else []),
        "arrows": json.dumps({str(k): list(v) for k, v in q.arrows.items()}),
        "superpotential": json.dumps([[s, list(c)] for s, c in q.W]),
    }


def build_toric_database(quivers, db_path: str):
    """quivers: iterable of toric.ToricQuiver.  Upserts one row per label."""
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(TORIC_SCHEMA)
        _ensure_columns(conn, "toric_quivers", _TORIC_5D_COLS)
        placeholders = ",".join("?" for _ in _TORIC_COLS)
        updates = ",".join(f"{c}=excluded.{c}" for c in _TORIC_COLS if c != "label")
        sql = (f"INSERT INTO toric_quivers ({','.join(_TORIC_COLS)}) "
               f"VALUES ({placeholders}) "
               f"ON CONFLICT(label) DO UPDATE SET {updates}")
        n = 0
        for q in quivers:
            rec = toric_record(q)
            conn.execute(sql, [rec[c] for c in _TORIC_COLS])
            n += 1
        conn.commit()
    finally:
        conn.close()
    return n


TORIC_DIAGRAM_SCHEMA = """
CREATE TABLE IF NOT EXISTS toric_diagrams (
    label           TEXT PRIMARY KEY,
    family          TEXT,
    description     TEXT,
    num_gauge       INTEGER,   -- # gauge groups = 2 * area of the toric diagram
    dim_conf        INTEGER,   -- dim_C M_conf = B - 1
    boundary_points INTEGER,   -- B
    interior_points INTEGER,   -- I
    norm_area2      INTEGER,   -- 2 * area
    edge_lengths    TEXT,      -- JSON: sorted primitive edge lengths
    diagram         TEXT,      -- JSON: CCW lattice-polygon vertices
    note            TEXT,
    rank_5d         INTEGER,   -- 5d SCFT rank (Coulomb-branch dim) = I
    flavor_rank_5d  INTEGER,   -- 5d flavor-symmetry rank = B - 3
    one_form_5d     TEXT,      -- 5d 1-form symmetry, e.g. 'Z_3' / 'trivial'
    defect_group_5d TEXT,      -- D = Gamma (+) Gamma, e.g. 'Z_3 (+) Z_3'
    pairing_5d      TEXT,      -- canonical Dirac pairing, e.g. '1/3' (isolated)
    n_global_forms_5d INTEGER  -- # polarizations of D (isolated; else NULL)
);
"""

_TORIC_DIAGRAM_COLS = ["label", "family", "description", "num_gauge", "dim_conf",
                       "boundary_points", "interior_points", "norm_area2",
                       "edge_lengths", "diagram", "note", "rank_5d",
                       "flavor_rank_5d", "one_form_5d", "defect_group_5d",
                       "pairing_5d", "n_global_forms_5d"]


def toric_diagram_record(d) -> dict:
    """Per-geometry record for a `toric.ToricDiagram`."""
    area2, B, I, edges = d.signature()
    hull = d.hull()
    return {
        "rank_5d": I,
        "flavor_rank_5d": fived.flavor_rank(B),
        "one_form_5d": fived.abelian_label(fived.one_form_symmetry(hull)),
        **_fived_defect_cols(hull),
        "label": d.label,
        "family": d.family,
        "description": d.description,
        "num_gauge": d.num_gauge_groups,
        "dim_conf": d.dim_conf(),
        "boundary_points": B,
        "interior_points": I,
        "norm_area2": area2,
        "edge_lengths": json.dumps(list(edges)),
        "diagram": json.dumps([list(v) for v in hull]),
        "note": d.note,
    }


def build_toric_diagram_database(diagrams, db_path: str):
    """diagrams: iterable of toric.ToricDiagram.  Upserts one row per label."""
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(TORIC_DIAGRAM_SCHEMA)
        _ensure_columns(conn, "toric_diagrams", _TORIC_5D_COLS)
        placeholders = ",".join("?" for _ in _TORIC_DIAGRAM_COLS)
        updates = ",".join(f"{c}=excluded.{c}" for c in _TORIC_DIAGRAM_COLS
                           if c != "label")
        sql = (f"INSERT INTO toric_diagrams ({','.join(_TORIC_DIAGRAM_COLS)}) "
               f"VALUES ({placeholders}) "
               f"ON CONFLICT(label) DO UPDATE SET {updates}")
        n = 0
        for d in diagrams:
            rec = toric_diagram_record(d)
            conn.execute(sql, [rec[c] for c in _TORIC_DIAGRAM_COLS])
            n += 1
        conn.commit()
    finally:
        conn.close()
    return n


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "quivers.db"
    n = build_database(default_entries(), path)
    from .toric import default_toric_library, default_toric_diagram_library
    m = build_toric_database(default_toric_library(), path)
    k = build_toric_diagram_database(default_toric_diagram_library(), path)
    print(f"wrote {n} orbifold rows + {m} toric-quiver rows + "
          f"{k} toric-diagram rows to {path}")
