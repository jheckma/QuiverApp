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


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "quivers.db"
    n = build_database(default_entries(), path)
    print(f"wrote {n} rows to {path}")
