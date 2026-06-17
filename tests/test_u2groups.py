"""Validation for the finite U(2) subgroups of SU(3) and the SQLite sweep DB."""

import os
import sqlite3
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conformalmanifold.groups import cyclic  # noqa: E402
from conformalmanifold.pipeline import run  # noqa: E402
from conformalmanifold.u2groups import (  # noqa: E402
    binary_dihedral, binary_polyhedral, embed_u2, u2_library)
from conformalmanifold.database import (  # noqa: E402
    build_database, default_entries, record_for)


def test_embed_u2_is_su3():
    """diag(g, 1/det g) is unitary with determinant 1 for any U(2) g."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        # random U(2): exp(i*Hermitian)
        H = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        H = H + H.conj().T
        g = np.linalg.matrix_power(np.eye(2) + 1j * H / 8, 1)  # near-unitary seed
        # orthonormalize columns -> exact U(2)
        q, _ = np.linalg.qr(g)
        M = embed_u2(q)
        assert np.allclose(M.conj().T @ M, np.eye(3), atol=1e-9)
        assert abs(np.linalg.det(M) - 1) < 1e-9


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_binary_dihedral_order(m):
    assert binary_dihedral(m).order == 4 * m


@pytest.mark.parametrize("kind,order", [("2T", 24), ("2O", 48), ("2I", 120)])
def test_binary_polyhedral_order(kind, order):
    assert binary_polyhedral(kind).order == order


def test_su2_embedding_fixes_third_coordinate():
    """A pure-SU(2) subgroup (det=1) fixes the 3rd coordinate, giving the
    clean sanity value dim_conf = |Gamma| + 1 (C^2/Gamma x C)."""
    for kind, order in [("2T", 24), ("2O", 48), ("2I", 120)]:
        r = run(binary_polyhedral(kind), verbose=False)
        assert r["dim_conf"] == order + 1


def test_cyclic_u2_matches_existing_cyclic():
    """diag(w,1) phase on C^2 = the existing Z_n(1,0,-1) cyclic orbifold."""
    from conformalmanifold.groups import closure
    n = 6
    w = np.exp(2j * np.pi / n)
    G = closure([embed_u2(np.diag([w, 1.0]))])
    a = run(G, verbose=False)["dim_conf"]
    b = run(cyclic(n, (1, 0, n - 1)), verbose=False)["dim_conf"]
    assert a == b


def test_phase_extension_is_genuine_u2():
    """The library's phase-extended groups must NOT fix the 3rd coordinate
    (they are genuine U(2), not SU(2))."""
    lib = u2_library()
    for name in ("BD_2.Z3", "2T.Z3", "2O.Z3", "2I.Z5"):
        rec = record_for(lib[name](), "u2-subgroup")
        assert rec["third_coord_fixed"] == 0, name
        assert rec["connected"] == 1


def test_build_database(tmp_path):
    db = str(tmp_path / "q.db")
    n = build_database(default_entries(), db)
    assert n >= 20
    conn = sqlite3.connect(db)
    try:
        (count,) = conn.execute("SELECT COUNT(*) FROM quivers").fetchone()
        assert count == n
        # every row has a computed dim_conf and a non-empty generators blob
        bad = conn.execute(
            "SELECT name FROM quivers WHERE dim_conf IS NULL OR generators IS NULL"
        ).fetchall()
        assert not bad, bad
        # idempotent: rebuild upserts, no duplication
        build_database(default_entries(), db)
        (count2,) = conn.execute("SELECT COUNT(*) FROM quivers").fetchone()
        assert count2 == count
    finally:
        conn.close()
