"""5d SCFT / AdS6 readings: the 1-form symmetry from the toric diagram.

Checks the 1-form symmetry (1-form part of the defect group) against standard
results from the literature.
"""

from conformalmanifold import api, fived as F, toric as T


def _one_form(points):
    hull = T.convex_hull([(int(x), int(y)) for (x, y) in points])
    return F.one_form_symmetry(hull)


def test_one_form_symmetry_known_geometries():
    # trivial cases
    assert _one_form([(0, 0), (1, 0), (0, 1)]) == []                    # C^3
    assert _one_form([(0, 0), (1, 0), (1, 1), (0, 1)]) == []            # conifold
    assert _one_form([(1, 0), (0, 1), (-1, -1), (0, -1)]) == []         # dP1
    assert _one_form([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]) == []  # dP3
    # nontrivial cases
    assert _one_form([(1, 0), (0, 1), (-1, -1)]) == [3]                 # dP0 = local P^2 (E0) -> Z_3
    assert _one_form([(-1, 0), (1, 0), (0, 1), (0, -1)]) == [2]         # F0 = local P^1xP^1 -> Z_2
    assert _one_form([(0, 0), (4, 0), (0, 1)]) == [4]                   # C^2/Z_4 x C -> Z_4
    assert _one_form([(0, 0), (2, 0), (0, 2)]) == [2, 2]               # C^3/(Z_2xZ_2)
    assert _one_form([(0, 0), (3, 0), (0, 3)]) == [3, 3]               # C^3/(Z_3xZ_3)


def test_one_form_symmetry_is_gl2z_invariant():
    # The group must not depend on a GL(2,Z) (lattice) change of basis of the diagram.
    base = [(1, 0), (0, 1), (-1, -1)]                                   # dP0 -> Z_3
    sheared = [(x + 2 * y, y) for (x, y) in base]                       # apply [[1,2],[0,1]]
    assert _one_form(base) == _one_form(sheared) == [3]


def test_flavor_rank():
    assert F.flavor_rank(3) == 0      # C^3 / dP0 = local P^2 (E_0): no flavor
    assert F.flavor_rank(4) == 1      # local F0 / dP1 (E_1): rank-1 flavor
    assert F.flavor_rank(6) == 3      # local dP3 (E_3): rank-3 flavor
    assert F.flavor_rank(2) == 0      # degenerate guard: never negative


def test_abelian_label():
    assert F.abelian_label([]) == "trivial"
    assert F.abelian_label([3]) == "Z_3"
    assert F.abelian_label([2, 2]) == "Z_2 x Z_2"


def test_api_exposes_fived_block():
    out = api.summarize_toric_web([(1, 0), (0, 1), (-1, -1)])           # dP0
    fived = out["fived"]
    assert fived["rank"] == 1                                           # 1 interior point
    assert fived["flavor_rank"] == 0                                    # E_0: no flavor (B=3)
    assert fived["one_form_factors"] == [3]
    assert fived["one_form_label"] == "Z_3"
    assert fived["note"] == ""                                          # dP0 is isolated


def test_non_isolated_note():
    # C^2/Z_4 x C: bottom edge has lattice length 4 -> non-isolated, note set
    out = api.summarize_toric_web([(0, 0), (4, 0), (0, 1)])
    assert "non-isolated" in out["fived"]["note"]
    # isolated geometries stay silent
    assert F.non_isolated_note([1, 1, 1, 1]) == ""
    assert F.non_isolated_note([1, 4, 1]) != ""


def test_database_carries_fived_columns(tmp_path):
    import sqlite3
    from conformalmanifold.database import (build_toric_database,
                                            build_toric_diagram_database)
    from conformalmanifold.toric import (default_toric_library,
                                         default_toric_diagram_library)
    db = str(tmp_path / "quivers.db")
    build_toric_database(default_toric_library(), db)
    build_toric_diagram_database(default_toric_diagram_library(), db)
    conn = sqlite3.connect(db)
    try:
        # explicit-quiver table: dP0 = local P^2 (E_0)
        row = conn.execute("SELECT rank_5d, flavor_rank_5d, one_form_5d "
                           "FROM toric_quivers WHERE label='dP0'").fetchone()
        assert row == (1, 0, "Z_3")
        # diagram-only table: dP3 (E_3)
        row = conn.execute("SELECT rank_5d, flavor_rank_5d, one_form_5d "
                           "FROM toric_diagrams WHERE label='dP3'").fetchone()
        assert row == (1, 3, "trivial")
        # the research query the columns exist for: rank-1, Z_3 1-form
        hits = conn.execute("SELECT label FROM toric_quivers WHERE rank_5d=1 "
                            "AND one_form_5d='Z_3'").fetchall()
        assert ("dP0",) in hits
    finally:
        conn.close()
    # migration path: rebuilding into the same file must not error (idempotent)
    build_toric_database(default_toric_library(), db)
