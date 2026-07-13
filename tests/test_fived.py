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


# ---------------------------------------------------------------------------
# defect group D = Gamma (+) Gamma, pairing, global forms
# ---------------------------------------------------------------------------
def test_smith_normal_form_properties():
    import random
    rng = random.Random(3)
    for _ in range(25):
        m, n = rng.randint(1, 4), rng.randint(1, 4)
        A = [[rng.randint(-6, 6) for _ in range(n)] for _ in range(m)]
        P, D, Q = F.smith_normal_form(A)
        # P @ A @ Q == D
        PA = [[sum(P[i][k] * A[k][j] for k in range(m)) for j in range(n)]
              for i in range(m)]
        PAQ = [[sum(PA[i][k] * Q[k][j] for k in range(n)) for j in range(n)]
               for i in range(m)]
        assert PAQ == D
        # diagonal, non-negative, divisibility chain
        ds = []
        for i in range(m):
            for j in range(n):
                if i != j:
                    assert D[i][j] == 0
                elif D[i][j]:
                    ds.append(D[i][j])
        assert all(d > 0 for d in ds)
        assert all(ds[k + 1] % ds[k] == 0 for k in range(len(ds) - 1))


def _lens_chain_gcd(hull):
    """ADEH eqs. (3.24)-(3.25): third independent route (isolated only)."""
    from math import gcd
    from functools import reduce
    v = len(hull)
    ns = [abs((hull[i][0] - hull[i - 1][0]) * (hull[(i + 1) % v][1] - hull[i - 1][1])
              - (hull[i][1] - hull[i - 1][1]) * (hull[(i + 1) % v][0] - hull[i - 1][0]))
          for i in range(v)]
    g = reduce(gcd, ns, 0)
    return [] if g <= 1 else [g]


def test_defect_group_anchors():
    # (points, factors, pairing, #global forms) -- anchors from ADEH 2005.12831
    cases = [
        ([(0, 0), (1, 0), (0, 1)], [], [], 1),                     # C^3
        ([(0, 0), (1, 0), (1, 1), (0, 1)], [], [], 1),             # conifold
        ([(1, 0), (0, 1), (-1, -1)], [3], ["1/3"], 4),             # dP0 = E0 (eq 6.10)
        ([(1, 0), (0, 1), (-1, -1), (0, -1)], [], [], 1),          # dP1
        ([(-1, 0), (1, 0), (0, 1), (0, -1)], [2], ["1/2"], 3),     # F0 = SU(2)_0
        ([(-1, 0), (0, 0), (1, 2), (0, 4)], [2], ["1/2"], 3),      # Y(4,2) (eq 4.4)
        ([(-1, 0), (0, 0), (1, 3), (0, 6)], [3], ["1/3"], 4),      # Y(6,3)
        ([(-1, 0), (0, 0), (1, 6), (0, 9)], [3], ["1/3"], 4),      # Y(9,3)
        ([(-1, 0), (0, 0), (1, 4), (0, 4)], [4], ["1/4"], 7),      # Y(4,0)
    ]
    for pts, factors, pairing, nforms in cases:
        hull = T.convex_hull(pts)
        dg = F.defect_group(hull)
        assert dg["isolated"] is True
        assert dg["factors"] == factors
        assert dg["pairing"] == pairing
        assert dg["num_global_forms"] == nforms
        assert dg["note"] == ""
        # three independent routes agree
        assert sorted(dg["factors"]) == sorted(F.one_form_symmetry(hull))
        assert sorted(dg["factors"]) == sorted(_lens_chain_gcd(hull))


def test_defect_group_random_isolated_matches_link():
    import random
    rng = random.Random(19)
    tested = 0
    while tested < 40:
        pts = [(rng.randint(-4, 4), rng.randint(-4, 4))
               for _ in range(rng.randint(3, 8))]
        hull = T.convex_hull(pts)
        if len(hull) < 3:
            continue
        _, _, _, edges = T.polygon_signature(hull)
        if any(e != 1 for e in edges):
            continue
        tested += 1
        dg = F.defect_group(hull)
        assert sorted(dg["factors"]) == sorted(F.one_form_symmetry(hull))
        assert sorted(dg["factors"]) == sorted(_lens_chain_gcd(hull))
        assert dg["note"] == ""


def test_defect_group_non_isolated_is_caveated():
    for pts, link in [([(0, 0), (4, 0), (0, 1)], [4]),        # C^2/Z_4 x C
                      ([(0, 0), (2, 0), (0, 2)], [2, 2])]:    # C^3/(Z_2 x Z_2)
        dg = F.defect_group(T.convex_hull(pts))
        assert dg["isolated"] is False
        assert dg["factors"] == link                          # link reading kept
        assert dg["pairing"] is None
        assert dg["num_global_forms"] is None
        assert "screening" in dg["note"]


def test_cubic_anomaly_anchors():
    # published closed forms (arXiv:2112.02092: eqs. 5.33, 5.42, 4.34 = 5.48)
    def cubic(pts):
        c = F.cubic_anomaly(T.convex_hull(pts))
        return (c["status"], c["coefficient"])
    # B_N family: B_3 = dP0 = E0 -> 1/9 (in both the library frame and table-4 frame)
    assert cubic([(1, 0), (0, 1), (-1, -1)]) == ("computed", "1/9")
    assert cubic([(2, 0), (1, 2), (0, 1)]) == ("computed", "1/9")
    assert cubic([(3, 0), (1, 3), (0, 1)]) == ("computed", "1/7")      # B_4
    assert cubic([(4, 0), (1, 4), (0, 1)]) == ("computed", "2/13")     # B_5
    # su(p)_q family: q p (p-1)(p-2) / (6 gcd^3); su(2)_0 and su(p)_0 vanish
    assert cubic([(-1, 0), (1, 0), (0, 1), (0, -1)]) == ("computed", "0")    # su(2)_0 = F0
    assert cubic([(-1, 0), (0, 0), (1, 4), (0, 4)]) == ("computed", "0")     # su(4)_0
    assert cubic([(-1, 0), (0, 0), (1, 2), (0, 4)]) == ("computed", "1")     # su(4)_2
    assert cubic([(-1, 0), (0, 0), (1, 3), (0, 6)]) == ("computed", "20/9")  # su(6)_3
    assert cubic([(-1, 0), (0, 0), (1, 4), (0, 6)]) == ("computed", "5")     # su(6)_2
    # trivial 1-form symmetry -> 0
    assert cubic([(0, 0), (1, 0), (1, 1), (0, 1)]) == ("computed", "0")      # conifold
    # honest fallbacks: Jeff's exact wording
    for pts in ([(0, 0), (2, 0), (0, 2)],                 # non-isolated
                [(-4, -4), (3, 0), (-3, -1)]):            # isolated Z_17, no family
        c = F.cubic_anomaly(T.convex_hull(pts))
        assert c["status"] == "not_computed"
        assert c["note"] == "not computed — no published closed form"
        assert c["coefficient"] is None


def test_sigma_divisors():
    assert [F.sigma_divisors(d) for d in (1, 2, 3, 4, 6, 12)] == [1, 3, 4, 7, 12, 28]


def test_api_exposes_fived_block():
    out = api.summarize_toric_web([(1, 0), (0, 1), (-1, -1)])           # dP0
    fived = out["fived"]
    assert fived["rank"] == 1                                           # 1 interior point
    assert fived["flavor_rank"] == 0                                    # E_0: no flavor (B=3)
    assert fived["one_form_factors"] == [3]
    assert fived["one_form_label"] == "Z_3"
    assert fived["note"] == ""                                          # dP0 is isolated
    dg = fived["defect_group"]
    assert dg["label"] == "Z_3 (+) Z_3"
    assert dg["pairing"] == ["1/3"]
    assert dg["num_global_forms"] == 4


def test_gravity_dual_poles_and_charge_conservation():
    """The D'Hoker-Gutperle-Uhlemann dual: L external 5-brane poles = polygon
    edges; each carries (p,q) = edge outward normal x lattice length; the web
    charges must sum to zero (the polygon closes)."""
    for pts, L in [([(1, 0), (0, 1), (-1, -1)], 3),                    # dP0
                   ([(-1, 0), (1, 0), (0, 1), (0, -1)], 4),            # F0
                   ([(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)], 6)]:  # dP3
        g = F.gravity_dual(T.convex_hull(pts))
        assert g["num_external_stacks"] == L
        assert len(g["five_brane_poles"]) == L
        assert g["charge_conservation"]["conserved"] is True
        assert g["charge_conservation"]["sum_pq"] == [0, 0]
        assert g["riemann_surface"]["genus"] == 0


def test_gravity_dual_free_energy_matches_TN():
    """F_{S^5} = -(9/8) zeta(3)/pi^2 * sum_{l<k}(p_l q_k - p_k q_l)^2, exact for
    a 3-leg web.  T_N (triangle, N branes per leg) must reproduce the published
    F = -27 zeta(3)/(8 pi^2) N^4 (Fluder-Uhlemann arXiv:1806.08374)."""
    import math
    z3 = 1.2020569031595942
    for N in (1, 2, 3):
        # T_N toric triangle (N,0),(0,N),(0,0): 3 legs, N branes each
        g = F.gravity_dual(T.convex_hull([(0, 0), (N, 0), (0, N)]))
        fe = g["free_energy"]
        assert fe["exact"] is True                    # 3-leg -> no tri-log term
        assert fe["wedge_sum_S"] == 3 * N ** 4
        expected = -27 * z3 / (8 * math.pi ** 2) * N ** 4
        assert abs(fe["value"] - round(expected, 6)) < 1e-6


def test_summarize_ads6_entry_point():
    d = api.summarize_ads6([(1, 0), (0, 1), (-1, -1)])                 # dP0 = E0
    assert d["available"] and d["rank"] == 1 and d["flavor_rank"] == 0
    assert d["one_form_label"] == "Z_3"
    g = d["gravity_dual"]
    assert g["family"].startswith("D'Hoker")
    assert g["num_external_stacks"] == 3
    assert g["free_energy"]["exact"] is True
    # a degenerate 1d input is refused
    assert api.summarize_ads6([(0, 0), (1, 0)])["available"] is False


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
        row = conn.execute("SELECT rank_5d, flavor_rank_5d, one_form_5d, "
                           "defect_group_5d, pairing_5d, n_global_forms_5d, "
                           "cubic_anomaly_5d "
                           "FROM toric_quivers WHERE label='dP0'").fetchone()
        assert row == (1, 0, "Z_3", "Z_3 (+) Z_3", "1/3", 4, "1/9")
        # diagram-only table: dP3 (E_3)
        row = conn.execute("SELECT rank_5d, flavor_rank_5d, one_form_5d, "
                           "defect_group_5d, pairing_5d, n_global_forms_5d "
                           "FROM toric_diagrams WHERE label='dP3'").fetchone()
        assert row == (1, 3, "trivial", "trivial", None, 1)
        # non-isolated rows carry the link reading but no pairing/global forms
        row = conn.execute("SELECT defect_group_5d, pairing_5d, n_global_forms_5d "
                           "FROM toric_diagrams WHERE label='Z(2,2)'").fetchone()
        assert row == ("Z_2 x Z_2 (+) Z_2 x Z_2", None, None)
        # the research query the columns exist for: rank-1, Z_3 1-form
        hits = conn.execute("SELECT label FROM toric_quivers WHERE rank_5d=1 "
                            "AND one_form_5d='Z_3'").fetchall()
        assert ("dP0",) in hits
    finally:
        conn.close()
    # migration path: rebuilding into the same file must not error (idempotent)
    build_toric_database(default_toric_library(), db)
