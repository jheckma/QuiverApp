"""The inverse algorithm: a toric diagram -> a quiver gauge theory.

The *forward* direction (quiver -> toric diagram) reads the Newton polygon off the
Kasteleyn determinant of a brane tiling.  This module does the **inverse**: given
only a toric diagram (a convex lattice polygon), it reconstructs a consistent
**brane tiling** (dimer model on T^2) and hence the quiver gauge theory --

  * gauge nodes  (the faces of the tiling, #= 2 * area of the diagram),
  * oriented bifundamental fields  (the edges of the tiling),
  * a *toric* (two-term-per-field) superpotential  (the vertices of the tiling,
    white = +, black = -),

for **any** convex lattice polygon, following Gulotta's "properly ordered dimer"
construction (arXiv:0807.3012).

Method (zig-zag / "alga" picture, Hanany-Kennaway, Stienstra).  A brane tiling is
recovered from its **zig-zag paths**: closed curves on T^2 whose homology winding
numbers are exactly the primitive outward normals of the toric-diagram edges --
i.e. the (p,q) legs of the dual 5-brane web (one zig-zag per primitive boundary
segment, so B of them for B boundary lattice points).  Drawn as straight geodesics
on T^2 the zig-zag curves form an arrangement that is precisely the **medial
graph** of the tiling:

  * each crossing of two zig-zag curves  <->  an edge of the tiling (a field);
    paths with windings w_k, w_l cross |det(w_k, w_l)| times, so the number of
    fields is  sum_{k<l} |det(w_k, w_l)|;
  * each face of the arrangement whose bounding arcs all run *forward* along their
    paths is a **white** tiling vertex; all *backward*, a **black** vertex; faces
    with alternating (forward/backward) arcs are the **gauge** faces.

"Properly ordered" = the curves are placed so the arrangement is the minimal,
non-degenerate one (no triple points; every crossing 4-valent; faces correctly
typed).  We realise it by searching base-point offsets until the extracted tiling
is consistent (correct counts, anomaly-free, genuinely two-term superpotential);
for the small polygons of interest a valid placement is found in a few tries.

Because the zig-zag windings are taken *directly* from the input polygon's edge
normals, a consistent tiling produced this way necessarily has the input toric
diagram as its Newton polygon -- the construction is self-certifying on the
geometry; the consistency checks certify that it is a valid dimer.

Pure standard library (Fraction-exact crossings; no numpy needed here).
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction as Fr
from math import gcd

from .toric import convex_hull, gl2z_equiv, normalized_area


# ===========================================================================
# zig-zag windings = primitive outward normals of the toric-diagram edges
# ===========================================================================
def zigzag_windings(hull):
    """One primitive (p,q) winding per primitive boundary segment of `hull`
    (CCW), equal to the outward edge normal -- the (p,q)-web legs."""
    ws = []
    n = len(hull)
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        g = gcd(abs(dx), abs(dy)) or 1
        ux, uy = dx // g, dy // g
        for _ in range(g):
            ws.append((uy, -ux))            # outward normal of a CCW polygon
    return ws


def _det(a, b):
    return a[0] * b[1] - a[1] * b[0]


# ===========================================================================
# torus arrangement of the zig-zag geodesics
# ===========================================================================
def _crossings(ws, bases):
    """All crossings of the zig-zag geodesics on T^2 = R^2/Z^2.

    Path k is { bases[k] + s*ws[k] (mod Z^2) }.  Returns a list of
    (k, l, s_k, s_l) with s in [0,1) Fractions; exactly |det(w_k,w_l)| per pair.
    """
    B = len(ws)
    out = []
    for k in range(B):
        for l in range(k + 1, B):
            D = _det(ws[k], ws[l])
            if D == 0:
                continue                    # parallel zig-zags never cross
            need = abs(D)
            pk, qk = ws[k]
            pl, ql = ws[l]
            detM = Fr(-D)                    # det[[pk,-pl],[qk,-ql]]
            found = {}
            rng = 4 + max(abs(pk), abs(qk), abs(pl), abs(ql))
            while len(found) < need:
                for m in range(-rng, rng + 1):
                    for nn in range(-rng, rng + 1):
                        rx = Fr(bases[l][0] - bases[k][0] + m)
                        ry = Fr(bases[l][1] - bases[k][1] + nn)
                        s = ((rx * (-ql) - (-pl) * ry) / detM) % 1
                        t = ((pk * ry - qk * rx) / detM) % 1
                        found[(s, t)] = (k, l, s, t)
                rng += 3
                if rng > 60:                 # safety; should never trigger
                    break
            out.extend(found.values())
    return out


@dataclass
class BraneTiling:
    """A consistent brane tiling (dimer) reconstructed from a toric diagram."""
    num_gauge: int
    num_fields: int
    num_white: int
    num_black: int
    adjacency: list                 # num_gauge x num_gauge, A[i][j] = #arrows i->j
    fields: list                    # [{"label","src","tgt"}], one per tiling edge
    superpotential: list            # [{"sign":+1/-1,"fields":[labels...]}]
    white_pos: list                 # [(x,y)] in [0,1)^2 (torus)
    black_pos: list
    field_pos: list                 # [(x,y)] one per field (the crossing site)
    tiling_edges: list              # [[white_idx, black_idx]] one per field
    checks: dict
    note: str = ""

    def adjacency_int(self):
        return [[int(v) for v in row] for row in self.adjacency]


def _try_tiling(hull, ws, bases):
    """One attempt at extracting the tiling for a given base-point placement;
    returns a BraneTiling or None if the placement is degenerate/inconsistent."""
    B = len(ws)
    a2 = normalized_area(hull)
    cr = _crossings(ws, bases)
    E = len(cr)

    # order each path's crossings by parameter; reject coincident params
    along = {k: [] for k in range(B)}
    for ci, (k, l, sk, sl) in enumerate(cr):
        along[k].append((sk, ci, "k"))
        along[l].append((sl, ci, "l"))
    for k in along:
        params = [x[0] for x in along[k]]
        if len(set(params)) != len(params):
            return None                      # triple point / coincidence
        along[k].sort()

    # combinatorial map: alpha (arc involution) and sigma (rotation at a crossing)
    alpha = {}
    for p in range(B):
        seq = along[p]
        L = len(seq)
        for idx in range(L):
            _, ci, side = seq[idx]
            _, ci2, side2 = seq[(idx + 1) % L]
            alpha[(ci, side, 1)] = (ci2, side2, -1)
            alpha[(ci2, side2, -1)] = (ci, side, 1)
    sigma = {}
    for ci, (k, l, sk, sl) in enumerate(cr):
        dk, dl = ws[k], ws[l]
        darts = [((ci, "k", 1), dk), ((ci, "k", -1), (-dk[0], -dk[1])),
                 ((ci, "l", 1), dl), ((ci, "l", -1), (-dl[0], -dl[1]))]
        darts.sort(key=lambda x: math.atan2(x[1][1], x[1][0]))
        order = [d[0] for d in darts]
        for i in range(4):
            sigma[order[i]] = order[(i + 1) % 4]

    def phi(d):
        return sigma[alpha[d]]

    # faces = orbits of phi
    seen = set()
    faces = []
    faceof = {}
    for d in alpha:
        if d in seen:
            continue
        orb = []
        x = d
        while x not in seen:
            seen.add(x)
            orb.append(x)
            x = phi(x)
        fid = len(faces)
        faces.append(orb)
        for x in orb:
            faceof[x] = fid

    # classify faces by the path-direction signs of their bounding darts
    kind = {}
    for i, f in enumerate(faces):
        signs = set(d[2] for d in f)
        kind[i] = "W" if signs == {1} else "B" if signs == {-1} else "G"
    gauge = [f for f in range(len(faces)) if kind[f] == "G"]
    nW = sum(1 for v in kind.values() if v == "W")
    nB = sum(1 for v in kind.values() if v == "B")
    if not (len(gauge) == a2 and nW == nB and nW + nB == E - a2):
        return None

    # quiver adjacency: each crossing borders 2 gauge faces; orient white->black
    gindex = {g: i for i, g in enumerate(gauge)}
    A = [[0] * len(gauge) for _ in range(len(gauge))]
    field_edges = {}            # ci -> (src gauge idx, tgt gauge idx)
    for ci in range(E):
        start = (ci, "k", 1)
        ring = [start]
        x = sigma[start]
        while x != start:
            ring.append(x)
            x = sigma[x]
            if len(ring) > 4:
                return None
        if len(ring) != 4:
            return None
        fring = [faceof[d] for d in ring]
        gi = [i for i, f in enumerate(fring) if kind[f] == "G"]
        if len(gi) != 2:
            return None
        i0 = gi[0]
        g1, nb, g2 = fring[i0], fring[(i0 + 1) % 4], fring[(i0 + 2) % 4]
        if kind[nb] == "W":
            src, tgt = gindex[g1], gindex[g2]
        else:
            src, tgt = gindex[g2], gindex[g1]
        A[src][tgt] += 1
        field_edges[ci] = (src, tgt)

    if any(sum(A[i]) != sum(A[j][i] for j in range(len(gauge)))
           for i in range(len(gauge))):
        return None                          # gauge anomaly

    # superpotential: white(+)/black(-) faces -> the fields (crossings) around them
    cp, cm = Counter(), Counter()
    W = []
    wfaces, bfaces = [], []
    for f in range(len(faces)):
        if kind[f] == "G":
            continue
        flds = [d[0] for d in faces[f]]
        sign = 1 if kind[f] == "W" else -1
        W.append({"sign": sign, "fields": [f"X{c}" for c in flds]})
        (wfaces if sign > 0 else bfaces).append(f)
        for c in flds:
            (cp if sign > 0 else cm)[c] += 1
    if not all(cp[c] == 1 and cm[c] == 1 for c in range(E)):
        return None                          # not a genuine two-term W

    # geometry for drawing: crossing site, and node = mean of its crossings
    def site(ci):
        k, l, sk, sl = cr[ci]
        x = (bases[k][0] + sk * ws[k][0]) % 1
        y = (bases[k][1] + sk * ws[k][1]) % 1
        return (float(x), float(y))

    field_pos = [site(ci) for ci in range(E)]

    def node_pos(f):
        cs = [d[0] for d in faces[f]]
        # average on the torus via mean of unit-circle angles per axis
        import cmath
        zx = sum(cmath.exp(2j * math.pi * field_pos[c][0]) for c in cs)
        zy = sum(cmath.exp(2j * math.pi * field_pos[c][1]) for c in cs)
        ax = (math.atan2(zx.imag, zx.real) / (2 * math.pi)) % 1
        ay = (math.atan2(zy.imag, zy.real) / (2 * math.pi)) % 1
        return (ax, ay)

    white_pos = [node_pos(f) for f in wfaces]
    black_pos = [node_pos(f) for f in bfaces]

    # each field (crossing) lies on exactly one white and one black tiling vertex
    w_of, b_of = {}, {}
    for wi, f in enumerate(wfaces):
        for d in faces[f]:
            w_of[d[0]] = wi
    for bi, f in enumerate(bfaces):
        for d in faces[f]:
            b_of[d[0]] = bi
    tiling_edges = [[w_of[ci], b_of[ci]] for ci in range(E)]

    # ----- exact edge-homology cochain (for the Kasteleyn Newton-polygon) -----
    # Each tiling edge ci gets a Z^2 class h_ci so that  det K(x,y), with entry
    # K[w][b] = sum_e x^{h_e.x} y^{h_e.y}, has Newton polygon = the toric diagram
    # (up to GL(2,Z)+translation).  Built exactly (Fraction arithmetic) by:
    #   (1) walking each white/black face's boundary, accumulating the arc
    #       displacements (param-step * winding) -> corners in the face's frame;
    #   (2) a spanning tree over the (white,black) graph that re-lifts each face
    #       so a shared crossing coincides -> one global frame;
    #   (3) h_ci = (black corner) - (white corner) of crossing ci in that frame
    #       (tree edges -> 0, cotree edges carry the winding).
    disp_plus, disp_minus = {}, {}
    for p in range(B):
        seq = along[p]
        L = len(seq)
        for idx in range(L):
            param_i, ci_i, side_i = seq[idx]
            param_j, ci_j, side_j = seq[(idx + 1) % L]
            ds = (param_j - param_i) % 1                 # Fraction in (0,1]
            v = (ds * ws[p][0], ds * ws[p][1])
            disp_plus[(ci_i, side_i)] = v
            disp_minus[(ci_j, side_j)] = (-v[0], -v[1])

    def _alpha_disp(d):
        return disp_plus[(d[0], d[1])] if d[2] == 1 else disp_minus[(d[0], d[1])]

    posd_face = {}                  # face id -> {dart: (Fr x, Fr y)} corner pos
    w_dart, b_dart = {}, {}
    for f in wfaces + bfaces:
        orb = faces[f]
        pos = {orb[0]: (Fr(0), Fr(0))}
        cur = (Fr(0), Fr(0))
        for i in range(len(orb)):
            dd = _alpha_disp(orb[i])
            cur = (cur[0] + dd[0], cur[1] + dd[1])
            j = (i + 1) % len(orb)
            if j != 0:
                pos[orb[j]] = cur
        posd_face[f] = pos
    for wi, f in enumerate(wfaces):
        for d in faces[f]:
            w_dart[d[0]] = d
    for bi, f in enumerate(bfaces):
        for d in faces[f]:
            b_dart[d[0]] = d

    def _local_white(ci):
        return posd_face[wfaces[w_of[ci]]][w_dart[ci]]

    def _local_black(ci):
        return posd_face[bfaces[b_of[ci]]][b_dart[ci]]

    # spanning tree over white/black faces, fields as edges -> global offsets
    fadj = {}
    for ci in range(E):
        u, v = ("W", w_of[ci]), ("B", b_of[ci])
        fadj.setdefault(u, []).append((v, ci))
        fadj.setdefault(v, []).append((u, ci))
    offset = {("W", 0): (Fr(0), Fr(0))}

    def _gpos(node, ci):
        loc = _local_white(ci) if node[0] == "W" else _local_black(ci)
        off = offset[node]
        return (loc[0] + off[0], loc[1] + off[1])

    stack = [("W", 0)]
    while stack:
        u = stack.pop()
        for (v, ci) in fadj[u]:
            if v in offset:
                continue
            gu = _gpos(u, ci)
            loc_v = _local_white(ci) if v[0] == "W" else _local_black(ci)
            offset[v] = (gu[0] - loc_v[0], gu[1] - loc_v[1])
            stack.append(v)

    homology = {}
    for ci in range(E):
        cw = _gpos(("W", w_of[ci]), ci)
        cb = _gpos(("B", b_of[ci]), ci)
        hx, hy = cb[0] - cw[0], cb[1] - cw[1]
        if hx.denominator != 1 or hy.denominator != 1:
            return None                          # inconsistent lift -> reject
        homology[ci] = (int(hx), int(hy))

    fields = [{"label": f"X{ci}", "src": field_edges[ci][0],
               "tgt": field_edges[ci][1],
               "white": w_of[ci], "black": b_of[ci],
               "zigzag": [cr[ci][0], cr[ci][1]],
               "homology": list(homology[ci])} for ci in range(E)]

    checks = {
        "gauge_eq_2area": len(gauge) == a2,
        "fields_eq_sum_det": E == sum(abs(_det(ws[k], ws[l]))
                                      for k in range(B) for l in range(k + 1, B)),
        "white_eq_black": nW == nB,
        "anomaly_free": True,
        "toric_superpotential": True,
        "euler_V_minus_E_plus_F": (nW + nB) - E + len(gauge),   # = 0 on T^2
    }
    return BraneTiling(
        num_gauge=len(gauge), num_fields=E, num_white=nW, num_black=nB,
        adjacency=A, fields=fields, superpotential=W,
        white_pos=white_pos, black_pos=black_pos, field_pos=field_pos,
        tiling_edges=tiling_edges, checks=checks,
    )


def inverse_quiver(vertices, max_attempts: int = 400, max_gauge: int = 60):
    """Reconstruct a quiver gauge theory + brane tiling from a toric diagram.

    `vertices` : lattice points / corners of the toric diagram (any order;
                 the convex hull is taken).
    Returns a `BraneTiling`, or raises `ValueError` if the diagram is degenerate
    or no consistent placement is found within `max_attempts`.
    """
    hull = convex_hull(vertices)
    if len(hull) < 3:
        raise ValueError("need at least 3 non-collinear lattice points")
    a2 = normalized_area(hull)
    if a2 > max_gauge:
        raise ValueError(f"toric diagram too large (2*area = {a2} gauge nodes > "
                         f"max_gauge = {max_gauge})")
    ws = zigzag_windings(hull)

    # deterministic pseudo-random base-point search (reproducible across runs)
    state = 0x9E3779B1
    P = 1_000_003

    def rnd():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state

    for attempt in range(max_attempts):
        bases = [(Fr(rnd() % P, P), Fr(rnd() % P, P)) for _ in range(len(ws))]
        res = _try_tiling(hull, ws, bases)
        if res is not None:
            res.note = (f"Gulotta inverse algorithm; properly-ordered placement "
                        f"found on attempt {attempt + 1}.")
            return res
    raise ValueError("no consistent brane tiling found; try again or report the "
                     "diagram (the placement search did not converge)")


def kasteleyn_newton_polygon(tiling: BraneTiling):
    """Newton polygon of the Kasteleyn determinant of `tiling`.

    Builds the (n_white x n_black) Kasteleyn matrix  K[w][b] = sum over edges
    w->b of  x^{h.x} y^{h.y}  (h = the exact edge-homology cochain stored on
    each field), forms det K as a Laurent polynomial, and returns the convex
    hull of its monomial exponents.

    For a consistent brane tiling this Newton polygon equals the input toric
    diagram up to GL(2,Z) + translation -- an *independent* certificate of the
    inverse construction (the extremal monomials come from unique perfect
    matchings, so they never cancel regardless of sign convention).
    """
    from itertools import permutations

    n = tiling.num_white
    K = [[dict() for _ in range(n)] for _ in range(n)]
    for f in tiling.fields:
        w, b = f["white"], f["black"]
        key = (f["homology"][0], f["homology"][1])
        K[w][b][key] = K[w][b].get(key, 0) + 1

    def _pmul(a, b):
        out = {}
        for (i, j), u in a.items():
            for (k, l), v in b.items():
                p = (i + k, j + l)
                out[p] = out.get(p, 0) + u * v
        return out

    det = {}
    for perm in permutations(range(n)):
        sign, seen = 1, [False] * n
        for i in range(n):
            if seen[i]:
                continue
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j, ln = perm[j], ln + 1
            if ln % 2 == 0:
                sign = -sign
        term, ok = {(0, 0): sign}, True
        for i in range(n):
            e = K[i][perm[i]]
            if not e:
                ok = False
                break
            term = _pmul(term, e)
        if not ok:
            continue
        for k, v in term.items():
            det[k] = det.get(k, 0) + v
    exps = [k for k, v in det.items() if v != 0]
    return convex_hull(exps)


def inverse_quiver_json(vertices, **kw) -> dict:
    """JSON-friendly inverse-algorithm result for the web API (or an error)."""
    try:
        t = inverse_quiver(vertices, **kw)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    # independent certificate: Newton polygon of det K == input toric diagram
    newton = kasteleyn_newton_polygon(t)
    checks = dict(t.checks)
    checks["kasteleyn_newton_matches"] = gl2z_equiv(newton, convex_hull(vertices))
    return {
        "available": True,
        "num_gauge": t.num_gauge,
        "num_fields": t.num_fields,
        "num_white": t.num_white,
        "num_black": t.num_black,
        "adjacency": t.adjacency_int(),
        "fields": t.fields,
        "superpotential": t.superpotential,
        "kasteleyn_newton": [list(v) for v in newton],
        "tiling": {
            "white": [[round(x, 5), round(y, 5)] for (x, y) in t.white_pos],
            "black": [[round(x, 5), round(y, 5)] for (x, y) in t.black_pos],
            "fields": [[round(x, 5), round(y, 5)] for (x, y) in t.field_pos],
            "edges": t.tiling_edges,
        },
        "checks": checks,
        "note": t.note,
    }
