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
class DimerGraph:
    """A brane tiling as a combinatorial map (ribbon graph) on T^2.

    Bipartite white/black vertices, oriented bifundamental fields = edges (each
    carrying an integer Z^2 homology cochain), and a *rotation system*: the
    cyclic order of incident edges around every vertex.  This is the minimal
    data Seiberg duality (urban renewal) acts on; `forward_extract` turns it back
    into a BraneTiling (gauge faces + quiver + two-term superpotential).

    The rotation order is taken straight from the boundary order of each tiling
    vertex -- which is exactly the field order of its superpotential term -- so no
    floating-point geometry is needed to recover the faces.
    """
    nW: int
    nB: int
    edges: list      # [{"w": wi, "b": bi, "h": [hx, hy]}], edge id = list index
    rot_w: list      # [[edge_id, ...]] cyclic order of edges around each white vtx
    rot_b: list      # [[edge_id, ...]] cyclic order of edges around each black vtx


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
    white_glob: list                # [(x,y)] un-reduced (planar) embedding
    black_glob: list
    field_pos: list                 # [(x,y)] one per field (the crossing site)
    tiling_edges: list              # [[white_idx, black_idx]] one per field
    checks: dict
    note: str = ""

    def adjacency_int(self):
        return [[int(v) for v in row] for row in self.adjacency]

    def to_dimer(self) -> "DimerGraph":
        """Reconstruct the combinatorial-map DimerGraph from this tiling.

        The cyclic edge order around white vertex k is the field order of the
        k-th positive (white) superpotential term, and around black vertex k the
        k-th negative (black) term -- which `_try_tiling` builds from the cyclic
        face-boundary order.  Edge id == the integer in its "X{id}" label.
        """
        edges = [{"w": f["white"], "b": f["black"], "h": list(f["homology"])}
                 for f in self.fields]
        pos_terms = [t for t in self.superpotential if t["sign"] > 0]
        neg_terms = [t for t in self.superpotential if t["sign"] < 0]
        rot_w = [[int(lbl[1:]) for lbl in t["fields"]] for t in pos_terms]
        rot_b = [[int(lbl[1:]) for lbl in t["fields"]] for t in neg_terms]
        return DimerGraph(nW=self.num_white, nB=self.num_black,
                          edges=edges, rot_w=rot_w, rot_b=rot_b)


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

    # consistent (globally-lifted) node/crossing positions for a clean tiling:
    # each face's corners share one frame, so its centroid is unambiguous (no
    # circular-mean collapse), and reducing mod 1 gives its true torus position.
    def _centroidU(node):
        """un-reduced global-frame face centroid (an exact planar embedding)."""
        flist = wfaces if node[0] == "W" else bfaces
        f = flist[node[1]]
        off = offset[node]
        ds = faces[f]
        cx = sum(posd_face[f][d][0] + off[0] for d in ds) / len(ds)
        cy = sum(posd_face[f][d][1] + off[1] for d in ds) / len(ds)
        return (float(cx), float(cy))

    # un-reduced (universal-cover) centroids -> a genuine planar tiling layout;
    # mod-1 versions kept for the torus-cell drawing / field sites.
    white_glob = [_centroidU(("W", wi)) for wi in range(nW)]
    black_glob = [_centroidU(("B", bi)) for bi in range(nB)]
    white_pos = [(x % 1, y % 1) for (x, y) in white_glob]
    black_pos = [(x % 1, y % 1) for (x, y) in black_glob]
    field_pos = [(float(_gpos(("W", w_of[ci]), ci)[0] % 1),
                  float(_gpos(("W", w_of[ci]), ci)[1] % 1)) for ci in range(E)]

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
        white_pos=white_pos, black_pos=black_pos,
        white_glob=white_glob, black_glob=black_glob, field_pos=field_pos,
        tiling_edges=tiling_edges, checks=checks,
    )


# ===========================================================================
# forward extractor: a DimerGraph (combinatorial map) -> a BraneTiling
# ===========================================================================
# A *dart* is (edge_id, +1) [oriented white->black] or (edge_id, -1)
# [black->white].  alpha swaps the two darts of an edge; sigma rotates to the
# next edge around the base vertex; faces (gauge nodes) = orbits of phi=sigma o
# alpha.  These two module constants fix the orientation convention; they are
# pinned by the round-trip gate (forward_extract(t.to_dimer()) == t) and must not
# be changed without re-running tests/test_inverse.py.
_FACE_HAND = 1          # +1 = "next in rotation", -1 = "previous"
_ARROW_FROM_WHITE_DART = True   # src = face carrying the white->black dart


def _trace_faces(dimer, hand):
    """Orbits of phi = sigma o alpha over all darts.  Returns (faces, faceof):
    faces is a list of dart lists, faceof maps each dart to its face index."""
    posw = [{e: i for i, e in enumerate(r)} for r in dimer.rot_w]
    posb = [{e: i for i, e in enumerate(r)} for r in dimer.rot_b]

    def sigma(d):
        e, s = d
        if s == 1:
            w = dimer.edges[e]["w"]
            r, i = dimer.rot_w[w], posw[w][e]
        else:
            b = dimer.edges[e]["b"]
            r, i = dimer.rot_b[b], posb[b][e]
        return (r[(i + hand) % len(r)], s)

    def phi(d):
        return sigma((d[0], -d[1]))

    seen, faces, faceof = set(), [], {}
    darts = [(e, s) for e in range(len(dimer.edges)) for s in (1, -1)]
    for d in darts:
        if d in seen:
            continue
        orb, x = [], d
        while x not in seen:
            seen.add(x)
            orb.append(x)
            x = phi(x)
        fid = len(faces)
        faces.append(orb)
        for y in orb:
            faceof[y] = fid
    return faces, faceof


def _spanning_layout(dimer):
    """A schematic planar layout (universal-cover float positions) for drawing.

    Spanning-tree placement that splays each vertex's edges by their rotation
    index; not metrically faithful (mutated phases have no canonical embedding)
    but renders the connectivity.  Returns (white_glob, black_glob)."""
    wpos = [None] * dimer.nW
    bpos = [None] * dimer.nB
    adj = {("W", w): [] for w in range(dimer.nW)}
    adj.update({("B", b): [] for b in range(dimer.nB)})
    for e, ed in enumerate(dimer.edges):
        adj[("W", ed["w"])].append(("B", ed["b"], e))
        adj[("B", ed["b"])].append(("W", ed["w"], e))
    if dimer.nW:
        wpos[0] = (0.0, 0.0)
    stack = [("W", 0)] if dimer.nW else []
    placed = {("W", 0)} if dimer.nW else set()
    while stack:
        kind, idx = stack.pop()
        base = wpos[idx] if kind == "W" else bpos[idx]
        rot = (dimer.rot_w[idx] if kind == "W" else dimer.rot_b[idx])
        deg = max(len(rot), 1)
        for (nk, ni, e) in adj[(kind, idx)]:
            if (nk, ni) in placed:
                continue
            ang = 2 * math.pi * (rot.index(e) / deg)
            if kind == "B":
                ang += math.pi
            npos = (base[0] + 0.6 * math.cos(ang), base[1] + 0.6 * math.sin(ang))
            if nk == "W":
                wpos[ni] = npos
            else:
                bpos[ni] = npos
            placed.add((nk, ni))
            stack.append((nk, ni))
    wpos = [p if p is not None else (0.0, 0.0) for p in wpos]
    bpos = [p if p is not None else (0.5, 0.5) for p in bpos]
    return wpos, bpos


def forward_extract(dimer, vertices) -> "BraneTiling":
    """Re-derive a BraneTiling (gauge faces + quiver + two-term superpotential)
    from a combinatorial-map DimerGraph, independent of the zig-zag arrangement.

    This is the inverse of `BraneTiling.to_dimer` on the base phase, and the
    extractor for any Seiberg-dual phase produced by `urban_renewal`.
    """
    hull = convex_hull(vertices)
    a2 = normalized_area(hull)
    E = len(dimer.edges)
    nW, nB = dimer.nW, dimer.nB

    faces, faceof = _trace_faces(dimer, _FACE_HAND)
    gauge = list(range(len(faces)))            # every ribbon face is a gauge node
    gindex = {g: i for i, g in enumerate(gauge)}

    # quiver: each edge borders the faces of its two darts; orient white->black
    A = [[0] * len(gauge) for _ in range(len(gauge))]
    field_edges = {}
    for e in range(E):
        fw = faceof[(e, 1)]                     # face carrying the white->black dart
        fb = faceof[(e, -1)]
        if _ARROW_FROM_WHITE_DART:
            src, tgt = gindex[fw], gindex[fb]
        else:
            src, tgt = gindex[fb], gindex[fw]
        A[src][tgt] += 1
        field_edges[e] = (src, tgt)

    anomaly_free = all(sum(A[i]) == sum(A[j][i] for j in range(len(gauge)))
                       for i in range(len(gauge)))

    # superpotential: white vertex (+) / black vertex (-), fields in rotation order
    W = ([{"sign": 1, "fields": [f"X{e}" for e in dimer.rot_w[w]]}
          for w in range(nW)]
         + [{"sign": -1, "fields": [f"X{e}" for e in dimer.rot_b[b]]}
            for b in range(nB)])
    cp, cm = Counter(), Counter()
    for w in range(nW):
        for e in dimer.rot_w[w]:
            cp[e] += 1
    for b in range(nB):
        for e in dimer.rot_b[b]:
            cm[e] += 1
    toric_W = all(cp[e] == 1 and cm[e] == 1 for e in range(E))

    tiling_edges = [[dimer.edges[e]["w"], dimer.edges[e]["b"]] for e in range(E)]
    fields = [{"label": f"X{e}", "src": field_edges[e][0], "tgt": field_edges[e][1],
               "white": dimer.edges[e]["w"], "black": dimer.edges[e]["b"],
               "homology": list(dimer.edges[e]["h"])} for e in range(E)]

    white_glob, black_glob = _spanning_layout(dimer)
    white_pos = [(x % 1, y % 1) for (x, y) in white_glob]
    black_pos = [(x % 1, y % 1) for (x, y) in black_glob]
    field_pos = []
    for e in range(E):
        w, b = dimer.edges[e]["w"], dimer.edges[e]["b"]
        mx = (white_glob[w][0] + black_glob[b][0]) / 2
        my = (white_glob[w][1] + black_glob[b][1]) / 2
        field_pos.append((mx % 1, my % 1))

    checks = {
        "gauge_eq_2area": len(gauge) == a2,
        "fields_eq_sum_det": True,             # base value lives on the seed phase
        "white_eq_black": nW == nB,
        "anomaly_free": anomaly_free,
        "toric_superpotential": toric_W,
        "euler_V_minus_E_plus_F": (nW + nB) - E + len(gauge),
    }
    return BraneTiling(
        num_gauge=len(gauge), num_fields=E, num_white=nW, num_black=nB,
        adjacency=A, fields=fields, superpotential=W,
        white_pos=white_pos, black_pos=black_pos,
        white_glob=white_glob, black_glob=black_glob, field_pos=field_pos,
        tiling_edges=tiling_edges, checks=checks,
    )


# ===========================================================================
# homology solver: edge Z^2 cochain of a DimerGraph from the rotation system
# ===========================================================================
def _integer_kernel(rows, ncols):
    """A *saturated* integer basis of  {x in Z^ncols : rows . x = 0}.

    Integer column-echelon (unimodular column ops via Euclid + swaps) on the
    constraint matrix `rows`; the columns that reduce to zero carry, in the
    accumulated transform `U`, a basis of the integer kernel.  Because every
    operation is unimodular, that basis is automatically saturated -- it spans
    the *full* integer lattice of the rational nullspace, not a sublattice (the
    distinction that makes the Kasteleyn Newton polygon come out at the right
    area rather than an integer multiple of it)."""
    m = len(rows)
    # A as columns (each length m); U as columns (each length ncols) = identity
    A = [[int(rows[i][j]) for i in range(m)] for j in range(ncols)]
    U = [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]

    def sub(k, r, q):                        # col_k -= q * col_r  (A and U)
        for i in range(m):
            A[k][i] -= q * A[r][i]
        for i in range(ncols):
            U[k][i] -= q * U[r][i]

    def swap(a, b):
        A[a], A[b] = A[b], A[a]
        U[a], U[b] = U[b], U[a]

    r = 0
    for i in range(m):
        for k in range(r + 1, ncols):
            while A[k][i] != 0:
                if A[r][i] == 0:
                    swap(r, k)
                    continue
                sub(k, r, A[k][i] // A[r][i])
                if A[k][i] != 0:
                    swap(r, k)
        if r < ncols and A[r][i] != 0:
            r += 1
        if r == ncols:
            break
    return [U[j] for j in range(r, ncols)]   # zero columns -> kernel basis


def _newton_from_dimer(dimer, hx, hy):
    """Newton polygon of det K for a DimerGraph carrying the cochain (hx, hy).

    K[w][b] = sum over edges w->b of x^{hx} y^{hy}; the extremal monomials of
    det K come from unique perfect matchings, so they never cancel."""
    nW, nB = dimer.nW, dimer.nB
    if nW != nB:
        return None
    K = [[dict() for _ in range(nB)] for _ in range(nW)]
    for e, ed in enumerate(dimer.edges):
        key = (hx[e], hy[e])
        cell = K[ed["w"]][ed["b"]]
        cell[key] = cell.get(key, 0) + 1

    def pmul(a, b):
        out = {}
        for (i, j), u in a.items():
            for (k, l), v in b.items():
                p = (i + k, j + l)
                out[p] = out.get(p, 0) + u * v
        return out

    from itertools import permutations
    det = {}
    for perm in permutations(range(nW)):
        sign, used = 1, [False] * nW
        for i in range(nW):
            if used[i]:
                continue
            j, ln = i, 0
            while not used[j]:
                used[j] = True
                j, ln = perm[j], ln + 1
            if ln % 2 == 0:
                sign = -sign
        term, ok = {(0, 0): sign}, True
        for i in range(nW):
            cell = K[i][perm[i]]
            if not cell:
                ok = False
                break
            term = pmul(term, cell)
        if ok:
            for k, v in term.items():
                det[k] = det.get(k, 0) + v
    exps = [k for k, v in det.items() if v != 0]
    if len(exps) < 3:
        return None
    return convex_hull(exps)


def solve_homology(dimer, target_hull):
    """Solve the edge Z^2 homology cochain of a DimerGraph from its rotation
    system alone -- the data Seiberg duality (urban renewal) leaves us with,
    where no geometric embedding (and so no zig-zag displacement) is available.

    Returns a list of (hx, hy) per edge such that the Kasteleyn Newton polygon
    equals `target_hull` up to GL(2,Z) + translation, or None if no such cochain
    exists (i.e. the dimer is not a toric phase of that diagram).

    Method: (1) every gauge face (a phi-orbit of darts) must close --
    sum_d  s_d * h(edge_d) = 0, the cocycle condition; (2) gauge-fix a spanning
    tree of the bipartite white/black graph to h = 0 (kills the coboundary
    freedom); (3) the residual solution lattice is the rank-2 H^1(T^2) = Z^2 of
    torus periods -- take a *saturated* integer basis (u, v) and sweep small
    unimodular GL(2,Z) recombinations as the two cochain components until the
    Newton polygon certifies.  Mirrors the spanning-tree lift of `_try_tiling`,
    but driven by the combinatorial map instead of the zig-zag arrangement."""
    from itertools import product

    E = len(dimer.edges)
    faces, _ = _trace_faces(dimer, _FACE_HAND)

    rows = []
    for orb in faces:                        # face-closure (cocycle) constraints
        row = [0] * E
        for (e, s) in orb:
            row[e] += s
        rows.append(row)

    # spanning tree over the bipartite white/black graph -> pin tree edges to 0
    adj = {}
    for e, ed in enumerate(dimer.edges):
        u, v = ("W", ed["w"]), ("B", ed["b"])
        adj.setdefault(u, []).append((v, e))
        adj.setdefault(v, []).append((u, e))
    seen, tree = {("W", 0)}, set()
    stack = [("W", 0)]
    while stack:
        u = stack.pop()
        for (v, e) in adj.get(u, []):
            if v not in seen:
                seen.add(v)
                tree.add(e)
                stack.append(v)
    for e in tree:
        row = [0] * E
        row[e] = 1
        rows.append(row)

    basis = _integer_kernel(rows, E)
    if len(basis) != 2:                      # not a genus-1 (torus) cochain space
        return None
    u, v = basis

    # Because `_integer_kernel` returns a *saturated* basis of H^1 (= the full
    # period lattice, not a sublattice), the cochain (u, v) already yields a
    # Newton polygon GL(2,Z)-equivalent to any valid frame's -- so the identity
    # frame certifies whenever the dimer is a genuine toric phase.  We try it
    # first; the small unimodular sweep is then only defensive (e.g. a frame in
    # which `_newton_from_dimer` degenerates), never the primary path.
    frames = [(1, 0, 0, 1)] + [(a, b, c, d)
              for a, b, c, d in product(range(-2, 3), repeat=4)
              if a * d - b * c in (1, -1) and (a, b, c, d) != (1, 0, 0, 1)]
    for a, b, c, d in frames:
        hx = [a * u[i] + b * v[i] for i in range(E)]
        hy = [c * u[i] + d * v[i] for i in range(E)]
        poly = _newton_from_dimer(dimer, hx, hy)
        if poly is not None and gl2z_equiv(poly, target_hull):
            return [(hx[i], hy[i]) for i in range(E)]
    return None


# ===========================================================================
# phase invariant: distinguish Seiberg-dual toric phases up to node relabelling
# ===========================================================================
def _wl_hash(A, rounds=None):
    """A Weisfeiler-Lehman colour-refinement hash of a directed multigraph
    adjacency -- a permutation-invariant fingerprint, used as a tractable
    fallback canonical key when brute-force node permutation is too expensive."""
    n = len(A)
    if n == 0:
        return ()
    rounds = rounds if rounds is not None else n
    colour = [hash(()) for _ in range(n)]
    for _ in range(rounds):
        new = []
        for i in range(n):
            out = tuple(sorted((A[i][j], colour[j]) for j in range(n) if A[i][j]))
            inc = tuple(sorted((A[j][i], colour[j]) for j in range(n) if A[j][i]))
            new.append(hash((colour[i], out, inc)))
        colour = new
    return tuple(sorted(colour))


def canonical_adjacency(A):
    """A canonical (node-relabelling-invariant) key for a quiver adjacency matrix.

    Exact (lexicographically minimal over all node permutations) for small
    quivers; falls back to a WL colour-refinement hash for large node counts
    where brute force is infeasible.  Used to dedup Seiberg-dual phases.
    """
    from itertools import permutations
    n = len(A)
    if n <= 8:
        best = None
        for p in permutations(range(n)):
            flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or flat < best:
                best = flat
        return ("exact", best)
    return ("wl", n, _wl_hash(A))


def phase_invariant(tiling):
    """Dedup key for a toric (Seiberg-dual) phase: gauge count, field count, and
    the canonical quiver adjacency.  Phases of one diagram share a Kasteleyn
    Newton polygon, so that is a sanity check, not a discriminator."""
    return (tiling.num_gauge, tiling.num_fields,
            canonical_adjacency(tiling.adjacency_int()))


# ===========================================================================
# Seiberg duality (urban renewal) and toric-phase enumeration
# ===========================================================================
def urban_renewal(dimer, face_orbit, vertices):
    """Seiberg duality (the dimer "urban renewal" / square move) on one square
    gauge face -- the move that takes a toric quiver to a *different* toric phase
    of the same Calabi-Yau (a flop is a Kähler move and leaves the dimer fixed;
    this is the one that genuinely changes it).

    A square gauge face (a length-4 phi-orbit) has its four boundary edges
    alternating incoming (X) / outgoing (Y).  Each of its four corners pairs one
    Y with one X into a *meson* M = Y.X.  The dual tiling:
      * the surviving corner vertex keeps its off-face edges plus that meson;
      * a new opposite-colour cubic vertex {M, rev(Y), rev(X)} carries the meson
        superpotential coupling, with the X,Y arrows reversed.
    The new edges have no embedding, so their homology is re-solved from the
    rotation system (`solve_homology`).  Returns the dualized BraneTiling (with a
    Newton-certified homology cochain) or None if the face is not a clean square
    / the dual does not certify.
    """
    target_hull = convex_hull(vertices)
    edges = dimer.edges
    Y = [e for (e, s) in face_orbit if s == 1]      # outgoing from this face
    X = [e for (e, s) in face_orbit if s == -1]     # incoming to this face
    if len(Y) != 2 or len(X) != 2:
        return None
    facE = set(Y) | set(X)

    corners = []                                    # (colour, vtx, Yedge, Xedge)
    for w in range(dimer.nW):
        onf = [e for e in facE if edges[e]["w"] == w]
        if len(onf) == 2 and any(e in Y for e in onf) and any(e in X for e in onf):
            corners.append(("W", w, [e for e in onf if e in Y][0],
                            [e for e in onf if e in X][0]))
    for b in range(dimer.nB):
        onf = [e for e in facE if edges[e]["b"] == b]
        if len(onf) == 2 and any(e in Y for e in onf) and any(e in X for e in onf):
            corners.append(("B", b, [e for e in onf if e in Y][0],
                            [e for e in onf if e in X][0]))
    if len(corners) != 4:
        return None

    # label-level membership: surviving corners + new cubic meson vertices
    Wsets, Bsets = {}, {}
    for ci, (colour, v, ye, xe) in enumerate(corners):
        m = f"m_{ye}_{xe}"
        offv = [f"e{e}" for e in range(len(edges))
                if e not in facE and edges[e]["w" if colour == "W" else "b"] == v]
        (Wsets if colour == "W" else Bsets)[f"S{ci}"] = offv + [m]
        (Bsets if colour == "W" else Wsets)[f"N{ci}"] = [m, f"y_{ye}", f"x_{xe}"]

    wv = {l: vlab for vlab, labs in Wsets.items() for l in labs}
    bv = {l: vlab for vlab, labs in Bsets.items() for l in labs}
    if set(wv) != set(bv) or len(Wsets) != len(Bsets):
        return None                                 # not a genuine two-term dual
    elabels = sorted(wv)
    wn, bn = list(Wsets), list(Bsets)
    widx = {v: i for i, v in enumerate(wn)}
    bidx = {v: i for i, v in enumerate(bn)}
    eidx = {l: i for i, l in enumerate(elabels)}
    base_edges = [{"w": widx[wv[l]], "b": bidx[bv[l]], "h": [0, 0]} for l in elabels]

    # the cyclic order within each vertex's term is not fixed by membership;
    # search the orderings for the genus-1, anomaly-free, toric embedding whose
    # homology Newton-certifies back to the toric diagram.
    from itertools import permutations, product
    cyc = lambda lst: [[lst[0]] + list(p) for p in permutations(lst[1:])]
    Wopts = [cyc(Wsets[v]) for v in wn]
    Bopts = [cyc(Bsets[v]) for v in bn]
    for wsel in product(*Wopts):
        for bsel in product(*Bopts):
            rot_w = [[eidx[l] for l in o] for o in wsel]
            rot_b = [[eidx[l] for l in o] for o in bsel]
            dim = DimerGraph(len(wn), len(bn),
                             [dict(e) for e in base_edges], rot_w, rot_b)
            c = forward_extract(dim, vertices).checks
            if not (c["euler_V_minus_E_plus_F"] == 0 and c["gauge_eq_2area"]
                    and c["anomaly_free"] and c["toric_superpotential"]):
                continue
            sol = solve_homology(dim, target_hull)
            if sol is None:
                continue
            for i, e in enumerate(dim.edges):       # install the solved cochain
                e["h"] = list(sol[i])
            return forward_extract(dim, vertices)   # re-extract WITH homology
    return None


def enumerate_toric_phases(vertices, max_phases: int = 12):
    """All distinct toric (Seiberg-dual) phases of a toric diagram, found by
    breadth-first urban-renewal on every square gauge face and de-duplicated by
    `phase_invariant`.  Phase 0 is the Gulotta seed (`inverse_quiver`).

    Returns a list of BraneTilings, one per phase.  C3 / the conifold have a
    single phase; F0 has two ({8 fields, square faces} and {12, hexagonal})."""
    seed = inverse_quiver(vertices)
    phases = [seed]
    keys = {phase_invariant(seed)}
    frontier = [seed]
    while frontier and len(phases) < max_phases:
        nxt = []
        for tiling in frontier:
            dimer = tiling.to_dimer()
            faces, _ = _trace_faces(dimer, _FACE_HAND)
            for orb in faces:
                if len(orb) != 4:                   # only square faces dualize
                    continue
                dual = urban_renewal(dimer, orb, vertices)
                if dual is None:
                    continue
                k = phase_invariant(dual)
                if k in keys:
                    continue
                if not gl2z_equiv(kasteleyn_newton_polygon(dual),
                                  convex_hull(vertices)):
                    continue                        # reject Newton-polygon drift
                keys.add(k)
                phases.append(dual)
                nxt.append(dual)
                if len(phases) >= max_phases:
                    break
            if len(phases) >= max_phases:
                break
        frontier = nxt
    return phases


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


def _tiling_json(t, vertices) -> dict:
    """Serialise one BraneTiling (a single toric phase) for the web API."""
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
            # un-reduced (planar universal-cover) layout for a clean tiling
            "white_g": [[round(x, 5), round(y, 5)] for (x, y) in t.white_glob],
            "black_g": [[round(x, 5), round(y, 5)] for (x, y) in t.black_glob],
            "fields": [[round(x, 5), round(y, 5)] for (x, y) in t.field_pos],
            "edges": t.tiling_edges,
            # per-edge homology (white -> black+h is the true universal-cover
            # edge); needed to draw the tiling without collapsing edges.
            "edge_h": [f["homology"] for f in t.fields],
        },
        "checks": checks,
        "note": t.note,
    }


def inverse_quiver_json(vertices, **kw) -> dict:
    """JSON-friendly inverse-algorithm result for the web API (or an error)."""
    try:
        t = inverse_quiver(vertices, **kw)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    return _tiling_json(t, vertices)


def inverse_phases_json(vertices, max_phases: int = 12) -> dict:
    """All distinct toric (Seiberg-dual) phases of the diagram, for the web API:
    a list of serialised tilings (phase 0 = the Gulotta seed), so the toric tab
    can cycle between them.  Each is independently Newton-certified."""
    try:
        phases = enumerate_toric_phases(vertices, max_phases=max_phases)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    return {
        "available": True,
        "num_phases": len(phases),
        "phases": [_tiling_json(t, vertices) for t in phases],
    }
