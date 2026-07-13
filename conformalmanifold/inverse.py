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

    Implementation note: the parameters are exact rationals, but doing the inner
    probe in `fractions.Fraction` was the dominant cost of the whole inverse
    algorithm (millions of Fraction ops).  Here the hot loop is pure integer
    arithmetic over a fixed common denominator per zig-zag pair, and Fractions are
    built only for the final |det| crossings -- bit-identical output, ~100x faster.
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
            # base-point offset bases[l]-bases[k], cleared to a common denom Q so
            # every probe stays integer; s,t are then numerator/(Q*-D) (mod 1).
            dx0 = bases[l][0] - bases[k][0]
            dy0 = bases[l][1] - bases[k][1]
            Q = dx0.denominator * dy0.denominator // gcd(dx0.denominator,
                                                         dy0.denominator)
            Ax = dx0.numerator * (Q // dx0.denominator)
            Ay = dy0.numerator * (Q // dy0.denominator)
            Dden = -D * Q                   # s = Ns/Dden, t = Nt/Dden (mod 1)
            Dabs = abs(Dden)
            sgn = 1 if Dden > 0 else -1
            found = {}                      # (s_int, t_int) -> kept (insertion order)
            rng = 4 + max(abs(pk), abs(qk), abs(pl), abs(ql))
            while len(found) < need:
                for m in range(-rng, rng + 1):
                    rxN = Ax + m * Q
                    for nn in range(-rng, rng + 1):
                        ryN = Ay + nn * Q
                        ns = -ql * rxN + pl * ryN
                        nt = pk * ryN - qk * rxN
                        key = ((sgn * ns) % Dabs, (sgn * nt) % Dabs)
                        if key not in found:
                            found[key] = True
                rng += 3
                if rng > 60:                 # safety; should never trigger
                    break
            for (s_int, t_int) in found:
                out.append((k, l, Fr(s_int, Dabs), Fr(t_int, Dabs)))
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
    # trace-face index -> displayed gauge-node id.  None = identity.  Set by
    # `urban_renewal`, which relabels the dual's nodes to match their heritage
    # (spectator nodes keep their ids; the dualized node keeps its id) so a
    # Seiberg duality on "node k" visibly acts AT node k.
    face_to_node: list = None

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
    the canonical quiver adjacency *up to charge conjugation* (reversing every
    arrow = transposing the adjacency is the same physical theory, and the
    standard phase counting -- e.g. dP2: 2 phases, dP3: 4 -- identifies the
    conjugate pair).  Phases of one diagram share a Kasteleyn Newton polygon,
    so that is a sanity check, not a discriminator."""
    A = tiling.adjacency_int()
    n = len(A)
    At = [[A[j][i] for j in range(n)] for i in range(n)]
    return (tiling.num_gauge, tiling.num_fields,
            min(canonical_adjacency(A), canonical_adjacency(At)))


# ===========================================================================
# Seiberg duality (urban renewal) and toric-phase enumeration
# ===========================================================================
def integrate_masses(dimer, return_map=False):
    """Integrate out massive fields: contract every 2-valent tiling vertex.

    A 2-valent vertex is a quadratic superpotential term  W ⊃ e1·e2  -- a mass
    pair for the two fields e1, e2.  Integrating them out (setting the F-term
    e1 ∝ e2 and substituting) deletes the vertex and both edges and *merges*
    the two opposite-colour neighbours into one vertex, splicing their rotation
    orders where the deleted edges sat:

        rot(merged) = rot(v1) cut open at e1  +  rot(v2) cut open at e2.

    Repeats until no 2-valent vertex remains (a contraction can create new mass
    pairs).  Returns the input dimer unchanged (homology intact) if there was
    nothing to do; otherwise a new DimerGraph with all edge homologies reset to
    (0,0) -- re-solve with `solve_homology`.  Returns None on a degenerate
    contraction (double edge at the mass vertex, digon neighbour pair, or a
    vertex left with fewer than 2 edges).

    With `return_map=True` returns (dimer, edge_map) instead, where edge_map
    maps each SURVIVING input edge index to its output index (integrated-out
    edges are absent) -- used by `urban_renewal` to track gauge-face heritage
    through the move."""
    edges = {i: dict(e) for i, e in enumerate(dimer.edges)}
    rot = {("W", i): list(r) for i, r in enumerate(dimer.rot_w)}
    rot.update({("B", i): list(r) for i, r in enumerate(dimer.rot_b)})
    changed = False
    while True:
        mv = next((v for v, r in rot.items() if len(r) == 2), None)
        if mv is None:
            break
        changed = True
        e1, e2 = rot[mv]
        if e1 == e2:
            return None                       # double edge: not a mass pair
        col = mv[0]                           # colour of the 2-valent vertex
        okey = "b" if col == "W" else "w"     # endpoint key of its neighbours
        ncol = "B" if col == "W" else "W"
        v1 = (ncol, edges[e1][okey])
        v2 = (ncol, edges[e2][okey])
        if v1 == v2:
            return None                       # digon: graph would degenerate
        r1, r2 = rot[v1], rot[v2]
        i1, i2 = r1.index(e1), r2.index(e2)
        merged = r1[i1 + 1:] + r1[:i1] + r2[i2 + 1:] + r2[:i2]
        if len(merged) < 2:
            return None
        for e in r2:                          # re-hang v2's edges on v1
            if e != e2:
                edges[e][okey] = v1[1]
        del rot[mv], rot[v2], edges[e1], edges[e2]
        rot[v1] = merged
    if not changed:
        return (dimer, {i: i for i in range(len(dimer.edges))}) if return_map \
            else dimer
    wids = sorted(i for (c, i) in rot if c == "W")
    bids = sorted(i for (c, i) in rot if c == "B")
    wmap = {old: new for new, old in enumerate(wids)}
    bmap = {old: new for new, old in enumerate(bids)}
    eids = sorted(edges)
    emap = {old: new for new, old in enumerate(eids)}
    new_edges = [{"w": wmap[edges[i]["w"]], "b": bmap[edges[i]["b"]],
                  "h": [0, 0]} for i in eids]
    rot_w = [[emap[e] for e in rot[("W", i)]] for i in wids]
    rot_b = [[emap[e] for e in rot[("B", i)]] for i in bids]
    out = DimerGraph(nW=len(wids), nB=len(bids),
                     edges=new_edges, rot_w=rot_w, rot_b=rot_b)
    return (out, emap) if return_map else out


def urban_renewal(dimer, face_orbit, vertices, node_of_face=None):
    """Seiberg duality (the dimer "urban renewal" / square move) on one square
    gauge face -- the move that takes a toric quiver to a *different* toric phase
    of the same Calabi-Yau (a flop is a Kähler move and leaves the dimer fixed;
    this is the one that genuinely changes it).

    Gauge nodes keep their IDENTITY through the move: every gauge face other
    than the dualized one persists (retaining its off-square boundary edges),
    so the dual's faces are matched back to the originals through the
    surviving edges and the returned tiling's nodes are relabelled to their
    heritage ids -- spectator nodes keep their labels, the dualized node keeps
    its label with reversed flavors, and the mesons connect its neighbours.
    (Without this the trace-order renumbering scrambled the labels: dualizing
    node 2 of F0 displayed the same LABELLED quiver as dualizing node 0.)
    `node_of_face` optionally maps this dimer's trace-face indices to display
    ids (composing heritage along a duality path); the result's map is stored
    in `BraneTiling.face_to_node`.

    A square gauge face (a length-4 phi-orbit) has its four boundary edges
    alternating incoming (X) / outgoing (Y).  Each of its four corners pairs one
    Y with one X into a *meson* M = Y.X.  The dual tiling:
      * the surviving corner vertex keeps its off-face edges plus that meson;
      * a new opposite-colour cubic vertex {M, rev(Y), rev(X)} carries the meson
        superpotential coupling, with the X,Y arrows reversed.
    The move is local, but the dual dimer carries the WHOLE graph: vertices not
    on the face ("spectators") keep their rotation orders verbatim, and each
    surviving corner gets the meson spliced into its rotation exactly where the
    adjacent (Y, X) pair sat.  Corners are the four *dart transitions* of the
    face orbit -- not four distinct vertices: a face may visit the same vertex
    twice (SPP's square face does), in which case that vertex gets two mesons.
    The cyclic orientation of the four new cubic vertices is fixed by the
    handedness of the rotation system (for _FACE_HAND = +1 the cubic at a white
    corner reads [M, rev X, rev Y], at a black corner [M, rev Y, rev X]); the
    derived orientation is tried first and the other 15 combinations remain
    only as a defensive fallback, gated by the Newton certificate.  If a
    surviving corner is left 2-valent the meson pair is *massive*;
    `integrate_masses` contracts it, so the returned phase is always the
    fully-integrated IR theory (this is what produces the extra dP2/dP3
    phases).  The new edges have no embedding, so their homology is re-solved
    from the rotation system (`solve_homology`).  Returns the dualized
    BraneTiling (with a Newton-certified homology cochain) or None if the face
    is not a clean square / the dual does not certify.
    """
    target_hull = convex_hull(vertices)
    edges = dimer.edges
    orb = list(face_orbit)
    facE = {e for (e, _) in orb}             # the four square-face edges
    if len(orb) != 4 or len(facE) != 4:
        return None
    if sorted(s for (_, s) in orb) != [-1, -1, 1, 1]:
        return None

    # corners = the four dart transitions of the phi-orbit.  phi((e, s)) turns
    # around the white endpoint of e when s = -1 and the black one when s = +1,
    # so the vertex shared by consecutive darts -- and its rotation-adjacent
    # (e1, e2) edge pair -- is read straight off the orbit.
    corners = []                             # (colour, vtx, e1, e2, ye, xe)
    for i in range(4):
        e1, s1 = orb[i]
        e2, s2 = orb[(i + 1) % 4]
        if s2 != -s1:
            return None                      # not an alternating gauge face
        colour = "W" if s1 == -1 else "B"
        v = edges[e1]["w" if colour == "W" else "b"]
        ye, xe = (e1, e2) if s1 == 1 else (e2, e1)
        corners.append((colour, v, e1, e2, ye, xe))

    # label-level rotations.  Spectator vertices are untouched by the move;
    # each corner occurrence has its adjacent (e1, e2) pair replaced -- in
    # place -- by that corner's meson (in the W term ...Y.X... becomes ...M...).
    occ = {}                                 # (colour, v) -> corner indices
    for ci, (colour, v, *_rest) in enumerate(corners):
        occ.setdefault((colour, v), []).append(ci)
    Wrots, Brots = {}, {}                    # vertex label -> cyclic label list
    for w in range(dimer.nW):
        if ("W", w) not in occ:
            Wrots[f"V_W{w}"] = [f"e{e}" for e in dimer.rot_w[w]]
    for b in range(dimer.nB):
        if ("B", b) not in occ:
            Brots[f"V_B{b}"] = [f"e{e}" for e in dimer.rot_b[b]]
    for (colour, v), cis in occ.items():
        rot = dimer.rot_w[v] if colour == "W" else dimer.rot_b[v]
        k = len(rot)
        if k - len(cis) < 2:                 # survivor would drop below 2-valent
            return None
        meson_at, skip = {}, set()
        for ci in cis:
            _, _, e1, e2, _, _ = corners[ci]
            p = rot.index(e1)
            if rot[(p + 1) % k] != e2:       # pair not rotation-adjacent here
                return None
            meson_at[p] = f"m{ci}"
            skip.add((p + 1) % k)
        surv = [meson_at.get(p, f"e{rot[p]}") for p in range(k)
                if p not in skip]
        (Wrots if colour == "W" else Brots)[f"S_{colour}{v}"] = surv
    cubics = [(f"N{ci}", "B" if colour == "W" else "W",
               f"m{ci}", f"y_{ye}", f"x_{xe}")
              for ci, (colour, v, e1, e2, ye, xe) in enumerate(corners)]

    # cubic vertex {M, rev(Y), rev(X)} orientation: the handedness-derived
    # combination first (o = 1 at white corners, 0 at black), then the
    # defensive sweep of the remaining 15.
    from itertools import product
    derived = tuple(1 if c[0] == "W" else 0 for c in corners)
    orients = [derived] + [o for o in product((0, 1), repeat=4) if o != derived]
    for orient in orients:
        Wr = dict(Wrots)
        Br = dict(Brots)
        for (lab, colour, m, y, x), o in zip(cubics, orient):
            order = [m, y, x] if o == 0 else [m, x, y]
            (Wr if colour == "W" else Br)[lab] = order

        wv = {l: vlab for vlab, labs in Wr.items() for l in labs}
        bv = {l: vlab for vlab, labs in Br.items() for l in labs}
        if set(wv) != set(bv) or len(Wr) != len(Br):
            return None                      # not a genuine two-term dual
        wn, bn = list(Wr), list(Br)
        widx = {v: i for i, v in enumerate(wn)}
        bidx = {v: i for i, v in enumerate(bn)}
        elabels = sorted(wv)
        eidx = {l: i for i, l in enumerate(elabels)}
        base_edges = [{"w": widx[wv[l]], "b": bidx[bv[l]], "h": [0, 0]}
                      for l in elabels]
        rot_w = [[eidx[l] for l in Wr[v]] for v in wn]
        rot_b = [[eidx[l] for l in Br[v]] for v in bn]
        dim = DimerGraph(len(wn), len(bn), base_edges, rot_w, rot_b)
        dim, emap = integrate_masses(dim, return_map=True)   # kill mass terms
        if dim is None:
            continue
        c = forward_extract(dim, vertices).checks
        if not (c["euler_V_minus_E_plus_F"] == 0 and c["gauge_eq_2area"]
                and c["anomaly_free"] and c["toric_superpotential"]):
            continue
        sol = solve_homology(dim, target_hull)
        if sol is None:
            continue
        for i, e in enumerate(dim.edges):    # install the solved cochain
            e["h"] = list(sol[i])
        t = forward_extract(dim, vertices)   # re-extract WITH homology
        _harmonic_torus_layout(t)            # true flat-torus embedding

        # --- gauge-node heritage: match dual faces to original faces --------
        # every surviving off-square edge keeps its two face-sides, so its
        # darts vote for the correspondence; the one dual face with no
        # surviving dart is the dualized node itself.
        ofaces, ofaceof = _trace_faces(dimer, _FACE_HAND)
        dfaces, dfaceof = _trace_faces(dim, _FACE_HAND)
        f0 = ofaceof[orb[0]]                 # the dualized face's original id
        match, ok = {}, True
        for n in range(len(edges)):
            if n in facE:
                continue
            j = emap.get(eidx.get(f"e{n}"))
            if j is None:
                continue                     # edge eaten by mass integration
            for s in (1, -1):
                of = ofaceof[(n, s)]
                df = dfaceof[(j, s)]
                if match.setdefault(df, of) != of:
                    ok = False               # conflicting votes: bail out
        missing_d = [j for j in range(len(dfaces)) if j not in match]
        missing_o = [i for i in range(len(ofaces)) if i not in match.values()]
        if ok and missing_d == [] and missing_o == []:
            pass
        elif ok and len(missing_d) == 1 and missing_o == [f0]:
            match[missing_d[0]] = f0
        else:
            ok = False
        if ok and len(set(match.values())) == len(dfaces) == len(ofaces):
            heritage = (node_of_face if node_of_face is not None
                        else list(range(len(ofaces))))
            perm = [heritage[match[j]] for j in range(len(dfaces))]
            n_g = t.num_gauge
            A = t.adjacency_int()
            A2 = [[0] * n_g for _ in range(n_g)]
            for i in range(n_g):
                for j in range(n_g):
                    A2[perm[i]][perm[j]] = A[i][j]
            t.adjacency = A2
            for f in t.fields:
                f["src"], f["tgt"] = perm[f["src"]], perm[f["tgt"]]
            t.face_to_node = perm
        # (on a failed match the tiling stays in trace order -- still a correct
        # phase, just without heritage labels)
        return t
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


def _trace_strands(dimer):
    """The ZIG-ZAG strands of a DimerGraph (turn to the next edge around
    black, to the previous around white): returns (windings, edge_strands)
    where windings[i] is strand i's homology winding and edge_strands[e] the
    (two) strands through edge e.  Each strand is the alga image of one
    primitive boundary segment of the toric diagram (Hanany-Vegh)."""
    posw = [{e: i for i, e in enumerate(r)} for r in dimer.rot_w]
    posb = [{e: i for i, e in enumerate(r)} for r in dimer.rot_b]
    seen = set()
    windings = []
    edge_strands = [[] for _ in dimer.edges]
    for e0 in range(len(dimer.edges)):
        for s0 in (1, -1):
            if (e0, s0) in seen:
                continue
            sid = len(windings)
            wx = wy = 0
            e, s = e0, s0
            while (e, s) not in seen:
                seen.add((e, s))
                edge_strands[e].append(sid)
                ed = dimer.edges[e]
                if s == 1:                       # traverse white -> black: +h
                    wx += ed["h"][0]
                    wy += ed["h"][1]
                    b = ed["b"]
                    r = dimer.rot_b[b]
                    e = r[(posb[b][e] + 1) % len(r)]     # next around black
                    s = -1
                else:                            # black -> white: -h
                    wx -= ed["h"][0]
                    wy -= ed["h"][1]
                    w = ed["w"]
                    r = dimer.rot_w[w]
                    e = r[(posw[w][e] - 1) % len(r)]     # previous around white
                    s = 1
            windings.append((wx, wy))
    return windings, edge_strands


def _strand_polygon(dimer):
    """The Newton polygon read off the zig-zag strand windings (rotated by 90
    degrees and chained in angular order).  Linear in the number of edges --
    an exact certificate at any size (unlike the O(n!) Kasteleyn permanent).
    Returns the polygon (up to translation) or None if inconsistent."""
    windings, _ = _trace_strands(dimer)
    windings = [w for w in windings if w != (0, 0)]
    if not windings:
        return None
    vecs = sorted(((-wy, wx) for (wx, wy) in windings),
                  key=lambda v: math.atan2(v[1], v[0]))
    if (sum(v[0] for v in vecs), sum(v[1] for v in vecs)) != (0, 0):
        return None
    pts, x, y = [(0, 0)], 0, 0
    for (dx, dy) in vecs:
        x += dx
        y += dy
        pts.append((x, y))
    return convex_hull(pts)


def _assign_strand_legs(hull, windings):
    """Map each zig-zag strand to a LEG index (0..B-1 in the hull's CCW
    boundary-leg order, the convention of `scft.toric_field_R_charges`).

    Both the hull legs and the strands trace the same polygon up to GL(2,Z),
    so their edge groups match cyclically; the alignment is pinned by the
    (GL(2,Z)-invariant, up to overall reflection) sequences of group sizes and
    consecutive-direction determinants.  Legs within one polygon edge are
    interchangeable (the arc gaps between them vanish).  Returns
    {strand_id: leg_index} or None."""
    n = len(hull)
    dirs, lens = [], []
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        g = gcd(abs(dx), abs(dy)) or 1
        dirs.append((dx // g, dy // g))
        lens.append(g)
    hull_dets = [dirs[i][0] * dirs[(i + 1) % n][1]
                 - dirs[i][1] * dirs[(i + 1) % n][0] for i in range(n)]
    starts = [sum(lens[:i]) for i in range(n)]

    # strand groups in CCW angular order of their polygon edge vectors
    order = sorted(range(len(windings)),
                   key=lambda i: math.atan2(windings[i][0], -windings[i][1]))
    groups = []                                  # [(edge_vec_dir, [strand ids])]
    for i in order:
        wx, wy = windings[i]
        v = (-wy, wx)
        g = gcd(abs(v[0]), abs(v[1])) or 1
        d = (v[0] // g, v[1] // g)
        if groups and groups[-1][0] == d:
            groups[-1][1].append(i)
        else:
            groups.append((d, [i]))
    if len(groups) != n:
        return None
    gdirs = [g[0] for g in groups]
    gdets = [gdirs[i][0] * gdirs[(i + 1) % n][1]
             - gdirs[i][1] * gdirs[(i + 1) % n][0] for i in range(n)]
    gsizes = [len(g[1]) for g in groups]

    for orient in (1, -1):                       # det(U) = -1 reverses traversal
        for r in range(n):
            idx = [(r + orient * i) % n for i in range(n)]
            if any(gsizes[idx[i]] != lens[i] for i in range(n)):
                continue
            seq = [gdirs[idx[i]][0] * gdirs[idx[(i + 1) % n]][1]
                   - gdirs[idx[i]][1] * gdirs[idx[(i + 1) % n]][0]
                   for i in range(n)]
            if seq != hull_dets and seq != [-d for d in hull_dets]:
                continue
            out = {}
            for i in range(n):
                for j, sid in enumerate(groups[idx[i]][1]):
                    out[sid] = starts[i] + j
            return out
    return None


def _orbifold_honeycomb(vertices):
    """Direct brane tiling for an abelian orbifold C^3/Gamma -- any lattice
    TRIANGLE diagram (|Gamma| = 2*area): the honeycomb dimer of C^3 on a
    quotient torus R^2/Lambda for an index-|Gamma| sublattice Lambda of Z^2.
    Exact -- no placement search (the random Gulotta search fails on e.g.
    Z3xZ3 / Z4xZ4), exact edge homology, geometric honeycomb embedding.
    Returns a BraneTiling or None (not a triangle / no sublattice certifies).

    Construction: the C^3 tiling is the honeycomb with one white + one black
    vertex per cell c of Z^2 and edges W(c)-B(c), W(c)-B(c+x), W(c)-B(c+y)
    (Kasteleyn det 1 + x + y = the unit triangle).  Quotienting by Lambda
    keeps one cell per coset of Z^2/Lambda; an edge into a neighbouring coset
    carries homology = the Lambda-coordinates of the wrap-around displacement.
    The right Lambda is found by enumerating the index-k sublattices (Hermite
    forms (a,0),(b,d), ad = k) and keeping the one whose zig-zag strand
    polygon (`_strand_polygon` -- cheap and exact at any size) matches the
    input triangle up to GL(2,Z)."""
    from fractions import Fraction as F
    from math import floor

    hull = convex_hull(vertices)
    if len(hull) != 3:
        return None
    k = normalized_area(hull)

    def attempt(v1, v2):
        det = v1[0] * v2[1] - v1[1] * v2[0]

        def coords(g):                   # g = s*v1 + t*v2 (exact rationals)
            s = F(g[0] * v2[1] - g[1] * v2[0], det)
            t = F(v1[0] * g[1] - v1[1] * g[0], det)
            return s, t

        def rep(g):                      # canonical coset representative
            s, t = coords(g)
            ds, dt = floor(s), floor(t)
            return (g[0] - ds * v1[0] - dt * v2[0],
                    g[1] - ds * v1[1] - dt * v2[1])

        xs = [0, v1[0], v2[0], v1[0] + v2[0]]
        ys = [0, v1[1], v2[1], v1[1] + v2[1]]
        cells = sorted({rep((x, y))
                        for x in range(min(xs), max(xs) + 1)
                        for y in range(min(ys), max(ys) + 1)})
        if len(cells) != k:
            return None
        idx = {c: i for i, c in enumerate(cells)}

        DELTAS = [(0, 0), (1, 0), (0, 1)]
        edges = []
        for i, g in enumerate(cells):
            for (dx, dy) in DELTAS:
                tgt = (g[0] + dx, g[1] + dy)
                r = rep(tgt)
                D = (tgt[0] - r[0], tgt[1] - r[1])       # in Lambda
                hs, ht = coords(D)                       # integer coordinates
                # sign convention: the true universal-cover edge runs from the
                # white vertex to the black vertex at (black_glob - h), which
                # is what the tiling renderer draws.
                edges.append({"w": i, "b": idx[r],
                              "h": [-int(hs), -int(ht)]})
        rot_w = [[3 * i + d for d in range(3)] for i in range(k)]

        for order in ((0, 2, 1), (0, 1, 2)):   # handedness fixed by the checks
            rot_b = [[] for _ in range(k)]
            for j, g in enumerate(cells):
                for d in order:                          # cyclic order around B
                    dx, dy = DELTAS[d]
                    src = rep((g[0] - dx, g[1] - dy))
                    rot_b[j].append(3 * idx[src] + d)
            dim = DimerGraph(nW=k, nB=k, edges=[dict(e) for e in edges],
                             rot_w=rot_w, rot_b=rot_b)
            poly = _strand_polygon(dim)
            if poly is None or not gl2z_equiv(poly, hull):
                continue                       # wrong sublattice / handedness
            t = forward_extract(dim, vertices)
            c = t.checks
            if not (c["euler_V_minus_E_plus_F"] == 0 and c["gauge_eq_2area"]
                    and c["anomaly_free"] and c["toric_superpotential"]):
                continue
            # geometric honeycomb embedding on the quotient torus: positions
            # are PARENT honeycomb positions (white at cell + (2/3,2/3), black
            # at cell + (1/3,1/3)) mapped through the torus frame -- a linear
            # image of the planar honeycomb, so edges never cross.
            def torus(px, py):
                s = F(px * v2[1] - py * v2[0], det)
                t2 = F(v1[0] * py - v1[1] * px, det)
                return float(s), float(t2)

            wglob = [torus(g[0] + F(2, 3), g[1] + F(2, 3)) for g in cells]
            bglob = [torus(g[0] + F(1, 3), g[1] + F(1, 3)) for g in cells]
            t.white_glob, t.black_glob = wglob, bglob
            t.white_pos = [(x % 1, y % 1) for (x, y) in wglob]
            t.black_pos = [(x % 1, y % 1) for (x, y) in bglob]
            fpos = []
            for e in dim.edges:
                wx, wy = t.white_glob[e["w"]]
                bx = t.black_glob[e["b"]][0] - e["h"][0]   # true edge endpoint
                by = t.black_glob[e["b"]][1] - e["h"][1]
                fpos.append((((wx + bx) / 2) % 1, ((wy + by) / 2) % 1))
            t.field_pos = fpos
            # per-field zig-zag leg pair (feeds the R-charge computation)
            windings, edge_strands = _trace_strands(dim)
            legof = _assign_strand_legs(hull, windings)
            if legof is not None:
                for e, f in enumerate(t.fields):
                    s1, s2 = edge_strands[e]
                    f["zigzag"] = sorted((legof[s1], legof[s2]))
            t.note = (f"abelian orbifold C^3/Gamma, |Gamma| = {k}: exact "
                      "honeycomb dimer on the quotient torus (no placement "
                      "search; zig-zag strand polygon certified).")
            return t
        return None

    def gauss_reduce(u, v):
        """Shortest basis of the lattice <u, v> (Lagrange/Gauss reduction) --
        same lattice, so same tiling, but a near-regular drawn honeycomb
        instead of a heavily sheared one."""
        u, v = tuple(u), tuple(v)
        while True:
            if u[0] * u[0] + u[1] * u[1] > v[0] * v[0] + v[1] * v[1]:
                u, v = v, u
            n = u[0] * u[0] + u[1] * u[1]
            if n == 0:
                return u, v
            q = round((u[0] * v[0] + u[1] * v[1]) / n)
            if q == 0:
                return u, v
            v = (v[0] - q * u[0], v[1] - q * u[1])

    # enumerate index-k sublattices in Hermite form: (a, 0), (b, d), ad = k
    for a in range(1, k + 1):
        if k % a:
            continue
        d = k // a
        for b in range(a):
            t = attempt(*gauss_reduce((a, 0), (b, d)))
            if t is not None:
                return t
    return None


def inverse_quiver(vertices, max_attempts: int = 400, max_gauge: int = 60):
    """Reconstruct a quiver gauge theory + brane tiling from a toric diagram.

    `vertices` : lattice points / corners of the toric diagram (any order;
                 the convex hull is taken).
    Triangle diagrams (= abelian orbifolds C^3/Gamma) are built EXACTLY via the
    quotient honeycomb (`_orbifold_honeycomb`); everything else goes through
    Gulotta's properly-ordered-dimer placement search.  Returns a
    `BraneTiling`, or raises `ValueError` if the diagram is degenerate or no
    consistent placement is found within `max_attempts`.
    """
    hull = convex_hull(vertices)
    if len(hull) < 3:
        raise ValueError("need at least 3 non-collinear lattice points")
    a2 = normalized_area(hull)
    if a2 > max_gauge:
        raise ValueError(f"toric diagram too large (2*area = {a2} gauge nodes > "
                         f"max_gauge = {max_gauge})")
    if len(hull) == 3:
        t = _orbifold_honeycomb(hull)
        if t is not None:
            return t
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


# The Kasteleyn Newton-polygon certificate expands a permanent over n_white!
# permutations -- a cheap, exact cross-check for small tilings, but it blows up
# (~90ms at 8 white nodes, ~0.8s at 9, then exponentially) and would stall the
# interactive toric tab on large blown-up diagrams.  It is only an INDEPENDENT
# certificate of the (already self-consistent) reconstruction, so skip it past
# this size; the quiver, superpotential and the tiling's own consistency checks
# are unaffected.
KAST_NEWTON_MAX_WHITE = 8


def _harmonic_torus_layout(t):
    """Install a HARMONIC (Tutte) embedding on the flat torus for a mutated
    tiling: its fields carry solved homology classes, and the renderer draws
    each universal-cover edge from the white vertex to (black - h), so
    minimising  sum_e |x_w - x_b + h_e|^2  (a Laplacian linear system with one
    vertex pinned) places every vertex at the centroid of its true neighbour
    images.  This is the natural flat-torus drawing -- crossing-free for the
    consistent tilings produced here -- and replaces the schematic
    spanning-tree layout, which ignored the homology and drew a tangle."""
    nW, nB = t.num_white, t.num_black
    n = nW + nB
    if n < 2:
        return
    A = [[0.0] * n for _ in range(n)]
    rx = [0.0] * n
    ry = [0.0] * n
    for f in t.fields:
        w, b = f["white"], nW + f["black"]
        hx, hy = f["homology"]
        A[w][w] += 1.0
        A[b][b] += 1.0
        A[w][b] -= 1.0
        A[b][w] -= 1.0
        rx[w] -= hx
        ry[w] -= hy
        rx[b] += hx
        ry[b] += hy

    # pin vertex 0 at the origin; solve the reduced system by Gaussian
    # elimination (two right-hand sides, same matrix)
    m = n - 1
    M = [[A[i + 1][j + 1] for j in range(m)] + [rx[i + 1], ry[i + 1]]
         for i in range(m)]
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12:
            return                            # disconnected: keep old layout
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(m):
            if r != col and M[r][col] != 0.0:
                f0 = M[r][col]
                M[r] = [a - f0 * b for a, b in zip(M[r], M[col])]
    xs = [0.0] + [M[i][m] for i in range(m)]
    ys = [0.0] + [M[i][m + 1] for i in range(m)]

    t.white_glob = [(xs[i], ys[i]) for i in range(nW)]
    t.black_glob = [(xs[nW + i], ys[nW + i]) for i in range(nB)]
    t.white_pos = [(x % 1, y % 1) for (x, y) in t.white_glob]
    t.black_pos = [(x % 1, y % 1) for (x, y) in t.black_glob]
    fpos = []
    for f in t.fields:
        wx, wy = t.white_glob[f["white"]]
        bx = t.black_glob[f["black"]][0] - f["homology"][0]
        by = t.black_glob[f["black"]][1] - f["homology"][1]
        fpos.append((((wx + bx) / 2) % 1, ((wy + by) / 2) % 1))
    t.field_pos = fpos


def _normalize_tiling(t, vertices):
    """Round-trip a BraneTiling through its combinatorial map so its gauge-node
    numbering is the canonical face-trace order (the order `urban_renewal` and
    `square_gauge_faces` use) -- while keeping the original geometric embedding
    (to_dimer/forward_extract preserve vertex and edge indexing)."""
    r = forward_extract(t.to_dimer(), vertices)
    r.white_pos, r.black_pos = t.white_pos, t.black_pos
    r.white_glob, r.black_glob = t.white_glob, t.black_glob
    r.field_pos = t.field_pos
    # keep extra per-field annotations (e.g. the seed's zig-zag pair, which
    # feeds the R-charge computation); edge indexing is preserved.
    for old, new in zip(t.fields, r.fields):
        for key, val in old.items():
            if key not in new:
                new[key] = val
    r.checks["fields_eq_sum_det"] = t.checks.get("fields_eq_sum_det", True)
    r.note = t.note
    return r


def square_gauge_faces(t):
    """The gauge nodes available for a Seiberg-duality (urban renewal) move:
    the DISPLAYED node ids whose gauge face is square (a square face <=> that
    node has N_f = 2 N_c, the self-dual case toric duality acts on).  Trace
    order is mapped through the tiling's heritage labels (`face_to_node`)."""
    dimer = t.to_dimer()
    faces, _ = _trace_faces(dimer, _FACE_HAND)
    lab = t.face_to_node or list(range(len(faces)))
    return sorted(lab[i] for i, orb in enumerate(faces) if len(orb) == 4)


def face_polygons(t):
    """Per gauge face: its displayed node id, whether it is square (an
    available urban-renewal move), and its boundary polygon in the SAME
    universal-cover coordinates the web tiling drawing uses (white vertex at
    `white_glob[w]`, black at `black_glob[b] - h_e`), so the frontend can
    overlay clickable face regions on the drawn dimer.

    The polygon walks the face's dart orbit (phi = sigma o alpha, the trace
    `square_gauge_faces`/`urban_renewal` use), taking each dart's START vertex.
    Consecutive darts share a vertex; matching the two drawn representatives
    fixes the integer translate of each edge: across a shared black vertex the
    translate shifts by h_next - h_prev, across a shared white vertex it is
    unchanged.  The walk closes because the face-cocycle condition
    (sum of signed homologies around a gauge face = 0) telescopes exactly to
    the accumulated shift.  Faces are returned sorted by node id."""
    dimer = t.to_dimer()
    faces, _ = _trace_faces(dimer, _FACE_HAND)
    lab = t.face_to_node or list(range(len(faces)))
    wg, bg = t.white_glob, t.black_glob
    eh = [f["homology"] for f in t.fields]
    out = []
    for tf, orb in enumerate(faces):
        poly = []
        tx, ty = 0, 0                       # integer translate of current edge
        prev = None
        for (e, s) in orb:
            if prev is not None:
                pe, ps = prev
                if ps == 1:                 # shared vertex is black
                    tx += eh[e][0] - eh[pe][0]
                    ty += eh[e][1] - eh[pe][1]
            f = t.fields[e]
            if s == 1:                      # start = white vertex
                x, y = wg[f["white"]][0] + tx, wg[f["white"]][1] + ty
            else:                           # start = black vertex (drawn at -h)
                x = bg[f["black"]][0] - eh[e][0] + tx
                y = bg[f["black"]][1] - eh[e][1] + ty
            poly.append([round(x, 5), round(y, 5)])
            prev = (e, s)
        out.append({"node": lab[tf], "square": len(orb) == 4, "poly": poly})
    return sorted(out, key=lambda d: d["node"])


def dualize_path(vertices, path, start_phase: int = 0):
    """Apply a SEQUENCE of user-chosen urban-renewal (Seiberg duality) moves.

    `path` is a list of gauge-NODE ids (as displayed); each must be a square
    face (N_f = 2 N_c) of the tiling reached so far (starting from the seed,
    or from `enumerate_toric_phases(vertices)[start_phase]`).  Node identity
    is preserved through every move (heritage relabelling in `urban_renewal`),
    so "dualize node k" acts at the displayed node k at every step.  Returns
    the resulting BraneTiling; raises ValueError if an id is not a square face
    or the dual fails to Newton-certify."""
    state = seiberg_path(vertices, path, start_phase)
    if state.tiling is None:
        raise ValueError(state.reason or "path left the dimer regime")
    return state.tiling


def quiver_seiberg(A, ranks, k):
    """General quiver-level Seiberg duality on node k (ranks in units of N):
    N_c -> N_f - N_c, all flavors at k reversed, mesons A[i][k]*A[k][j] added
    between its neighbours, massive vector-like pairs integrated out.  Returns
    (A', ranks').  Raises ValueError if the dual rank would be < 1, or if node
    k carries an adjoint (self-loop): ordinary SU(N) Seiberg duality is not
    defined for a gauge node with adjoint matter (that is a Kutasov-type
    duality), and the naive meson/mass rule would otherwise emit a nonsensical
    quiver with negative arrow multiplicities."""
    n = len(A)
    if A[k][k] != 0:
        raise ValueError(
            f"node {k} has an adjoint (self-loop, multiplicity {A[k][k]}); "
            "ordinary Seiberg duality is undefined for a node with adjoint "
            "matter (this would be a Kutasov-type duality)")
    nf = sum(A[i][k] * ranks[i] for i in range(n) if i != k)
    nf_out = sum(A[k][j] * ranks[j] for j in range(n) if j != k)
    if nf != nf_out:
        raise ValueError(f"node {k} is anomalous (N_f in {nf} != out {nf_out})")
    rk = nf - ranks[k]
    if rk < 1:
        raise ValueError(
            f"Seiberg duality on node {k} gives rank {rk}N <= 0 "
            f"(N_f = {nf}N < 2 N_c) -- no dual gauge theory")
    B = [row[:] for row in A]
    for j in range(n):
        B[k][j], B[j][k] = A[j][k], A[k][j]          # reverse flavors at k
    # mesons M[ij] = A[i][k]*A[k][j] between k's neighbours (i,j != k).  Keep the
    # diagonal i==j: M[ii] is a genuine ADJOINT of node i when k has an arrow both
    # from and to i -- dropping it (an old `i != j` clause) under-counted fields.
    mes = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != k and j != k:
                mes[i][j] = A[i][k] * A[k][j]
                B[i][j] += mes[i][j]
    # integrate out mass terms: ONLY a NEW meson can be massive, and only against
    # an OPPOSITE existing arrow (this is what the superpotential Delta-term pairs
    # up).  The old blanket min over every pair wrongly cancelled pre-existing
    # vector-like bifundamentals that have no mass term (e.g. spectator conifold
    # pairs, the C^3/(Z2xZ2) dual).  Never touch the diagonal (an adjoint is not
    # integrated out here, and the blanket loop turned M[ii] into -M[ii]).
    # DWZ (`wmutation.mutate`) remains authoritative when W is tracked; this is
    # the safer adjacency-only fallback for when it is not.
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            c = min(mes[i][j], B[j][i])
            B[i][j] -= c
            B[j][i] -= c
    ranks2 = list(ranks)
    ranks2[k] = rk
    return B, ranks2


class DualState:
    """The state of an interactive Seiberg-duality path: always a quiver
    (adjacency + ranks in units of N); a brane tiling too while every move is
    an urban-renewal square move (all ranks equal, N_f = 2N_c node).  The
    superpotential `W` (list of (coeff, cyclic word)) is tracked through
    general moves by DWZ mutation (`wmutation`); None if a reduction step was
    beyond the implementation."""

    def __init__(self, tiling, adjacency, ranks, reason="", W=None):
        self.tiling = tiling                 # BraneTiling or None
        self.adjacency = adjacency
        self.ranks = ranks
        self.reason = reason                 # why the dimer is unavailable
        self.W = W                           # potential, or None (untracked)


def seiberg_path(vertices, path, start_phase: int = 0) -> DualState:
    """Apply a sequence of Seiberg dualities on ARBITRARY quiver nodes.

    Each move on a node with N_f = 2N_c (a square gauge face, all ranks equal)
    is performed as dimer urban renewal -- the result stays a toric phase with
    a brane tiling.  A move on any OTHER node is genuine Seiberg duality with
    N_c -> N_f - N_c: the ranks become unequal, the theory leaves the toric /
    dimer regime, and only the quiver (with ranks) is tracked from then on.
    Node identity is preserved throughout.  Raises ValueError on an ill-defined
    move (unknown node, dual rank <= 0)."""
    from . import wmutation as WM

    if start_phase:
        phases = enumerate_toric_phases(vertices)
        if not (0 <= start_phase < len(phases)):
            raise ValueError(f"no toric phase #{start_phase}")
        t = _normalize_tiling(phases[start_phase], vertices)
    else:
        t = _normalize_tiling(inverse_quiver(vertices), vertices)
    A = t.adjacency_int()
    ranks = [1] * t.num_gauge
    try:                                     # potential, DWZ-tracked past the
        arrows, W = WM.tiling_potential(t)   # dimer regime
    except WM.WMutationError:
        arrows, W = None, None
    reason = ""
    for f in path:
        if not (0 <= f < len(A)):
            raise ValueError(f"no gauge node {f} in the current quiver")
        if t is not None:
            dimer = t.to_dimer()
            faces, _ = _trace_faces(dimer, _FACE_HAND)
            lab = t.face_to_node or list(range(len(faces)))
            tf = lab.index(f)                # displayed node id -> trace face
            is_square = len(faces[tf]) == 4
            nxt = (urban_renewal(dimer, faces[tf], vertices, node_of_face=lab)
                   if is_square else None)
            if nxt is not None:
                t = nxt
                A = t.adjacency_int()
                try:
                    arrows, W = WM.tiling_potential(t)   # the dimer's exact W
                except WM.WMutationError:
                    arrows, W = None, None
                continue
            # genuine Seiberg duality beyond the dimer regime
            t = None
            reason = ((f"Seiberg duality on node {f} has N_f != 2N_c: the "
                       "dual ranks are unequal, so this is a non-toric phase "
                       "with no brane tiling / dimer description")
                      if not is_square else
                      (f"urban renewal on node {f} did not certify; the "
                       "quiver is tracked, but no dimer is available"))
        A, ranks = quiver_seiberg(A, ranks, f)
        if W is not None:
            try:                             # DWZ mutation with potential
                arrows, W = WM.mutate(arrows, W, f)
                A_w = WM.adjacency_of(arrows, len(A))
                if A_w != A:
                    # the potential-aware reduction is the physical one: some
                    # vector-like pair had no mass term and must NOT cancel
                    A = A_w
            except WM.WMutationError as exc:
                arrows, W = None, None
                reason = (reason + "  Superpotential no longer tracked: "
                          + str(exc)).strip()
    return DualState(t, A, ranks, reason, W=W)


def dualize_path_json(vertices, path, start_phase: int = 0) -> dict:
    """JSON for an interactive Seiberg-duality path (see `seiberg_path`)."""
    try:
        st = seiberg_path(vertices, path, start_phase)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    if st.tiling is not None:
        out = _tiling_json(st.tiling, vertices)
        out["ranks"] = st.ranks
        out["dimer_available"] = True
    else:
        from . import wmutation as WM
        n = len(st.adjacency)
        anomaly = all(
            sum(st.adjacency[i][k] * st.ranks[i] for i in range(n))
            == sum(st.adjacency[k][j] * st.ranks[j] for j in range(n))
            for k in range(n))
        out = {
            "available": True,
            "dimer_available": False,
            "reason": st.reason,
            "num_gauge": n,
            "num_fields": sum(sum(r) for r in st.adjacency),
            "adjacency": st.adjacency,
            "ranks": st.ranks,
            "checks": {"anomaly_free": anomaly},
            "square_faces": [],              # no dimer moves from here
            # the superpotential, tracked through the general dualities by
            # DWZ mutation (None when a reduction was beyond the tracker)
            "superpotential_w": WM.w_json(st.W) if st.W is not None else None,
            "note": st.reason,
        }
    out["dual_path"] = list(path)
    out["start_phase"] = start_phase
    return out


def _tiling_json(t, vertices) -> dict:
    """Serialise one BraneTiling (a single toric phase) for the web API."""
    checks = dict(t.checks)
    if t.num_white <= KAST_NEWTON_MAX_WHITE:
        newton = kasteleyn_newton_polygon(t)
        checks["kasteleyn_newton_matches"] = gl2z_equiv(newton, convex_hull(vertices))
        newton_out = [list(v) for v in newton]
    else:
        newton = None                       # O(n!) certificate skipped (too large)
        checks["kasteleyn_newton_matches"] = None
        newton_out = None
    return {
        "available": True,
        "num_gauge": t.num_gauge,
        "num_fields": t.num_fields,
        "num_white": t.num_white,
        "num_black": t.num_black,
        "adjacency": t.adjacency_int(),
        "fields": t.fields,
        "superpotential": t.superpotential,
        "kasteleyn_newton": newton_out,
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
            # gauge-face boundary polygons in the same coordinates, so the
            # frontend can label faces and make SQUARE faces clickable
            # Seiberg-duality (urban renewal) moves.
            "faces": face_polygons(t),
        },
        "checks": checks,
        # gauge nodes with N_f = 2N_c (square faces): the available
        # Seiberg-duality / urban-renewal moves on THIS tiling.
        "square_faces": square_gauge_faces(t),
        "note": t.note,
    }


def inverse_quiver_json(vertices, **kw) -> dict:
    """JSON-friendly inverse-algorithm result for the web API (or an error).
    The tiling is normalized so gauge-node numbering matches the face indices
    used by the interactive urban-renewal moves (`square_faces`)."""
    try:
        t = _normalize_tiling(inverse_quiver(vertices, **kw), vertices)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    return _tiling_json(t, vertices)


def inverse_phases_json(vertices, max_phases: int = 12) -> dict:
    """All distinct toric (Seiberg-dual) phases of the diagram, for the web API:
    a list of serialised tilings (phase 0 = the seed), so the toric tab can
    cycle between them.  Each is independently Newton-certified and normalized
    so its gauge-node numbering matches the interactive urban-renewal face
    indices (`square_faces`)."""
    try:
        phases = enumerate_toric_phases(vertices, max_phases=max_phases)
    except ValueError as exc:
        return {"available": False, "error": str(exc)}
    return {
        "available": True,
        "num_phases": len(phases),
        "phases": [_tiling_json(_normalize_tiling(t, vertices), vertices)
                   for t in phases],
    }
