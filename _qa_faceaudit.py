"""Audit: (1) is the quiver after a face-click a CORRECT Seiberg dual?
(2) do the clickable face polygons actually sit over the right faces?"""
import json
from conformalmanifold.inverse import (
    inverse_quiver, inverse_phases_json, dualize_path_json,
    _normalize_tiling, square_gauge_faces, face_polygons)

CASES = [
    ("C3",       [(0, 0), (1, 0), (0, 1)]),
    ("conifold", [(0, 0), (1, 0), (1, 1), (0, 1)]),
    ("dP0",      [(1, 0), (0, 1), (-1, -1)]),
    ("F0",       [(-1, 0), (1, 0), (0, 1), (0, -1)]),
    ("dP1",      [(1, 0), (0, 1), (-1, -1), (0, -1)]),
    ("dP2",      [(1, 0), (0, 1), (-1, 0), (-1, -1), (0, -1)]),
    ("dP3",      [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]),
    ("Z2xZ2",    [(0, 0), (2, 0), (0, 2)]),
    ("SPP",      [(0, 0), (2, 0), (1, 1), (0, 1)]),
]

def point_in_poly(px, py, poly):
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % n]
        if (y1 > py) != (y2 > py):
            if px < x1 + (py - y1) * (x2 - x1) / (y2 - y1):
                inside = not inside
    return inside

fails = []

for label, verts in CASES:
    t = _normalize_tiling(inverse_quiver(verts), verts)
    A = t.adjacency_int()
    n = len(A)
    squares = square_gauge_faces(t)

    for k in squares:
        d = dualize_path_json(verts, [k])
        if not d.get("available"):
            fails.append((label, k, "dualize failed", d.get("error")))
            continue
        A2 = d["adjacency"]
        # (a) arrows through the dualized node are exactly reversed
        for j in range(n):
            if j == k:
                continue
            if A2[k][j] != A[j][k] or A2[j][k] != A[k][j]:
                fails.append((label, k, "reversal at dual node broken",
                              {"j": j, "A": [A[k][j], A[j][k]],
                               "A2": [A2[k][j], A2[j][k]]}))
        # (b) chiral flux between spectators: A + meson composition
        for i in range(n):
            for j in range(n):
                if k in (i, j) or i == j:
                    continue
                want = (A[i][j] + A[i][k] * A[k][j]) - (A[j][i] + A[j][k] * A[k][i])
                got = A2[i][j] - A2[j][i]
                if want != got:
                    fails.append((label, k, "spectator chirality broken",
                                  {"i": i, "j": j, "want": want, "got": got}))
        # (c) anomaly freedom of the dual
        for m in range(n):
            if sum(A2[i][m] for i in range(n)) != sum(A2[m][j] for j in range(n)):
                fails.append((label, k, "dual not anomaly-free", {"node": m}))
        # (d) involution: same node twice returns the seed exactly
        d2 = dualize_path_json(verts, [k, k])
        if not d2.get("available") or d2["adjacency"] != A:
            fails.append((label, k, "involution broken",
                          d2.get("adjacency")))

    # (e) label placement: each face's centroid must lie inside ITS polygon
    #     and inside NO other face polygon (over 3x3 lattice translates)
    ph = inverse_phases_json(verts)
    phases = ph["phases"] if ph.get("available") else \
        [json.loads(json.dumps({"tiling": {"faces": face_polygons(t)}, "num_gauge": n}))]
    for pk, p in enumerate(phases):
        fs = p["tiling"]["faces"]
        for f in fs:
            cx = sum(q[0] for q in f["poly"]) / len(f["poly"])
            cy = sum(q[1] for q in f["poly"]) / len(f["poly"])
            hits = []
            for g in fs:
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        gp = [[q[0] + di, q[1] + dj] for q in g["poly"]]
                        if point_in_poly(cx, cy, gp):
                            hits.append(g["node"])
            if hits != [f["node"]]:
                fails.append((label, f"phase{pk}", "label/click placement",
                              {"face": f["node"], "centroid_hits": hits}))

if fails:
    print(f"{len(fails)} FAILURES")
    for f in fails:
        print(" ", f)
else:
    print("all clean: duality correct at every square face; every face label "
          "sits in exactly its own polygon, all cases, all phases")
