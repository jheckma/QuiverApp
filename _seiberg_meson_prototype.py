"""DEV PROTOTYPE (not wired into the app) — Seiberg duality / urban renewal.

Status 2026-06-22: the meson construction below produces the CORRECT F0 phase II
*quiver* (4 nodes, 12 fields, anomaly-free, valid two-term/toric superpotential,
genus-1 torus dimer), confirmed via `forward_extract` + `phase_invariant`.  The
OPEN problem is the edge **homology** of the new (meson / reversed) edges: the
naive composite rule below does NOT make the Kasteleyn Newton polygon certify
back to the F0 square, so the dual tiling is not yet Newton-verified.  See
HANDOFF_seiberg_duality.md.

Run:  python _seiberg_meson_prototype.py
"""
from itertools import permutations, product

from conformalmanifold import inverse as inv
from conformalmanifold import toric as T

F0 = [(-1, 0), (1, 0), (0, 1), (0, -1)]


def main():
    t = inv.inverse_quiver(F0)
    d = t.to_dimer()
    H = {i: tuple(e["h"]) for i, e in enumerate(d.edges)}

    add = lambda a, b: (a[0] + b[0], a[1] + b[1])
    neg = lambda a: (-a[0], -a[1])

    # --- the Seiberg/meson move on node0 of F0 (forward_extract frame) ----------
    # node0 edges: incoming X={e3,e4}, outgoing Y={e0,e7}; corners (a-edge pairs):
    #   b1:(e0,e3)  w0:(e3,e7)  b0:(e7,e4)  w1:(e4,e0)
    # mesons m_{Yi Xj} (one per corner): m37=Y7X3@w0, m40=Y0X4@w1,
    #                                     m74=Y7X4@b0, m03=Y0X3@b1
    # reversed fields: x3,x4 (in->out of node0), y0,y7 (out->in)
    # OPEN: homology rule below (composite = sum, reversed = -h) fails the Newton
    # certificate; needs a face-closure / torus-embedding solve.
    hm = {"e1": H[1], "e2": H[2], "e5": H[5], "e6": H[6],
          "x3": neg(H[3]), "x4": neg(H[4]), "y0": neg(H[0]), "y7": neg(H[7]),
          "m37": add(H[7], H[3]), "m40": add(H[0], H[4]),
          "m74": add(H[7], H[4]), "m03": add(H[0], H[3])}
    # vertex (superpotential-term) memberships -- VALIDATED (give phase-II quiver):
    Wsets = {"W0": ["e1", "m37", "e5"], "W1": ["m40", "e2", "e6"],
             "W2": ["m74", "y7", "x4"], "W3": ["m03", "y0", "x3"]}
    Bsets = {"B0": ["e1", "e2", "m74"], "B1": ["m03", "e6", "e5"],
             "B2": ["m37", "y7", "x3"], "B3": ["m40", "y0", "x4"]}

    wn, bn = list(Wsets), list(Bsets)
    elabels = list(hm)
    eid = {l: i for i, l in enumerate(elabels)}
    wv = {l: i for i, w in enumerate(wn) for l in Wsets[w]}
    bv = {l: i for i, b in enumerate(bn) for l in Bsets[b]}
    edges = [{"w": wv[l], "b": bv[l], "h": list(hm[l])} for l in elabels]

    # cyclic order within each 3-field term is not fixed by the membership; search
    # the orderings for the genus-1 (F=4), anomaly-free, toric embedding.
    cyc = lambda lst: [[lst[0]] + list(p) for p in permutations(lst[1:])]
    Wopts = [cyc(Wsets[w]) for w in wn]
    Bopts = [cyc(Bsets[b]) for b in bn]
    hull = T.convex_hull(F0)
    n_genus1 = n_newton = 0
    example = None
    for wsel in product(*Wopts):
        for bsel in product(*Bopts):
            rot_w = [[eid[l] for l in o] for o in wsel]
            rot_b = [[eid[l] for l in o] for o in bsel]
            r = inv.forward_extract(inv.DimerGraph(4, 4, edges, rot_w, rot_b), F0)
            c = r.checks
            if c["euler_V_minus_E_plus_F"] == 0 and c["gauge_eq_2area"] \
                    and c["anomaly_free"] and c["toric_superpotential"]:
                n_genus1 += 1
                if example is None:
                    example = r
                if T.gl2z_equiv(inv.kasteleyn_newton_polygon(r), hull):
                    n_newton += 1
    print(f"genus-1 + anomaly-free + toric orderings : {n_genus1}")
    print(f"...of which Newton-certify to F0 square   : {n_newton}  (OPEN: want >0)")
    if example is not None:
        print(f"example dual quiver: gauge={example.num_gauge} fields={example.num_fields}")
        print(f"  phase_invariant = {inv.phase_invariant(example)}")
        print(f"  (F0 phase I invariant = {inv.phase_invariant(t)})")


if __name__ == "__main__":
    main()
