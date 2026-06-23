"""Toric resolutions: triangulations of the toric diagram and their (p,q) webs.

A crepant resolution of a toric Calabi-Yau three-fold singularity is a
triangulation of its toric diagram into *unimodular* (area-1/2, lattice-point-
free) triangles, using ALL lattice points of the diagram (boundary + interior)
as vertices.  Each such triangulation is a smooth toric phase; different
triangulations are related by **flops** (swapping the diagonal of the quadri-
lateral formed by two adjacent triangles).

The dual of a triangulation is the (p,q) 5-brane web.  It is the **tropical /
Legendre dual** of the triangulation: lift the lattice points to a strictly
convex height function nu whose lower hull projects to exactly this triangulation
(a "regular" triangulation), and send each triangle T to the gradient grad(nu|_T)
of the affine piece of the lift over T.  This junction placement makes the web a
genuine (p,q) web:

    * one trivalent **junction** per triangle, at grad(nu|_T);
    * one finite internal **edge** per internal triangulation edge, joining the
      two adjacent junctions.  Because the two affine pieces agree along the
      shared toric edge, the segment between their gradients is **exactly
      perpendicular** to that toric edge and carries (p,q) charge equal to the
      edge normal -- so the legs meeting at every junction balance to zero (local
      5-brane charge conservation), and a flop visibly restructures the web
      (e.g. the resolved conifold's internal edge has positive length and rotates
      by 90 degrees under its flop, instead of collapsing to a point);
    * one semi-infinite external **leg** per boundary segment, leaving the
      adjacent junction along the outward edge-normal with **exact** (p,q) charge
      equal to that normal.

Contrast with two simpler-but-wrong junction choices kept here only for
reference: the **centroid** is never degenerate but its internal edges are not
perpendicular to the toric edges (junctions do not balance); the **circumcenter**
gives the metrically perpendicular web for the *Delaunay* phase but collapses to a
point for cocircular quads such as the conifold square, hiding the flop.  The
Legendre dual fixes both.

Heights are computed by Motzkin relaxation (floats); junction gradients are exact
over the CCW unimodular triangles (determinant 1).  Pure stdlib (no numpy needed).
"""

from __future__ import annotations

from math import gcd

from .toric import convex_hull, normalized_area, polygon_signature


# ----------------------------------------------------------------------------
# exact integer predicates
# ----------------------------------------------------------------------------
def _orient(a, b, c):
    """>0 iff a,b,c make a left turn (CCW); 0 iff collinear."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _in_circumcircle(a, b, c, d):
    """For CCW triangle a,b,c: >0 iff d is strictly inside its circumcircle."""
    ax, ay = a[0] - d[0], a[1] - d[1]
    bx, by = b[0] - d[0], b[1] - d[1]
    cx, cy = c[0] - d[0], c[1] - d[1]
    return ((ax * ax + ay * ay) * (bx * cy - cx * by)
            - (bx * bx + by * by) * (ax * cy - cx * ay)
            + (cx * cx + cy * cy) * (ax * by - bx * ay))


def point_in_convex(p, hull):
    """True iff p is inside or on the boundary of the CCW convex polygon `hull`."""
    n = len(hull)
    for i in range(n):
        if _orient(hull[i], hull[(i + 1) % n], p) < 0:
            return False
    return True


# ----------------------------------------------------------------------------
# lattice points of the toric diagram
# ----------------------------------------------------------------------------
def lattice_points(hull):
    """All lattice points (boundary + interior) of the polygon, in canonical
    (sorted) order -- this order indexes triangulation vertices."""
    hull = convex_hull(hull)
    xs = [p[0] for p in hull]
    ys = [p[1] for p in hull]
    pts = []
    for x in range(min(xs), max(xs) + 1):
        for y in range(min(ys), max(ys) + 1):
            if point_in_convex((x, y), hull):
                pts.append((x, y))
    return sorted(pts)


def interior_lattice_points(hull):
    """Lattice points strictly inside the polygon (no boundary)."""
    hull = convex_hull(hull)
    n = len(hull)
    out = []
    for p in lattice_points(hull):
        if all(_orient(hull[i], hull[(i + 1) % n], p) > 0 for i in range(n)):
            out.append(p)
    return out


def boundary_lattice_points_ordered(hull):
    """Boundary lattice points (a point is on the boundary iff it lies on some
    edge).  Returned sorted (canonical)."""
    hull = convex_hull(hull)
    inter = set(interior_lattice_points(hull))
    return [p for p in lattice_points(hull) if p not in inter]


# ----------------------------------------------------------------------------
# triangulation: fan + point insertion (always valid) + Lawson Delaunay flips
# ----------------------------------------------------------------------------
def _ccw(tri, pts):
    a, b, c = (pts[i] for i in tri)
    return tri if _orient(a, b, c) > 0 else (tri[0], tri[2], tri[1])


def _insert_point(tris, pts, pi):
    """Insert vertex index pi into the current triangulation (split the triangle
    that contains it, or the triangles sharing the edge it lies on)."""
    P = pts[pi]
    inside = None
    edge_tris = []
    for ti, (a, b, c) in enumerate(tris):
        A, B, C = pts[a], pts[b], pts[c]
        o1, o2, o3 = _orient(A, B, P), _orient(B, C, P), _orient(C, A, P)
        if o1 > 0 and o2 > 0 and o3 > 0:
            inside = ti
            break
        if o1 >= 0 and o2 >= 0 and o3 >= 0:
            zeros = (o1 == 0, o2 == 0, o3 == 0)
            if sum(zeros) == 1:
                edge_tris.append((ti, zeros))
    if inside is not None:
        a, b, c = tris[inside]
        tris[inside] = (a, b, pi)
        tris.append((b, c, pi))
        tris.append((c, a, pi))
        return
    new, remove = [], set()
    for ti, zeros in edge_tris:
        a, b, c = tris[ti]
        remove.add(ti)
        if zeros[0]:          # P on edge AB, apex C
            new += [(a, pi, c), (pi, b, c)]
        elif zeros[1]:        # P on edge BC, apex A
            new += [(b, pi, a), (pi, c, a)]
        else:                 # P on edge CA, apex B
            new += [(c, pi, b), (pi, a, b)]
    tris[:] = [t for i, t in enumerate(tris) if i not in remove] + new


def _delaunay_flips(tris, pts, max_iter=10000):
    """Lawson edge flips toward the Delaunay triangulation (nicer default web).
    Only strictly-improving flips on convex quads are taken, so it terminates;
    any outcome is still a valid unimodular triangulation."""
    for _ in range(max_iter):
        em = edge_map(tris)
        flipped = False
        for edge, owners in em.items():
            if len(owners) != 2:
                continue
            i, j = edge
            t1, t2 = owners
            k = _apex(tris[t1], i, j)
            l = _apex(tris[t2], i, j)
            if k is None or l is None:
                continue
            A, Bp, C, D = pts[i], pts[j], pts[k], pts[l]
            tri = _ccw((i, j, k), pts)
            a, b, c = (pts[tri[0]], pts[tri[1]], pts[tri[2]])
            if _in_circumcircle(a, b, c, D) > 0 and _convex_quad(pts, i, j, k, l):
                tris[t1] = _ccw((i, k, l), pts)
                tris[t2] = _ccw((j, k, l), pts)
                flipped = True
                break
        if not flipped:
            return


def triangulate(hull_in):
    """Return (pts, tris): the canonical lattice-point list and a unimodular
    triangulation (CCW index-triples into pts) of the toric diagram."""
    hull = convex_hull(hull_in)
    if len(hull) < 3:
        raise ValueError("need a 2-dimensional toric diagram (>=3 corners)")
    pts = lattice_points(hull)
    idx = {p: i for i, p in enumerate(pts)}
    tris = []
    c0 = idx[hull[0]]
    for k in range(1, len(hull) - 1):           # fan the convex corners
        tris.append(_ccw((c0, idx[hull[k]], idx[hull[k + 1]]), pts))
    corners = set(hull)
    for p in pts:
        if p in corners:
            continue
        _insert_point(tris, pts, idx[p])
    _delaunay_flips(tris, pts)
    return pts, [tuple(t) for t in tris]


# ----------------------------------------------------------------------------
# edges, apexes, flops
# ----------------------------------------------------------------------------
def is_valid_triangulation(pts, tris, hull):
    """True iff `tris` (index-triples into `pts`) is a genuine unimodular
    triangulation of `hull`: right number of triangles (= 2*area), each CCW and
    lattice-area-1/2, every lattice point used, and edge-manifold (no edge shared
    by >2 triangles).  Guards the API against a hand-crafted bogus `tri=` query."""
    hull = convex_hull(hull)
    n = len(pts)
    if len(tris) != normalized_area(hull):
        return False
    used = set()
    for t in tris:
        if len(t) != 3 or len(set(t)) != 3 or any(not (0 <= v < n) for v in t):
            return False
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        o = _orient(a, b, c)
        if o <= 0 or o != 1:           # CCW and unimodular (2*area == 1)
            return False
        used.update(t)
    if used != set(range(n)):
        return False
    return all(len(o) <= 2 for o in edge_map(tris).values())


def edge_map(tris):
    """{(i,j) sorted: [triangle indices owning that edge]}."""
    em = {}
    for ti, (a, b, c) in enumerate(tris):
        for (u, v) in ((a, b), (b, c), (c, a)):
            em.setdefault((min(u, v), max(u, v)), []).append(ti)
    return em


def _apex(tri, i, j):
    """The vertex of triangle `tri` other than i and j (or None)."""
    rest = [v for v in tri if v != i and v != j]
    return rest[0] if len(rest) == 1 else None


def _convex_quad(pts, i, j, k, l):
    """True iff the quad with diagonal (i,j) and apexes k,l is strictly convex,
    i.e. the flop to diagonal (k,l) is geometrically valid."""
    # k,l must lie on opposite sides of line ij, and i,j on opposite sides of kl
    s1 = _orient(pts[i], pts[j], pts[k])
    s2 = _orient(pts[i], pts[j], pts[l])
    s3 = _orient(pts[k], pts[l], pts[i])
    s4 = _orient(pts[k], pts[l], pts[j])
    return (s1 * s2 < 0) and (s3 * s4 < 0)


def flippable_edges(pts, tris):
    """Internal edges whose two triangles form a convex quad (can be flopped)."""
    out = []
    for (i, j), owners in edge_map(tris).items():
        if len(owners) != 2:
            continue
        k = _apex(tris[owners[0]], i, j)
        l = _apex(tris[owners[1]], i, j)
        if k is not None and l is not None and _convex_quad(pts, i, j, k, l):
            out.append([i, j])
    return out


def flop(pts, tris, edge):
    """Flop the internal edge `edge=(i,j)`: replace diagonal (i,j) with (k,l).
    Returns a new triangulation list.  Raises ValueError if not flippable."""
    i, j = (min(edge), max(edge))
    owners = edge_map(tris).get((i, j), [])
    if len(owners) != 2:
        raise ValueError(f"edge {(i, j)} is not an internal edge")
    k = _apex(tris[owners[0]], i, j)
    l = _apex(tris[owners[1]], i, j)
    if k is None or l is None or not _convex_quad(pts, i, j, k, l):
        raise ValueError(f"edge {(i, j)} is not flippable (non-convex quad)")
    new = []
    for ti, t in enumerate(tris):
        if ti in owners:
            continue
        new.append(t)
    new.append(_ccw((i, k, l), pts))
    new.append(_ccw((j, k, l), pts))
    return new


# ----------------------------------------------------------------------------
# blow-up: chamfer a (singular) corner of the toric diagram
# ----------------------------------------------------------------------------
def blowup_corner(points, corner):
    """Blow up the toric diagram at the convex-hull corner `corner=(x,y)`.

    Geometrically this *chamfers* the corner: it cuts the corner vertex off with
    a new boundary edge, stepping one primitive lattice unit inward along each of
    the two edges meeting at the corner.  Physically it (partially) resolves the
    local singularity sitting at that corner -- the C^2/Z_k x C corner is replaced
    by an exceptional edge -- which changes the diagram, its (p,q) web, the gauge-
    group count (= 2 x area) and the whole reconstructed quiver.

    A *smooth* (unimodular) corner -- both adjacent edges of lattice length 1, as
    in C^3 or the conifold square -- has nothing to resolve; blowing it up raises
    ValueError so the UI can mark only genuinely singular corners as actionable.

    `points` : the diagram's lattice points (any of the user's placed points).
    `corner` : an (x, y) that must be a convex-hull vertex of `points`.
    Returns the new diagram point set, sorted -- the corner removed, the two
    chamfer points added -- ready to feed back through the pipeline."""
    cx, cy = int(round(corner[0])), int(round(corner[1]))
    pset = {(int(round(x)), int(round(y))) for (x, y) in points}
    hull = convex_hull(pset)
    if len(hull) < 3:
        raise ValueError("need a 2-dimensional toric diagram (>=3 corners) to "
                         "blow up a corner")
    if (cx, cy) not in hull:
        raise ValueError(f"({cx},{cy}) is not a corner (convex-hull vertex) of "
                         "the toric diagram")
    n = len(hull)
    i = hull.index((cx, cy))
    V = hull[i]
    P = hull[(i - 1) % n]                 # CCW-previous neighbouring corner
    Nx = hull[(i + 1) % n]                # CCW-next neighbouring corner

    def step(frm, to):
        """primitive lattice step from `frm` toward `to`, and the edge's length."""
        dx, dy = to[0] - frm[0], to[1] - frm[1]
        g = gcd(abs(dx), abs(dy)) or 1
        return (dx // g, dy // g), g

    aP, gP = step(V, P)
    aN, gN = step(V, Nx)
    if gP == 1 and gN == 1:
        raise ValueError(f"corner ({cx},{cy}) is smooth (unimodular) -- there is "
                         "no singularity there to blow up")
    w1 = (V[0] + aP[0], V[1] + aP[1])     # one step toward P (== P iff gP == 1)
    w2 = (V[0] + aN[0], V[1] + aN[1])     # one step toward N (== N iff gN == 1)
    pset.discard(V)
    pset.add(w1)
    pset.add(w2)
    new_hull = convex_hull(pset)
    if len(new_hull) < 3 or normalized_area(new_hull) == 0:
        raise ValueError(f"blowing up corner ({cx},{cy}) degenerates the diagram")
    return sorted(pset)


# ----------------------------------------------------------------------------
# the dual (p,q) web
# ----------------------------------------------------------------------------
def circumcenter(a, b, c):
    """Circumcenter of triangle a,b,c as a float (x,y) -- the metrically exact
    (perpendicular) web junction, kept for reference (see `dual_web`)."""
    ax, ay = a
    bx, by = b
    cx, cy = c
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    a2 = ax * ax + ay * ay
    b2 = bx * bx + by * by
    c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return [ux, uy]


def centroid(a, b, c):
    """Centroid of triangle a,b,c -- a never-degenerate but non-perpendicular
    junction (see module docstring); kept as a drawing fallback."""
    return [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0]


# ----------------------------------------------------------------------------
# tropical / Legendre dual: strictly-convex lift -> perpendicular (p,q) web
# ----------------------------------------------------------------------------
def _affine_coeffs(pts, i, j, k, l):
    """Barycentric coords (a,b,c) of point l w.r.t. the affine frame (i,j,k):
    l == a*i + b*j + c*k with a+b+c == 1.  c<0 iff l is across edge (i,j) from k
    (the configuration of two triangles sharing an internal edge)."""
    xi, yi = pts[i]
    xj, yj = pts[j]
    xk, yk = pts[k]
    xl, yl = pts[l]
    det = xi * (yj - yk) - xj * (yi - yk) + xk * (yi - yj)
    a = (xl * (yj - yk) - xj * (yl - yk) + xk * (yl - yj)) / det
    b = (xi * (yl - yk) - xl * (yi - yk) + xk * (yi - yl)) / det
    return a, b, 1.0 - a - b


def convex_heights(pts, tris, margin=1.0, max_iter=20000):
    """A strictly-convex height nu (list, one float per lattice point) whose lower
    hull projects to exactly the triangulation `tris`.  For every internal edge
    (i,j) with neighbour apexes k,l the lift of l must sit strictly above the
    affine plane through the lifted (i,j,k) -- a strict linear inequality on nu.
    Solve the system by Motzkin relaxation (projection onto each violated
    constraint); the paraboloid lift seeds it near the Delaunay phase.  Returns
    None if it fails to converge (a non-regular triangulation -- not produced by
    flopping the default phase of the small toric diagrams handled here)."""
    h = [float(x * x + y * y) for (x, y) in pts]      # paraboloid seed
    cons = []
    for (i, j), owners in edge_map(tris).items():
        if len(owners) != 2:
            continue
        k = _apex(tris[owners[0]], i, j)
        l = _apex(tris[owners[1]], i, j)
        if k is None or l is None:
            continue
        a, b, c = _affine_coeffs(pts, i, j, k, l)     # value at l = a*hi+b*hj+c*hk
        coeff = {}                                    # want h_l - (a hi+b hj+c hk) >= margin
        for v, co in ((l, 1.0), (i, -a), (j, -b), (k, -c)):
            coeff[v] = coeff.get(v, 0.0) + co
        cons.append(coeff)
    for _ in range(max_iter):
        updated = False
        for coeff in cons:
            r = sum(co * h[v] for v, co in coeff.items())
            if r < margin - 1e-9:
                nn = sum(co * co for co in coeff.values())
                step = (margin - r) / nn
                for v, co in coeff.items():
                    h[v] += step * co
                updated = True
        if not updated:
            return h
    return None


def _triangle_gradient(pts, h, tri):
    """Gradient (mx,my) of the affine interpolant of heights `h` over CCW
    unimodular triangle `tri` -- the (p,q)-web junction for that triangle.
    The frame determinant is +1 (CCW, area 1/2), so this is exact up to `h`."""
    a, b, c = tri
    x0, y0 = pts[a]
    x1, y1 = pts[b]
    x2, y2 = pts[c]
    dx1, dy1, dx2, dy2 = x1 - x0, y1 - y0, x2 - x0, y2 - y0
    det = dx1 * dy2 - dy1 * dx2                        # == 1 for CCW unimodular
    z1, z2 = h[b] - h[a], h[c] - h[a]
    mx = (z1 * dy2 - z2 * dy1) / det
    my = (dx1 * z2 - dx2 * z1) / det
    return [mx, my]


def tropical_junctions(pts, tris):
    """Legendre-dual junction positions (one [x,y] per triangle), or None if no
    strictly-convex lift induces this triangulation."""
    h = convex_heights(pts, tris)
    if h is None:
        return None
    return [_triangle_gradient(pts, h, t) for t in tris]


def dual_web(pts, tris, hull):
    """The (p,q) web dual to triangulation `tris` of polygon `hull`.

    Returns {"junctions": [...], "internal_edges": [{"tris":[t1,t2], "pq":[p,q]}],
             "external_legs": [{"junction": t, "pq": [p,q], "base": [x,y]}, ...],
             "junction_kind": "tropical"|"centroid"}.
    `junctions[t]` is the tropical/Legendre-dual junction of triangle t (see the
    module docstring); it falls back to the centroid only if no convex lift is
    found, in which case `junction_kind == "centroid"` and the web is schematic
    (internal edges no longer perpendicular)."""
    hull = convex_hull(hull)
    junctions = tropical_junctions(pts, tris)
    junction_kind = "tropical"
    if junctions is None:                              # non-regular fallback
        junctions = [centroid(pts[a], pts[b], pts[c]) for (a, b, c) in tris]
        junction_kind = "centroid"
    # which (i,j) lattice segments are on the polygon boundary?
    bset = set()
    n = len(hull)
    for e in range(n):
        a, b = hull[e], hull[(e + 1) % n]
        dx, dy = b[0] - a[0], b[1] - a[1]
        g = gcd(abs(dx), abs(dy)) or 1
        ux, uy = dx // g, dy // g
        for t in range(g):
            p0 = (a[0] + ux * t, a[1] + uy * t)
            p1 = (a[0] + ux * (t + 1), a[1] + uy * (t + 1))
            bset.add((p0, p1, (uy, -ux)))   # segment + outward normal (CCW)
    idx = {p: i for i, p in enumerate(pts)}
    bnorm = {}
    for (p0, p1, nrm) in bset:
        key = (min(idx[p0], idx[p1]), max(idx[p0], idx[p1]))
        bnorm[key] = nrm

    internal_edges, external_legs = [], []
    for (i, j), owners in edge_map(tris).items():
        if len(owners) == 2:
            t1, t2 = owners
            # charge = primitive normal of the shared toric edge (i,j); the
            # tropical junctions make J[t2]-J[t1] parallel to it -- orient to
            # match so the label points along the drawn edge.
            dx, dy = pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]
            g = gcd(abs(dx), abs(dy)) or 1
            px, qy = dy // g, -dx // g                 # primitive edge normal
            wx, wy = junctions[t2][0] - junctions[t1][0], junctions[t2][1] - junctions[t1][1]
            if px * wx + qy * wy < 0:
                px, qy = -px, -qy
            internal_edges.append({"tris": [t1, t2], "pq": [px, qy]})
        elif len(owners) == 1 and (i, j) in bnorm:
            t = owners[0]
            mid = [(pts[i][0] + pts[j][0]) / 2.0, (pts[i][1] + pts[j][1]) / 2.0]
            external_legs.append({"junction": t, "pq": list(bnorm[(i, j)]),
                                  "base": mid})
    return {"junctions": junctions,
            "internal_edges": internal_edges,
            "external_legs": external_legs,
            "junction_kind": junction_kind}
