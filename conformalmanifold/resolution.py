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


def triangulate(hull_in, active=None):
    """Return (pts, tris): the canonical lattice-point list and a triangulation
    (CCW index-triples into pts) of the toric diagram.

    `active` (optional): iterable of lattice points to USE as triangulation
    vertices -- the *blown-up* exceptional divisors of a partial resolution of
    the singular cone.  The hull corners are always used.  With `active=None`
    (default) every lattice point is used and the triangulation is unimodular
    (the fully resolved phase); with a subset, cells of 2*area > 1 remain --
    the RESIDUAL ORBIFOLD SINGULARITIES of the partial resolution.  `pts` is
    always the full canonical lattice-point list (stable indexing)."""
    hull = convex_hull(hull_in)
    if len(hull) < 3:
        raise ValueError("need a 2-dimensional toric diagram (>=3 corners)")
    pts = lattice_points(hull)
    idx = {p: i for i, p in enumerate(pts)}
    corners = set(hull)
    if active is None:
        use = set(pts)
    else:
        use = {(int(round(x)), int(round(y))) for (x, y) in active}
        use = (use & set(pts)) | corners        # corners are not optional
    tris = []
    c0 = idx[hull[0]]
    for k in range(1, len(hull) - 1):           # fan the convex corners
        tris.append(_ccw((c0, idx[hull[k]], idx[hull[k + 1]]), pts))
    for p in pts:
        if p in corners or p not in use:
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
    return is_valid_subdivision(pts, tris, hull, set(pts))


def is_valid_subdivision(pts, tris, hull, active):
    """True iff `tris` is a valid triangulation of `hull` using EXACTLY the
    vertex set `active` (a partial resolution: cells of 2*area > 1 are the
    residual singularities).  Checks: every triangle CCW with positive lattice
    area, total 2*area equals the hull's, vertices used == active, and
    edge-manifold.  `is_valid_triangulation` is the active == all-points case."""
    hull = convex_hull(hull)
    n = len(pts)
    active_idx = {i for i, p in enumerate(pts) if p in active}
    total = 0
    used = set()
    for t in tris:
        if len(t) != 3 or len(set(t)) != 3 or any(not (0 <= v < n) for v in t):
            return False
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        o = _orient(a, b, c)
        if o <= 0:                      # must be CCW, non-degenerate
            return False
        total += o
        used.update(t)
    if total != normalized_area(hull) or used != active_idx:
        return False
    return all(len(o) <= 2 for o in edge_map(tris).values())


def residual_cells(pts, tris):
    """The non-unimodular cells of a (partial-resolution) triangulation: each is
    a leftover orbifold singularity of the CY cone.  Returns
    [{"tri": t, "area2": n, "vertices": [[x,y],[x,y],[x,y]]}] sorted by cell
    index; empty iff the triangulation is fully resolved (all unimodular)."""
    out = []
    for t, (a, b, c) in enumerate(tris):
        n = _orient(pts[a], pts[b], pts[c])
        if n > 1:
            out.append({"tri": t, "area2": n,
                        "vertices": [list(pts[a]), list(pts[b]), list(pts[c])]})
    return out


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
# blow-up / blow-down of the BASE SURFACE (toric star subdivision / -1-curve)
# ----------------------------------------------------------------------------
#
# These are the genuine geometric operations, defined exactly when the diagram
# is the cone over a compact toric surface S: the polygon has EXACTLY ONE
# interior lattice point O, and the fan of S is the set of rays nu_i from O to
# the boundary lattice points (in cyclic order).
#
#   * BLOW UP a smooth torus-fixed point of S = star-subdivide an adjacent ray
#     pair nu_i, nu_{i+1} with det(nu_i, nu_{i+1}) = +-1 by the new ray
#     nu_i + nu_{i+1}; on the diagram this adds the boundary lattice point
#           W  =  A + B - O            (A = O + nu_i,  B = O + nu_{i+1}),
#     growing 2*area by one and adding one external (p,q) leg / gauge node.
#     Iterating from dP0 = P^2:  dP0 -> dP1 -> dP2 -> dP3 -> weak dPs ...
#     (dP_n has n+3 boundary points).
#   * BLOW DOWN a -1-curve of S = contract a boundary point whose ray is the
#     sum of its neighbours' rays, nu = nu_prev + nu_next (with unimodular
#     det(nu_prev, nu_next)) -- the exact inverse.  dP1 -> dP0;  F0 = P^1xP^1
#     correctly has NO blow-down sites (it is minimal).
#
# Diagrams with no interior lattice point (C^3, the conifold, ...) are not
# cones over a compact surface: there is no divisor to blow up and NO sites
# are offered.  (Free-form diagram editing -- adding or removing arbitrary
# lattice points -- is still available by clicking the dots directly; that is
# a change of the CY, not a blow-up.)


def _boundary_cycle(hull):
    """All boundary lattice points of the polygon, in cyclic (hull) order."""
    out = []
    n = len(hull)
    for i in range(n):
        A, B = hull[i], hull[(i + 1) % n]
        dx, dy = B[0] - A[0], B[1] - A[1]
        g = gcd(abs(dx), abs(dy))
        for t in range(g):
            out.append((A[0] + t * (dx // g), A[1] + t * (dy // g)))
    return out


def _surface_base(points):
    """(hull, O, boundary_cycle, lat, area0) when the diagram is a cone over a
    toric surface (exactly one interior lattice point O), else None."""
    hull = convex_hull({(int(round(x)), int(round(y))) for (x, y) in points})
    if len(hull) < 3:
        return None
    inter = interior_lattice_points(hull)
    if len(inter) != 1:
        return None
    return (hull, inter[0], _boundary_cycle(hull),
            set(lattice_points(hull)), normalized_area(hull))


def surface_blowup_candidates(points):
    """The blow-up sites of the base surface (see the section comment above).

    Returns a list (sorted by new point) of dicts:
        {"edge": [[ax,ay],[bx,by]], "apex": [ox,oy], "new_point": [wx,wy]}
    one per adjacent boundary-lattice-point pair A, B whose rays from the
    unique interior point O are unimodular (a smooth fixed point of the
    surface): `new_point` W = A + B - O is the inserted exceptional ray.
    Empty when the diagram is not a cone over a toric surface."""
    base = _surface_base(points)
    if base is None:
        return []
    hull, O, bpts, lat, area0 = base
    m = len(bpts)
    cands = {}
    for i in range(m):
        A, B = bpts[i], bpts[(i + 1) % m]
        ux, uy = A[0] - O[0], A[1] - O[1]
        vx, vy = B[0] - O[0], B[1] - O[1]
        if abs(ux * vy - uy * vx) != 1:   # singular fixed point: no simple blow-up
            continue
        W = (A[0] + B[0] - O[0], A[1] + B[1] - O[1])
        if W in lat or W in cands:
            continue
        if normalized_area(convex_hull(lat | {W})) != area0 + 1:
            continue                      # safety: must add exactly one triangle
        cands[W] = {"edge": [list(A), list(B)], "apex": list(O),
                    "new_point": [W[0], W[1]]}
    return [cands[w] for w in sorted(cands)]


def surface_blowup(points, new_point):
    """Blow up the base surface: add the exceptional ray `new_point=(x,y)`.

    `new_point` must be one of the sites reported by
    `surface_blowup_candidates(points)` (the star subdivision W = A + B - O of a smooth
    fixed point of the base surface); this is enforced so the UI can only
    trigger genuine blow-ups.  Returns the new diagram point set, sorted --
    the same points plus the new exceptional ray."""
    wx, wy = int(round(new_point[0])), int(round(new_point[1]))
    pset = {(int(round(x)), int(round(y))) for (x, y) in points}
    legal = {tuple(c["new_point"]) for c in surface_blowup_candidates(points)}
    if (wx, wy) not in legal:
        raise ValueError(
            f"({wx},{wy}) is not a blow-up site of this diagram; blow-ups need "
            "a cone over a toric surface (exactly one interior lattice point) "
            "and insert the ray sum of an adjacent smooth pair (W = A + B - O)")
    area_before = normalized_area(convex_hull(pset))
    pset.add((wx, wy))
    new_hull = convex_hull(pset)
    if normalized_area(new_hull) != area_before + 1:
        # A genuine blow-up adds exactly one unimodular triangle.
        raise ValueError(f"blowing up at ({wx},{wy}) does not yield a minimal "
                         "exceptional divisor")
    return sorted(pset)


def surface_blowdown_candidates(points):
    """The blow-down sites of the base surface -- the exact inverse of
    `surface_blowup` (see the section comment above): boundary lattice points W whose
    ray from the unique interior point O is the sum of its neighbours' rays,
    nu_W = nu_prev + nu_next with det(nu_prev, nu_next) = +-1 -- an
    exceptional (-1)-curve of the surface, contractible.  dP1 -> dP0; F0 has
    none (minimal surface).  Empty when the diagram is not a cone over a
    toric surface.

    Returns a list (sorted) of dicts {"corner": [wx,wy], "neighbors": [P, N]}
    (P, N = the boundary lattice points adjacent to W); `corner` is the
    lattice point removed by `surface_blowdown`."""
    base = _surface_base(points)
    if base is None:
        return []
    hull, O, bpts, lat, area0 = base
    m = len(bpts)
    if m < 4:                             # must keep >= 3 boundary points
        return []
    out = {}
    for i in range(m):
        P, W, N = bpts[(i - 1) % m], bpts[i], bpts[(i + 1) % m]
        px, py = P[0] - O[0], P[1] - O[1]
        nx, ny = N[0] - O[0], N[1] - O[1]
        if (W[0] - O[0], W[1] - O[1]) != (px + nx, py + ny):
            continue                      # ray is not the neighbour sum
        if abs(px * ny - py * nx) != 1:
            continue                      # not a -1-curve (surface singular)
        nh = convex_hull(lat - {W})
        if len(nh) < 3 or normalized_area(nh) != area0 - 1:
            continue                      # safety: removes exactly one triangle
        out[W] = {"corner": [W[0], W[1]], "neighbors": [list(P), list(N)]}
    return [out[w] for w in sorted(out)]


def surface_blowdown(points, corner):
    """Blow down a (-1)-curve of the base surface: contract the boundary point
    `corner=(x,y)` (the exact inverse of `surface_blowup`), shrinking the diagram by
    one unimodular triangle (one fewer gauge group / external leg).  `corner`
    must be one of the sites reported by `surface_blowdown_candidates(points)`.
    Returns the new diagram point set, sorted -- all lattice points of the
    shrunk polygon, so a corner of a long-edged diagram (placed with 3 clicks)
    contracts without degenerating."""
    cx, cy = int(round(corner[0])), int(round(corner[1]))
    hull = convex_hull({(int(round(x)), int(round(y))) for (x, y) in points})
    legal = {tuple(c["corner"]) for c in surface_blowdown_candidates(points)}
    if (cx, cy) not in legal:
        raise ValueError(
            f"({cx},{cy}) is not a contractible (-1)-curve of the base "
            "surface; choose a highlighted corner (its ray must be the sum of "
            "its neighbours' rays)")
    lat = set(lattice_points(hull))
    lat.discard((cx, cy))
    return sorted(lat)


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

    internal_edges, external_legs = [], []
    for (i, j), owners in edge_map(tris).items():
        dx, dy = pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]
        g = gcd(abs(dx), abs(dy)) or 1
        px, qy = dy // g, -dx // g                     # +- primitive edge normal
        if len(owners) == 2:
            t1, t2 = owners
            # charge = primitive normal of the shared toric edge (i,j); the
            # tropical junctions make J[t2]-J[t1] parallel to it -- orient to
            # match so the label points along the drawn edge.
            wx, wy = junctions[t2][0] - junctions[t1][0], junctions[t2][1] - junctions[t1][1]
            if px * wx + qy * wy < 0:
                px, qy = -px, -qy
            internal_edges.append({"tris": [t1, t2], "pq": [px, qy],
                                   "mult": g})
        elif len(owners) == 1:
            # boundary edge of the subdivision.  In a partial resolution it can
            # span g > 1 primitive segments (inactive boundary points on it):
            # one external leg of multiplicity g (g coincident 5-branes).
            t = owners[0]
            k = _apex(tris[t], i, j)
            vx, vy = pts[k][0] - pts[i][0], pts[k][1] - pts[i][1]
            if px * vx + qy * vy > 0:                  # orient away from cell
                px, qy = -px, -qy
            mid = [(pts[i][0] + pts[j][0]) / 2.0, (pts[i][1] + pts[j][1]) / 2.0]
            external_legs.append({"junction": t, "pq": [px, qy],
                                  "base": mid, "mult": g})
    return {"junctions": junctions,
            "internal_edges": internal_edges,
            "external_legs": external_legs,
            "junction_kind": junction_kind}
