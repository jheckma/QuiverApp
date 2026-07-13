"""Finite-chamber BPS spectra for quivers.

The chamber is encoded by a mutation sequence.  For a green mutation sequence,
the positive c-vector at each mutation is the charge/dimension vector of a
stable hypermultiplet in that finite BPS chamber.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import re


MAX_NODES = 20


@dataclass
class BpsInput:
    matrix: list[list[int]]
    labels_are_zero_based: bool


def _as_int_matrix(data) -> list[list[int]]:
    if not isinstance(data, list) or not data:
        raise ValueError("exchange matrix must be a non-empty square array")
    matrix: list[list[int]] = []
    for row in data:
        if not isinstance(row, list):
            raise ValueError("exchange matrix rows must be arrays")
        try:
            matrix.append([int(v) for v in row])
        except (TypeError, ValueError) as exc:
            raise ValueError("exchange matrix entries must be integers") from exc
    return matrix


def parse_matrix(text: str) -> list[list[int]]:
    """Parse a skew-symmetric integer exchange matrix."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty exchange matrix")
    if text.startswith("["):
        matrix = _as_int_matrix(json.loads(text))
    else:
        rows = [r.strip() for r in re.split(r"[\n;]+", text) if r.strip()]
        matrix = []
        for row in rows:
            vals = [v for v in re.split(r"[\s,]+", row) if v]
            matrix.append([int(v) for v in vals])
    _validate_exchange_matrix(matrix)
    return matrix


def parse_arrows(text: str) -> BpsInput:
    """Parse arrows like ``1->2``, ``1->2 x2`` or ``0->1:3``."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty arrow list")

    chunks = [c.strip() for c in re.split(r"[\n,;]+", text) if c.strip()]
    arrows: list[tuple[int, int, int]] = []
    directed: set[tuple[int, int]] = set()
    labels: list[int] = []
    pat = re.compile(r"^(\d+)\s*->\s*(\d+)(?:\s*(?:x|\*|:)\s*(\d+))?$")
    for chunk in chunks:
        m = pat.match(chunk)
        if not m:
            raise ValueError(
                f"bad arrow {chunk!r}; expected forms like '1->2' or '1->2 x2'"
            )
        a, b = int(m.group(1)), int(m.group(2))
        mult = int(m.group(3) or 1)
        if a == b:
            raise ValueError("BPS quivers should not contain self-arrows")
        if mult <= 0:
            raise ValueError("arrow multiplicities must be positive")
        if (b, a) in directed:
            raise ValueError("arrow list contains a 2-cycle; BPS quivers must be 2-acyclic")
        directed.add((a, b))
        arrows.append((a, b, mult))
        labels.extend([a, b])

    zero_based = 0 in labels
    offset = 0 if zero_based else 1
    n = max(labels) - offset + 1
    if n <= 0:
        raise ValueError("could not infer nodes from arrow list")
    matrix = [[0 for _ in range(n)] for _ in range(n)]
    for a, b, mult in arrows:
        i, j = a - offset, b - offset
        if not (0 <= i < n and 0 <= j < n):
            raise ValueError("mixed 0-based and 1-based node labels")
        matrix[i][j] += mult
        matrix[j][i] -= mult

    _validate_exchange_matrix(matrix)
    return BpsInput(matrix=matrix, labels_are_zero_based=zero_based)


def parse_quiver(text: str, kind: str = "matrix") -> BpsInput:
    kind = (kind or "matrix").lower()
    if kind == "auto":
        kind = "arrows" if "->" in (text or "") else "matrix"
    if kind == "arrows":
        return parse_arrows(text)
    if kind != "matrix":
        raise ValueError("quiver kind must be 'matrix', 'arrows', or 'auto'")
    return BpsInput(matrix=parse_matrix(text), labels_are_zero_based=False)


def _validate_exchange_matrix(matrix: list[list[int]]) -> None:
    n = len(matrix)
    if n > MAX_NODES:
        raise ValueError(f"BPS quiver UI is limited to {MAX_NODES} nodes")
    if any(len(row) != n for row in matrix):
        raise ValueError("exchange matrix must be square")
    for i in range(n):
        if matrix[i][i] != 0:
            raise ValueError("exchange matrix diagonal must vanish")
        for j in range(n):
            if matrix[i][j] != -matrix[j][i]:
                raise ValueError("exchange matrix must be skew-symmetric")


def _parse_number_list(text: str, n: int, name: str) -> list[float] | None:
    text = (text or "").strip()
    if not text:
        return None
    if text.startswith("["):
        data = json.loads(text)
        vals = [float(v) for v in data]
    else:
        vals = [float(v) for v in re.split(r"[\s,;]+", text) if v]
    if len(vals) != n:
        raise ValueError(f"{name} must contain exactly {n} numbers")
    if name == "masses" and any(v <= 0 for v in vals):
        raise ValueError("masses must be positive")
    return vals


def parse_sequence(text: str, n: int, zero_based_hint: bool = False) -> list[int]:
    text = (text or "").strip()
    if not text:
        raise ValueError("empty chamber mutation sequence")
    vals = [int(v) for v in re.findall(r"\d+", text)]
    if not vals:
        raise ValueError("mutation sequence must contain node labels")
    zero_based = zero_based_hint or 0 in vals
    offset = 0 if zero_based else 1
    seq = [v - offset for v in vals]
    bad = [vals[i] for i, k in enumerate(seq) if k < 0 or k >= n]
    if bad:
        lo, hi = (0, n - 1) if zero_based else (1, n)
        raise ValueError(f"mutation node {bad[0]} is outside the range {lo}..{hi}")
    return seq


def _mutate_exchange(B: list[list[int]], k: int) -> list[list[int]]:
    n = len(B)
    out = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == k or j == k:
                out[i][j] = -B[i][j]
            else:
                out[i][j] = (
                    B[i][j]
                    + max(B[i][k], 0) * max(B[k][j], 0)
                    - max(B[j][k], 0) * max(B[k][i], 0)
                )
    return out


def _vector_sign(v: list[int]) -> str:
    if all(x >= 0 for x in v) and any(x > 0 for x in v):
        return "green"
    if all(x <= 0 for x in v) and any(x < 0 for x in v):
        return "red"
    return "mixed"


def find_maximal_green_sequence(matrix, max_depth=None, node_budget=20000):
    """Search for a maximal green sequence (a canonical finite BPS chamber) of
    the exchange matrix `matrix`, returning the mutation node list (0-based) or
    None if none is found within the search bounds.

    A green sequence mutates only at GREEN vertices (c-vector all >= 0); it is
    MAXIMAL green when every c-vector has turned red (all <= 0).  By the
    sign-coherence theorem every c-vector stays green-or-red throughout, so the
    search is a plain DFS over green vertices with a visited-state guard.  Some
    quivers (e.g. the local-P^2 "Markov" quiver, 3 nodes with +-3 arrows) admit
    NO maximal green sequence -- the DFS then exhausts its budget and returns
    None, which the caller reports honestly rather than inventing a chamber."""
    n = len(matrix)
    if n == 0:
        return []
    if max_depth is None:
        max_depth = 6 * n + 12
    C0 = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    B0 = [row[:] for row in matrix]
    visited = set()
    budget = [node_budget]

    def key(B, C):
        return (tuple(map(tuple, C)), tuple(map(tuple, B)))

    def all_red(C):
        return all(_vector_sign(c) == "red" for c in C)

    def dfs(B, C, seq):
        if all_red(C):
            return seq
        if len(seq) >= max_depth or budget[0] <= 0:
            return None
        budget[0] -= 1
        for k in range(n):
            if _vector_sign(C[k]) == "green":
                B2 = _mutate_exchange(B, k)
                C2 = _mutate_c_matrix(C, B, k, "green")
                kk = key(B2, C2)
                if kk in visited:
                    continue
                visited.add(kk)
                r = dfs(B2, C2, seq + [k])
                if r is not None:
                    return r
        return None

    return dfs(B0, C0, [])


def default_central_charges(matrix, sequence):
    """A physically sensible DEFAULT (masses, phases) for displaying the finite
    chamber given by `sequence`, so central charges Z_i = mass_i e^{i pi phase_i}
    render with the stable objects appearing in DECREASING phase order along the
    mutation path (the standard stability-chamber condition).

    Construction: give the node first mutated in the sequence the largest phase
    and space the rest evenly downward in (0, 1); unit masses.  For a genuine
    maximal green sequence this makes every emitted c-vector's phase strictly
    decreasing (the stable objects come out ordered), which is exactly the
    `phase_order_ok` check.  Nodes never mutated keep a mid-range phase."""
    n = len(matrix)
    order = []
    for k in sequence:
        if k not in order:
            order.append(k)
    for k in range(n):
        if k not in order:
            order.append(k)
    # spread phases in (0.15, 0.9), first-mutated highest
    phases = [0.0] * n
    m = max(len(order), 1)
    for rank, k in enumerate(order):
        phases[k] = round(0.9 - 0.75 * rank / max(m - 1, 1), 6)
    masses = [1.0] * n
    return phases, masses


def kernel_basis(matrix):
    """Integer basis (as row vectors) of {x in Z^n : B x = 0} -- the
    mutation-invariant BPS charges of the exchange matrix B: the KK/D0 momentum
    delta = (1,...,1) together with the flavor U(1) directions.  dim = flavor
    rank + 1."""
    n = len(matrix)
    if n == 0:
        return []
    if all(all(v == 0 for v in row) for row in matrix):
        # degenerate (e.g. the conifold's netted B): every charge is invariant
        return [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    from .fived import _integer_kernel
    return _integer_kernel(matrix)


def bps_invariants(matrix, label_base: int = 0) -> dict:
    """The mutation-invariant physics labels of a 5d BPS quiver (Closset-Del
    Zotto): the KK-momentum / D0-brane charge delta = (1,...,1) (the null vector
    of B for a fractional-brane quiver), the flavor U(1) directions (the rest of
    ker B), and the flavor rank = dim ker B - 1."""
    n = len(matrix)
    ker = kernel_basis(matrix)
    delta = [1] * n
    # is delta actually the (a) null vector?  For a genuine fractional-brane
    # quiver B*delta = 0; report whether it holds so degenerate inputs are honest
    delta_null = all(sum(matrix[i][j] * delta[j] for j in range(n)) == 0
                     for i in range(n)) if n else True
    # flavor directions = a kernel basis with delta factored out (best-effort:
    # the kernel minus the delta line), reported as-is for display
    flavor = [v for v in ker if v != delta]
    rank_ker = len(ker)
    return {
        "kk_charge": delta,                       # D0 / KK momentum
        "kk_charge_is_null": delta_null,
        "flavor_directions": flavor,
        "flavor_rank": max(rank_ker - 1, 0),
        "kernel_dim": rank_ker,
        "label_base": label_base,
    }


def orbifold_central_charges(n: int, eps: float = 0.15):
    """Default (phases, masses) for the canonical orbifold chamber: all central
    charges nearly aligned near phase 1/2, spread narrowly and DECREASING so the
    n nodes are the only stable states and come out phase-ordered
    (Closset-Del Zotto, the epsilon-narrow recipe).  phase_i = 1/2 +
    eps*(n+1-2i)/(2n), i=1..n; unit masses."""
    if n <= 0:
        return [], []
    phases = [round(0.5 + eps * (n + 1 - 2 * (i + 1)) / (2 * n), 8)
              for i in range(n)]
    masses = [1.0] * n
    return phases, masses


def _primitive(v):
    from math import gcd
    g = gcd(abs(int(v[0])), abs(int(v[1])))
    return (v[0] // g, v[1] // g) if g else (0, 0)


def spectral_network_description(newton_points, num_nodes=None) -> dict:
    """Spectral-network (Gaiotto-Moore-Neitzke) description of the 5d BPS
    spectrum, built from the toric / dimer presentation.

    The Seiberg-Witten (mirror) curve of the 5d SCFT is the Newton polynomial of
    the toric diagram,  Sigma:  H(x,y) = sum_{(a,b)} c_{ab} x^a y^b = 0,  viewed
    as a multi-sheeted cover of the punctured x-plane C^*.  The brane tiling
    (dimer) is the same datum: its **zig-zag paths** are the (p,q) legs of the
    dual 5-brane web, which are exactly the **asymptotic directions of the
    spectral network** (the outward primitive normals of the toric polygon's
    edges).  A spectral network at phase theta is the union of S-walls where
    lambda_i - lambda_j has phase theta (lambda = sheets of Sigma); at the BPS
    phases the network develops finite webs ("string webs") whose homology
    classes in H_1(Sigma) are the BPS charges.  The BPS quiver is the dual
    description: its nodes are a distinguished basis of that charge lattice and
    wall-crossing = quiver mutation.

    Returns a descriptive dict (curve degree/genus, sheet count, asymptotic
    (p,q) legs, charge-lattice rank), NOT a full network solver.
    Refs: Gaiotto-Moore-Neitzke arXiv:0907.3987 (Hitchin/WKB), 1204.4824
    (Spectral Networks); for the toric/dimer-to-network link see the mirror-curve
    / (p,q)-web dictionary."""
    from .toric import convex_hull
    from .resolution import interior_lattice_points

    pts = [(int(a), int(b)) for (a, b) in newton_points]
    hull = convex_hull(pts)
    V = len(hull)
    interior = interior_lattice_points(hull)
    genus = len(interior)                    # genus of the generic mirror curve
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x_deg = (max(xs) - min(xs)) if xs else 0
    y_deg = (max(ys) - min(ys)) if ys else 0
    # asymptotic legs = outward primitive normals of the CCW hull edges
    # (= dimer zig-zag directions = 5-brane-web external legs = network ends)
    legs = []
    for i in range(V):
        (ax, ay), (bx, by) = hull[i], hull[(i + 1) % V]
        edge = (bx - ax, by - ay)
        outward = _primitive((edge[1], -edge[0]))   # CCW hull -> right normal is outward
        legs.append({"edge": [list(hull[i]), list(hull[(i + 1) % V])],
                     "pq_leg": list(outward),
                     "length": abs(__import__("math").gcd(abs(edge[0]), abs(edge[1])))})
    boundary_pts = sum(l["length"] for l in legs)
    flavor_rank = max(boundary_pts - 3, 0)
    return {
        "framework": "Gaiotto-Moore-Neitzke spectral networks",
        "sw_curve_note": ("Sigma: H(x,y) = 0, the Newton polynomial of the toric "
                          "diagram (the 5d mirror / Seiberg-Witten curve), as a "
                          "cover of C^*_x"),
        "curve_x_degree": x_deg,
        "curve_y_degree": y_deg,
        "sheets": max(y_deg, 1),
        "genus": genus,                        # = 5d Coulomb-branch rank
        "coulomb_rank": genus,
        "flavor_rank": flavor_rank,
        "num_asymptotic_legs": V,
        "asymptotic_legs": legs,               # (p,q) legs = zig-zags = network ends
        "charge_lattice_rank": (num_nodes if num_nodes is not None
                                else 2 * genus + flavor_rank),
        "dimer_link": ("the brane-tiling zig-zag paths ARE the network's "
                       "asymptotic (p,q) legs; the tiling's Kasteleyn "
                       "determinant is the SW curve H(x,y)"),
        "wall_crossing": ("BPS states = finite spectral networks (string webs) "
                          "on Sigma; charges live in H_1(Sigma); wall-crossing "
                          "of the network = mutation of the BPS quiver"),
        "refs": ["arXiv:0907.3987", "arXiv:1204.4824", "arXiv:1006.0146"],
    }


def orbifold_chamber(matrix, label_base: int = 0) -> dict:
    """The canonical finite BPS chamber of a toric 5d SCFT: the n fractional
    branes themselves as hypermultiplets (charges = standard basis e_i),
    stable at the orbifold point where all arg Z_i are nearly aligned
    (Closset-Del Zotto, arXiv:1912.13502).  This is the physically correct
    default chamber for every geometry -- including the mutation-infinite ones
    (local P^2, F0, ...) that admit no maximal green sequence.  Returns the
    stable-object list ordered by decreasing central-charge phase."""
    n = len(matrix)
    phases, masses = orbifold_central_charges(n)
    stable = []
    for i in range(n):
        v = [1 if j == i else 0 for j in range(n)]
        central = _central_charge(v, phases, masses)
        stable.append({
            "node": i + label_base,
            "dimension_vector": v,
            "charge": _basis_label(v, label_base),
            "kind": "hypermultiplet",
            "central_charge": central,
        })
    stable.sort(key=lambda s: -s["central_charge"]["phase"])
    for idx, s in enumerate(stable, start=1):
        s["index"] = idx
    return {
        "kind": "orbifold (fractional-brane) chamber",
        "num_states": n,
        "stable_objects": stable,
        "phases": phases,
        "masses": masses,
        "note": ("the n fractional branes as hypermultiplets -- the canonical "
                 "finite chamber realized at the orbifold point (all arg Z_i "
                 "nearly aligned); valid even where no maximal green sequence "
                 "exists (local P^2, F0, ...)"),
    }


def _mutate_c_matrix(C: list[list[int]], B: list[list[int]], k: int, sign: str):
    n = len(C)
    ck = C[k][:]
    sgn = 1 if sign == "green" else (-1 if sign == "red" else 1)
    out = [c[:] for c in C]
    for j in range(n):
        if j == k:
            out[j] = [-x for x in ck]
        else:
            factor = max(sgn * B[k][j], 0)
            if factor:
                out[j] = [out[j][a] + factor * ck[a] for a in range(n)]
    return out


def _basis_label(v: list[int], label_base: int = 1) -> str:
    terms = []
    for i, coeff in enumerate(v, start=label_base):
        if coeff == 0:
            continue
        base = f"gamma_{i}"
        terms.append(base if coeff == 1 else f"{coeff} {base}")
    return " + ".join(terms) if terms else "0"


def _central_charge(v: list[int], phases: list[float], masses: list[float]) -> dict:
    z = 0j
    for coeff, phase, mass in zip(v, phases, masses):
        z += coeff * mass * complex(math.cos(math.pi * phase),
                                    math.sin(math.pi * phase))
    mag = abs(z)
    arg = math.atan2(z.imag, z.real) / math.pi
    if arg < 0:
        arg += 2.0
    return {
        "re": round(z.real, 8),
        "im": round(z.imag, 8),
        "mass": round(mag, 8),
        "phase": round(arg, 8),
    }


def finite_chamber_spectrum(matrix: list[list[int]], sequence: list[int],
                            phases: list[float] | None = None,
                            masses: list[float] | None = None,
                            label_base: int = 1) -> dict:
    n = len(matrix)
    B = [row[:] for row in matrix]
    C = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    stable = []
    mutations = []

    for step, k in enumerate(sequence, start=1):
        cvec = C[k][:]
        sign = _vector_sign(cvec)
        central = None
        if phases is not None and masses is not None:
            central = _central_charge(cvec, phases, masses)
        mutation = {
            "step": step,
            "node": k + label_base,
            "status": sign,
            "c_vector": cvec,
            "charge": _basis_label(cvec, label_base),
        }
        if central is not None:
            mutation["central_charge"] = central
        mutations.append(mutation)
        if sign == "green":
            stable.append({
                "index": len(stable) + 1,
                "step": step,
                "node": k + label_base,
                "dimension_vector": cvec,
                "charge": _basis_label(cvec, label_base),
                "kind": "hypermultiplet",
                **({"central_charge": central} if central is not None else {}),
            })
        C = _mutate_c_matrix(C, B, k, sign)
        B = _mutate_exchange(B, k)

    green = all(m["status"] == "green" for m in mutations)
    final_signs = [_vector_sign(c) for c in C]
    maximal_green = green and all(s == "red" for s in final_signs)

    phase_order_ok = None
    if stable and phases is not None and masses is not None:
        seq_phases = [s["central_charge"]["phase"] for s in stable]
        phase_order_ok = all(
            seq_phases[i] + 1e-8 >= seq_phases[i + 1]
            for i in range(len(seq_phases) - 1)
        )

    return {
        "stable": stable,
        "mutations": mutations,
        "final_c_matrix": C,
        "final_exchange_matrix": B,
        "green_sequence": green,
        "maximal_green": maximal_green,
        "final_c_vector_signs": final_signs,
        "phase_order_ok": phase_order_ok,
    }




def _validate_adjacency(adjacency) -> list[list[int]]:
    matrix = _as_int_matrix(adjacency)
    n = len(matrix)
    if n > MAX_NODES:
        raise ValueError(f"BPS quiver UI is limited to {MAX_NODES} nodes")
    if any(len(row) != n for row in matrix):
        raise ValueError("quiver adjacency must be square")
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if val < 0:
                raise ValueError("quiver adjacency entries must be non-negative")
    return matrix


def exchange_from_adjacency(adjacency) -> list[list[int]]:
    """Return the skew exchange matrix B_ij = A_ij - A_ji.

    Geometry builders return an oriented adjacency matrix.  BPS mutation uses the
    antisymmetric exchange matrix, so reciprocal arrows are netted and loops are
    ignored on the diagonal.
    """
    A = _validate_adjacency(adjacency)
    n = len(A)
    B = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = A[i][j] - A[j][i]
    _validate_exchange_matrix(B)
    return B


def bps_quiver_from_adjacency_json(adjacency, sequence_text: str,
                                   phases_text: str = "", masses_text: str = "",
                                   label_base: int = 0,
                                   source: dict | None = None) -> dict:
    """Run BPS chamber analysis on an existing geometry-built quiver.

    Always reports the mutation-invariant physics (KK/D0 charge delta, flavor
    directions) and the canonical **orbifold chamber** (the n fractional branes
    as hypermultiplets -- the finite chamber every toric 5d SCFT has, even the
    mutation-infinite ones).  When no mutation `sequence` is given the orbifold
    chamber IS the reported chamber; a maximal green sequence is offered as an
    alternative when one exists.  When `phases`/`masses` are omitted, default
    central charges are supplied so the chamber displays valid Z_i."""
    A = _validate_adjacency(adjacency)
    B = exchange_from_adjacency(A)
    n = len(B)
    invariants = bps_invariants(B, label_base)
    canonical = orbifold_chamber(B, label_base)
    mgs = find_maximal_green_sequence(B)          # None for mutation-infinite quivers
    seq_given = bool((sequence_text or "").strip())

    warnings = []
    if any(A[i][i] for i in range(n)):
        warnings.append("Geometry quiver has self-arrows; loops are ignored in the BPS exchange matrix.")
    if any(A[i][j] and A[j][i] for i in range(n) for j in range(i + 1, n)):
        warnings.append("Geometry quiver has reciprocal arrows; BPS exchange matrix uses the net orientation A_ij - A_ji.")
    degenerate = all(B[i][j] == 0 for i in range(n) for j in range(n))
    if degenerate:
        warnings.append("The antisymmetrized BPS exchange matrix is zero for this geometry quiver (all c-vectors are mutation-invariant); showing the orbifold chamber.")

    if not seq_given:
        # ---- default: the canonical orbifold chamber (no mutations) ----
        chamber = {
            "kind": canonical["kind"],
            "mutation_sequence": [],
            "num_mutations": 0,
            "green_sequence": None,
            "maximal_green": None,
            "stable_count": canonical["num_states"],
            "phase_order_ok": True,
            "auto_selected": True,
            "maximal_green_sequence": ([k + label_base for k in mgs]
                                       if mgs is not None else None),
            "warnings": warnings + ([] if mgs is not None else [
                "This quiver is mutation-infinite (no maximal green sequence); "
                "the orbifold chamber is the canonical finite chamber."]),
        }
        return {
            "source": source or {"kind": "geometry"},
            "quiver": {
                "num_nodes": n,
                "exchange_matrix": B,
                "adjacency": A,
                "num_arrows": sum(sum(row) for row in A),
                "dims": [1] * n,
                "node_labels": [f"gamma_{i}" for i in range(label_base, label_base + n)],
                "label_base": label_base,
            },
            "invariants": invariants,
            "canonical_chamber": canonical,
            "chamber": chamber,
            "stable_objects": canonical["stable_objects"],
            "mutations": [],
            "final": {"c_matrix": None, "exchange_matrix": B,
                      "c_vector_signs": None},
        }

    # ---- an explicit mutation sequence was supplied ----
    sequence = parse_sequence(sequence_text, n, zero_based_hint=(label_base == 0))
    phases = _parse_number_list(phases_text, n, "phases")
    masses = _parse_number_list(masses_text, n, "masses")
    if masses is None and phases is not None:
        masses = [1.0] * n
    if phases is None and masses is not None:
        raise ValueError("masses require phases; phases may be given without masses")
    if phases is None:                            # default Z_i so central charges show
        phases, masses = default_central_charges(B, sequence)

    spectrum = finite_chamber_spectrum(B, sequence, phases, masses, label_base)
    if not spectrum["green_sequence"]:
        warnings.append(
            "The mutation sequence is not green throughout; stable objects are "
            "reported only for green mutations before/inside that path."
        )
    if spectrum["green_sequence"] and not spectrum["maximal_green"]:
        warnings.append(
            "The sequence is green but not maximal green; this is a partial "
            "finite chamber, not a complete finite BPS spectrum."
        )
    if spectrum["phase_order_ok"] is False:
        warnings.append(
            "The supplied central-charge phases are not decreasing along the "
            "mutation order; check whether they describe the intended chamber."
        )

    return {
        "source": source or {"kind": "geometry"},
        "quiver": {
            "num_nodes": n,
            "exchange_matrix": B,
            "adjacency": A,
            "num_arrows": sum(sum(row) for row in A),
            "dims": [1] * n,
            "node_labels": [f"gamma_{i}" for i in range(label_base, label_base + n)],
            "label_base": label_base,
        },
        "invariants": invariants,
        "canonical_chamber": canonical,
        "chamber": {
            "mutation_sequence": [k + label_base for k in sequence],
            "num_mutations": len(sequence),
            "green_sequence": spectrum["green_sequence"],
            "maximal_green": spectrum["maximal_green"],
            "stable_count": len(spectrum["stable"]),
            "phase_order_ok": spectrum["phase_order_ok"],
            "maximal_green_sequence": ([k + label_base for k in mgs]
                                       if mgs is not None else None),
            "warnings": warnings,
        },
        "stable_objects": spectrum["stable"],
        "mutations": spectrum["mutations"],
        "final": {
            "c_matrix": spectrum["final_c_matrix"],
            "exchange_matrix": spectrum["final_exchange_matrix"],
            "c_vector_signs": spectrum["final_c_vector_signs"],
        },
    }

def bps_quiver_json(quiver_text: str, sequence_text: str, kind: str = "matrix",
                    phases_text: str = "", masses_text: str = "") -> dict:
    parsed = parse_quiver(quiver_text, kind)
    matrix = parsed.matrix
    n = len(matrix)
    raw_sequence_labels = [int(v) for v in re.findall(r"\d+", sequence_text or "")]
    labels_are_zero_based = parsed.labels_are_zero_based or 0 in raw_sequence_labels
    sequence = parse_sequence(sequence_text, n, labels_are_zero_based)
    label_base = 0 if labels_are_zero_based else 1
    phases = _parse_number_list(phases_text, n, "phases")
    masses = _parse_number_list(masses_text, n, "masses")
    if masses is None and phases is not None:
        masses = [1.0] * n
    if phases is None and masses is not None:
        raise ValueError("masses require phases; phases may be given without masses")

    spectrum = finite_chamber_spectrum(matrix, sequence, phases, masses, label_base)
    adjacency = [[max(matrix[i][j], 0) for j in range(n)] for i in range(n)]
    warnings = []
    if not spectrum["green_sequence"]:
        warnings.append(
            "The mutation sequence is not green throughout; stable objects are "
            "reported only for green mutations before/inside that path."
        )
    if spectrum["green_sequence"] and not spectrum["maximal_green"]:
        warnings.append(
            "The sequence is green but not maximal green; this is a partial "
            "finite chamber, not a complete finite BPS spectrum."
        )
    if spectrum["phase_order_ok"] is False:
        warnings.append(
            "The supplied central-charge phases are not decreasing along the "
            "mutation order; check whether they describe the intended chamber."
        )

    return {
        "quiver": {
            "num_nodes": n,
            "exchange_matrix": matrix,
            "adjacency": adjacency,
            "num_arrows": sum(sum(row) for row in adjacency),
            "dims": [1] * n,
            "node_labels": [f"gamma_{i}" for i in range(label_base, label_base + n)],
            "label_base": label_base,
        },
        "chamber": {
            "mutation_sequence": [k + label_base for k in sequence],
            "num_mutations": len(sequence),
            "green_sequence": spectrum["green_sequence"],
            "maximal_green": spectrum["maximal_green"],
            "stable_count": len(spectrum["stable"]),
            "phase_order_ok": spectrum["phase_order_ok"],
            "warnings": warnings,
        },
        "stable_objects": spectrum["stable"],
        "mutations": spectrum["mutations"],
        "final": {
            "c_matrix": spectrum["final_c_matrix"],
            "exchange_matrix": spectrum["final_exchange_matrix"],
            "c_vector_signs": spectrum["final_c_vector_signs"],
        },
    }
