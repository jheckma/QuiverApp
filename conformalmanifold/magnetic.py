"""Magnetic quivers for 4d N=2 (and 5d/6d) theories.

A *magnetic quiver* for a hyper-Kaehler space X is a 3d N=4 quiver gauge theory Q
whose **Coulomb branch** equals X.  For a 4d N=2 theory T the Higgs branch is a
symplectic singularity that is realized as the Coulomb branch of such a Q, so the
whole Higgs-branch structure (Hilbert series, global symmetry, Hasse diagram) can be
computed combinatorially from Q with the **monopole formula** (Cremonesi-Hanany-
Zaffaroni, arXiv:1309.2657).

This module provides:
  * a ``Quiver`` data model (unitary / orthosymplectic gauge nodes, flavors, edges),
  * ``monopole_hilbert_series`` -- the Coulomb-branch Hilbert series as a lattice sum,
  * ``magnetic_quiver_from_partition`` -- the linear quiver M(lambda) whose Coulomb
    branch is the closure of a nilpotent orbit O_lambda of a classical group,
  * JSON entry points consumed by ``conformalmanifold.api``.

Conventions (the "2Delta" convention, as in most implementations):
  * a half-hypermultiplet scalar has dimension 1/2 and contributes t^1, so Hilbert
    series are power series in t with integer exponents;
  * HS(t) = sum_{m dominant} t^{2 Delta(m)} * P_residual(t; m), with
        2 Delta(m) = - 2 sum_{alpha in Delta+} |alpha(m)|
                     + sum_{hypers} sum_{rho in R} |rho(m)| ;
  * for a U(N) node whose flux m has value-multiplicities {n_k} (residual gauge group
    prod_k U(n_k)) the dressing factor is
        P = prod_k prod_{j=1}^{n_k} 1 / (1 - t^{2j}) .

Orthosymplectic (O/USp) nodes and half-integer/spinor lattices are handled by the
gated code path in ``magnetic_ortho`` (Phase 5); the unitary path below is exact and
locked by the fixtures in ``tests/test_magnetic.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from . import partitions as P


MAX_NODES = 24
MAX_RANK = 12
MAX_LATTICE_CUTOFF = 40           # hard cap on |m_i| when summing the monopole lattice
MAX_MONOPOLE_LEAVES = 400_000     # work budget: total flux configurations summed (~2s)


# ---------------------------------------------------------------------------
# quiver data model
# ---------------------------------------------------------------------------

@dataclass
class Quiver:
    """A 3d N=4 quiver.

    ``gauge`` is a list of ``(kind, rank)`` with ``kind`` in ``{'U','SU','O','USp'}``
    (``'SU'`` currently only for rank 2; ``'O'``/``'USp'`` via the gated ortho path).
    ``flavors[i]`` is the number of fundamental hypermultiplets on gauge node ``i``.
    ``edges`` maps an unordered pair ``(i, j)`` with ``i < j`` to the bifundamental
    multiplicity (number of hypers).  ``half_edges`` records edges that are
    half-hypermultiplets (relevant only for orthosymplectic nodes).
    """

    gauge: list          # [(kind:str, rank:int)]
    flavors: list        # [int] per gauge node
    edges: dict          # {(i,j): multiplicity}, i < j
    half_edges: set = field(default_factory=set)   # subset of edges that are half-hypers
    node_labels: list = field(default_factory=list)

    # -- structure ---------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.gauge)

    def rank(self, i: int) -> int:
        return self.gauge[i][1]

    def kind(self, i: int) -> str:
        return self.gauge[i][0]

    def neighbors(self, i: int):
        """Yield (j, multiplicity) for every gauge node joined to node i."""
        for (a, b), m in self.edges.items():
            if a == i:
                yield b, m
            elif b == i:
                yield a, m

    def balance(self, i: int) -> int:
        """e_i = (adjacent gauge ranks + flavors) - 2 * rank, per unitary convention."""
        incoming = self.flavors[i] + sum(m * self.rank(j) for j, m in self.neighbors(i))
        return incoming - 2 * self.rank(i)

    def coulomb_dim(self) -> int:
        """Quaternionic dimension of the Coulomb branch = sum of gauge-group ranks.

        For OSp/SU the stored ``rank`` is the defining-rep size, so the Cartan rank is
        U(r)->r, SU(r)->r-1, USp(2n)->n, SO(N)->floor(N/2).
        """
        tot = 0
        for k, r in self.gauge:
            if k == "U":
                tot += r
            elif k == "SU":
                tot += r - 1
            else:                      # USp(2n) / SO(N): r is the defining-rep size
                tot += r // 2
        return tot

    def is_unitary(self) -> bool:
        return all(k in ("U", "SU") for k, _r in self.gauge)

    def to_adjacency(self) -> list:
        """Symmetric gauge-gauge adjacency matrix (bifundamental multiplicities)."""
        n = self.n
        adj = [[0] * n for _ in range(n)]
        for (a, b), m in self.edges.items():
            adj[a][b] = adj[b][a] = m
        return adj

    def labels(self) -> list:
        if self.node_labels:
            return list(self.node_labels)
        out = []
        for k, r in self.gauge:
            out.append(str(r) if k == "U" else f"{k}({r})")
        return out


def linear_quiver(ranks, flavors=None, kind="U") -> Quiver:
    """Build a linear (A-type) unitary quiver from a list of gauge ranks."""
    ranks = [int(r) for r in ranks]
    n = len(ranks)
    if flavors is None:
        flavors = [0] * n
    flavors = [int(f) for f in flavors]
    if len(flavors) != n:
        raise ValueError("flavors must have one entry per gauge node")
    edges = {(i, i + 1): 1 for i in range(n - 1)}
    gauge = [(kind, r) for r in ranks]
    return Quiver(gauge=gauge, flavors=flavors, edges=edges)


# ---------------------------------------------------------------------------
# power-series arithmetic (integer coefficients, truncated at `order`)
# ---------------------------------------------------------------------------

def _series_geom(degree: int, order: int) -> list:
    """1 / (1 - t^degree) truncated at t^order."""
    out = [0] * (order + 1)
    k = 0
    while k <= order:
        out[k] = 1
        k += degree
    return out


def _series_mul(a: list, b: list, order: int) -> list:
    out = [0] * (order + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        if i > order:
            break
        for j in range(0, order - i + 1):
            bj = b[j]
            if bj:
                out[i + j] += ai * bj
    return out


def _dressing_series(degrees, order: int) -> list:
    """Product of 1/(1 - t^{2d}) over residual-group Casimir degrees d."""
    out = [0] * (order + 1)
    out[0] = 1
    for d in degrees:
        out = _series_mul(out, _series_geom(2 * d, order), order)
    return out


# ---------------------------------------------------------------------------
# per-node group primitives (unitary path)
# ---------------------------------------------------------------------------

def _dominant_u_fluxes(N: int, cutoff: int):
    """Non-increasing integer tuples of length N with entries in [-cutoff, cutoff]."""
    def rec(remaining, hi):
        if remaining == 0:
            yield ()
            return
        for v in range(hi, -cutoff - 1, -1):
            for rest in rec(remaining - 1, v):
                yield (v,) + rest
    yield from rec(N, cutoff)


def _u_root_term(m) -> int:
    """sum_{i<j} |m_i - m_j| for a U(N)/SU(N) flux (positive roots e_i - e_j)."""
    total = 0
    for i in range(len(m)):
        for j in range(i + 1, len(m)):
            total += abs(m[i] - m[j])
    return total


def _u_dressing_degrees(m):
    """Residual-group Casimir degrees for a U(N) flux: prod U(n_k) -> degrees 1..n_k."""
    degrees = []
    run = 1
    for i in range(1, len(m) + 1):
        if i < len(m) and m[i] == m[i - 1]:
            run += 1
        else:
            degrees.extend(range(1, run + 1))
            run = 1
    return degrees


# --- orthosymplectic per-node primitives (gated, validated cases only) ------
#
# 2Delta conventions (Cremonesi-Hanany-Zaffaroni 1309.2657 classical specialization;
# dressing tables Hanany-Kalveks 1807.02521):
#   * USp(2n)=C_n and SO(2n+1)=B_n share the integer magnetic lattice m_1>=..>=m_n>=0.
#   * Positive roots (both B and C): {e_i +- e_j : i<j}; plus the long root 2e_i (C) or
#     the short root e_i (B).
#   * The vector rep has weights +-e_i (and one 0 for SO odd), so the fundamental
#     weight sum is 2*sum|m_i|.
#   * Dressing: each nonzero flux value of multiplicity k gives a U(k) factor (degrees
#     1..k); the p zeros give C_p / B_p (degrees 2,4,..,2p).
#   * SO(2n)=D_n (even) needs m_n possibly negative and the Pfaffian degree -> GATED.

def _ortho_family(kind: str, size: int) -> str:
    """Cartan family B/C/D for an orthosymplectic node of given defining-rep size."""
    if kind == "USp":
        return "C"                              # USp(2n) = C_n
    return "B" if size % 2 == 1 else "D"         # SO/O(2n+1)=B_n, SO/O(2n)=D_n


def _dominant_bc_fluxes(n: int, cutoff: int):
    """Non-increasing non-negative integer n-tuples m_1>=..>=m_n>=0 (B_n/C_n chamber)."""
    def rec(remaining, hi):
        if remaining == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(remaining - 1, v):
                yield (v,) + rest
    yield from rec(n, cutoff)


def _dominant_d_fluxes(n: int, cutoff: int):
    """D_n dominant chamber m_1>=..>=m_{n-1}>=|m_n|: the last charge may be negative."""
    if n == 1:
        for m in range(cutoff, -cutoff - 1, -1):
            yield (m,)
        return

    def head(remaining, hi):                    # length n-1, non-increasing, in [0, hi]
        if remaining == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in head(remaining - 1, v):
                yield (v,) + rest
    for h in head(n - 1, cutoff):
        last = h[-1]
        for mn in range(last, -last - 1, -1):
            yield h + (mn,)


def _bc_root_term(family: str, m) -> int:
    """sum_{alpha>0} |alpha(m)| for the B/C/D family."""
    total = 0
    for i in range(len(m)):
        for j in range(i + 1, len(m)):
            total += abs(m[i] - m[j]) + abs(m[i] + m[j])      # roots e_i +- e_j
    if family == "C":            # USp: long roots 2 e_i
        total += sum(2 * abs(x) for x in m)
    elif family == "B":          # SO(odd): short roots e_i
        total += sum(abs(x) for x in m)
    # D: only e_i +- e_j
    return total


def _bc_dressing_degrees(family: str, m):
    """Residual Casimir degrees: U(k) on each nonzero value block, B_p/C_p/D_p on zeros.

    For D-type the last charge may be negative, but the residual group depends on the
    flux only through the multiplicities of |m_i| (m_N and -m_N give the same block),
    so multiplicities are counted in (m_1,...,m_{N-1},|m_N|) (2505.03875 eq. 5.7)."""
    if family == "D" and m:
        m = tuple(m[:-1]) + (abs(m[-1]),)
    degrees = []
    nonzero = [x for x in m if x != 0]
    zeros = len(m) - len(nonzero)
    run, prev = 0, None
    for x in nonzero:                            # U(k) blocks on equal nonzero values
        if x == prev:
            run += 1
        else:
            if prev is not None:
                degrees.extend(range(1, run + 1))
            prev, run = x, 1
    if prev is not None:
        degrees.extend(range(1, run + 1))
    if family == "D":            # SO(even): degrees 2,4,..,2(p-1) and the Pfaffian p
        degrees.extend(2 * j for j in range(1, zeros))
        if zeros >= 1:
            degrees.append(zeros)
    else:                         # B/C: degrees 2,4,..,2p
        degrees.extend(2 * j for j in range(1, zeros + 1))
    return degrees


def _fund_weight_sum(kind: str, m) -> int:
    """sum_rho |rho(m)| for one fundamental (vector) hyper on a node of the given kind."""
    if kind in ("U", "SU"):
        return sum(abs(x) for x in m)          # weights e_i
    return 2 * sum(abs(x) for x in m)          # USp/SO vector: weights +-e_i (0 -> 0)


def _signed_vals(kind: str, size: int, m):
    """Weight-values of the (self-dual for OSp) vector/fundamental rep on flux m."""
    if kind in ("U", "SU"):
        return list(m)                          # fundamental: weights e_i
    vals = []
    for x in m:
        vals.append(x)
        vals.append(-x)                          # vector: weights +-e_i
    if kind in ("SO", "O") and size % 2 == 1:
        vals.append(0)                           # SO(odd) vector has a zero weight
    return vals


def _bifund_weight_sum(ka, sa, ma, kb, sb, mb, half: bool) -> int:
    """sum_rho |rho(m)| for a bifundamental hyper between two nodes (half-hyper halves it).

    Uses R_a (x) conj(R_b); for U this reproduces the |m_i - n_j| differences, for the
    self-dual orthosymplectic vector the +- signs give the standard OSp weight sums.
    """
    A = _signed_vals(ka, sa, ma)
    B = _signed_vals(kb, sb, mb)
    if kb in ("U", "SU"):
        B = [-x for x in B]                      # conjugate fundamental of the U node
    s = sum(abs(a + b) for a in A for b in B)
    if half:
        if s % 2:
            raise ValueError("half-hypermultiplet weight sum is not an integer")
        s //= 2
    return s


def _node_fluxes(kind: str, rank: int, cutoff: int):
    """Iterate (flux, root_term, dressing_degrees, weight) for one gauge node.

    ``rank`` stores the defining-rep size (U(N)->N, USp(2n)->2n, SO(N)->N)."""
    if kind == "U":
        for m in _dominant_u_fluxes(rank, cutoff):
            yield m, _u_root_term(m), _u_dressing_degrees(m), 1
    elif kind in ("USp", "SO", "O"):
        fam = _ortho_family(kind, rank)
        n = rank // 2
        gen = _dominant_d_fluxes(n, cutoff) if fam == "D" else _dominant_bc_fluxes(n, cutoff)
        for m in gen:
            yield m, _bc_root_term(fam, m), _bc_dressing_degrees(fam, m), 1
    elif kind == "SU" and rank == 2:
        # SU(2): the traceless flux (m, -m) with m >= 0 spans the coweight lattice;
        # the positive root e_1 - e_2 gives alpha(m) = 2m and the fundamental doublet
        # has weights +-1.  Residual is SU(2) (Casimir degree 2) at m = 0 and U(1)
        # (degree 1) for m > 0.
        for mm in range(0, cutoff + 1):
            m = (mm, -mm)
            degrees = [2] if mm == 0 else [1]
            yield m, _u_root_term(m), degrees, 1
    else:
        raise ValueError(
            f"unitary monopole path supports U(N) and SU(2); got {kind}({rank})"
        )


# ---------------------------------------------------------------------------
# the monopole formula (Coulomb-branch Hilbert series)
# ---------------------------------------------------------------------------

def _matter_2delta(quiver: Quiver, fluxes) -> int:
    """sum over hypers of sum_rho |rho(m)|  (the +matter part of 2 Delta)."""
    total = 0
    # fundamental flavors (kind-aware: U/SU weights e_i; USp/SO vector weights +-e_i)
    for i, m in enumerate(fluxes):
        w = quiver.flavors[i]
        if w:
            total += w * _fund_weight_sum(quiver.kind(i), m)
    # bifundamentals (kind-aware; OSp-OSp edges default to half-hypermultiplets)
    for (a, b), mult in quiver.edges.items():
        if mult == 0:
            continue
        half = (a, b) in quiver.half_edges
        total += mult * _bifund_weight_sum(
            quiver.kind(a), quiver.rank(a), fluxes[a],
            quiver.kind(b), quiver.rank(b), fluxes[b], half)
    return total


def _twice_delta(quiver: Quiver, fluxes, root_terms) -> int:
    root = sum(root_terms)                       # sum_alpha |alpha(m)| over all nodes
    return -2 * root + _matter_2delta(quiver, fluxes)


def monopole_hilbert_series(quiver: Quiver, order: int = 12):
    """Coulomb-branch Hilbert series as an integer coefficient list [c0, c1, ..., c_order].

    Returns ``(coeffs, info)`` where ``info`` carries diagnostics: the minimal 2Delta
    over non-trivial fluxes (good/ugly/bad classification), whether the truncated sum
    converged, the lattice cutoff reached, and any warnings.
    """
    if quiver.n > MAX_NODES:
        raise ValueError(f"too many gauge nodes ({quiver.n} > {MAX_NODES})")
    warnings = []
    osp = ("USp", "SO", "O")
    has_o = any(k == "O" for k, _r in quiver.gauge)
    has_sosp = any(k in ("USp", "SO") for k, _r in quiver.gauge)
    has_ortho_bifund = any(
        quiver.kind(a) in osp and quiver.kind(b) in osp for (a, b) in quiver.edges)
    for k, r in quiver.gauge:
        if k in ("U", "SU") and r > MAX_RANK:
            raise ValueError(f"gauge rank {r} exceeds MAX_RANK={MAX_RANK}")
        if k == "SU" and r != 2:
            raise ValueError("only SU(2) is supported among special-unitary nodes")
    if has_o:
        warnings.append(
            "O(N) is EXPERIMENTAL: only the connected SO lattice is summed; the "
            "charge-conjugation +/- sectors, the O-projection and the gauged-flavor "
            "flux phase (arXiv:2505.03875) are not implemented -- treat as exploratory."
        )
    if has_sosp:
        warnings.append(
            "orthosymplectic node (USp(2n)=C, SO(2n+1)=B, SO(2n)=D): the monopole "
            "lattice and dressing are validated against arXiv:2505.03875 "
            "(USp(2)+3f = C^2/Z4; SO(4)+3v; USp(2)=SU(2))."
        )
    if has_ortho_bifund:
        warnings.append(
            "SO/USp bifundamentals are treated as half-hypermultiplets (the pseudoreal "
            "convention, validated); full-quiver results with gauged O(N) flavor may "
            "still need the arXiv:2505.03875 flux-phase correction."
        )
    node_specs = [(quiver.kind(i), quiver.rank(i)) for i in range(quiver.n)]

    # classification (good/ugly/bad) from the minimal monopole -- cheap unit-flux sample
    min_2delta = _sample_min_two_delta(quiver, node_specs)
    classification = _classify(min_2delta)
    if classification == "bad":
        warnings.append(
            "some non-trivial monopole has 2*Delta <= 0: this is a bad theory and the "
            "Coulomb-branch Hilbert series is not reliable."
        )

    def sum_at_cutoff(per_node):
        coeffs = [0] * (order + 1)

        def recurse(idx, chosen_roots, chosen_degrees, fluxes):
            if idx == quiver.n:
                two_delta = _twice_delta(quiver, fluxes, chosen_roots)
                if two_delta < 0 or two_delta > order:
                    return
                dress = [0] * (order + 1)
                dress[0] = 1
                for dd in chosen_degrees:
                    dress = _series_mul(dress, _dressing_series(dd, order), order)
                for p in range(0, order + 1 - two_delta):
                    if dress[p]:
                        coeffs[two_delta + p] += dress[p]
                return
            for m, rt, dd, _wt in per_node[idx]:
                recurse(idx + 1, chosen_roots + [rt], chosen_degrees + [dd], fluxes + [m])

        recurse(0, [], [], [])
        return coeffs

    # adaptive cutoff: grow until the coefficients up to `order` stabilise, but bail out
    # (gracefully) once the flux lattice would exceed the work budget.  The cost of a
    # cutoff is exactly the product of the per-node flux-box sizes, so we can estimate it
    # before doing the expensive per-leaf work.
    cutoff = 1
    prev = None
    converged = False
    coeffs = [0] * (order + 1)
    spent = 0
    too_expensive = False
    while cutoff <= MAX_LATTICE_CUTOFF:
        per_node = [list(_node_fluxes(k, r, cutoff)) for k, r in node_specs]
        leaves = 1
        for pn in per_node:
            leaves *= len(pn)
        if spent + leaves > MAX_MONOPOLE_LEAVES:
            too_expensive = True
            break
        spent += leaves
        coeffs = sum_at_cutoff(per_node)
        if prev is not None and coeffs == prev:
            converged = True
            break
        prev = coeffs
        cutoff += 1

    if too_expensive:
        info = {
            "min_two_delta": min_2delta,
            "classification": classification,
            "converged": False,
            "too_expensive": True,
            "cutoff": cutoff,
            "reason": (
                "the monopole lattice for this quiver is too large to sum in the "
                "browser (gauge ranks too high); the quiver, global symmetry and Hasse "
                "diagram are still exact. Try a smaller orbit or a lower series order."
            ),
            "warnings": warnings,
        }
        return None, info

    if not converged:
        warnings.append(
            "monopole sum did not stabilise within the lattice cutoff; the theory may "
            "be bad (Coulomb branch not a cone) or the requested order is too high."
        )
    info = {
        "min_two_delta": min_2delta,
        "classification": classification,
        "converged": converged,
        "too_expensive": False,
        "cutoff": cutoff,
        "warnings": warnings,
    }
    return coeffs, info


def _sample_min_two_delta(quiver: Quiver, node_specs):
    """Minimal 2*Delta over unit (|m_i| <= 1) non-trivial fluxes -- fixes the
    good/ugly/bad classification (the minimal monopole is a unit flux)."""
    per_node = [list(_node_fluxes(k, r, 1)) for k, r in node_specs]
    best = [None]

    def recurse(idx, roots, fluxes):
        if idx == quiver.n:
            if any(any(v != 0 for v in m) for m in fluxes):
                td = _twice_delta(quiver, fluxes, roots)
                if best[0] is None or td < best[0]:
                    best[0] = td
            return
        for m, rt, _dd, _wt in per_node[idx]:
            recurse(idx + 1, roots + [rt], fluxes + [m])

    recurse(0, [], [])
    return best[0]


def _classify(min_2delta):
    if min_2delta is None:
        return "trivial"
    if min_2delta <= 0:
        return "bad"
    if min_2delta == 1:
        return "ugly"
    return "good"


# ---------------------------------------------------------------------------
# global symmetry from balanced nodes
# ---------------------------------------------------------------------------

def balanced_nodes(quiver: Quiver) -> list:
    return [i for i in range(quiver.n) if quiver.kind(i) in ("U", "SU")
            and quiver.balance(i) == 0]


def coulomb_global_symmetry(quiver: Quiver) -> dict:
    """Best-effort Coulomb-branch global symmetry from the balanced-node subgraph.

    Connected chains of balanced unitary nodes enhance to SU(len+1); each gauge node
    also carries a topological U(1).  Returns a structured description plus a display
    string.  (Non-simply-laced / exceptional enhancements are not detected.)
    """
    bal = set(balanced_nodes(quiver))
    # connected components of the balanced subgraph
    seen, comps = set(), []
    for start in bal:
        if start in seen:
            continue
        stack, comp = [start], []
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            comp.append(v)
            for j, _m in quiver.neighbors(v):
                if j in bal and j not in seen:
                    stack.append(j)
        comps.append(sorted(comp))

    nonabelian = []
    for comp in comps:
        # a linear chain of k balanced nodes -> SU(k+1)
        sub_edges = sum(1 for (a, b) in quiver.edges
                        if a in comp and b in comp)
        if sub_edges == len(comp) - 1:
            nonabelian.append(("SU", len(comp) + 1))
        else:
            nonabelian.append(("?", len(comp)))
    n_u1 = quiver.n - len(bal)
    factors = [f"SU({r})" if kind == "SU" else f"G[{r}]" for kind, r in nonabelian]
    if n_u1 > 0:
        factors.append(f"U(1)^{n_u1}" if n_u1 > 1 else "U(1)")
    display = " x ".join(factors) if factors else "trivial"
    return {
        "balanced_nodes": sorted(bal),
        "nonabelian_factors": nonabelian,
        "u1_count": n_u1,
        "display": display,
    }


# ---------------------------------------------------------------------------
# magnetic quiver of a nilpotent orbit closure (partition -> quiver)
# ---------------------------------------------------------------------------

def _a_type_ranks(part) -> list:
    """Balanced ranks N_i = sum_j (min(i,j) - i j / N) w_j for the A-type orbit quiver.

    w_i = lambda_i - lambda_{i+1}; the A_{N-1} inverse-Cartan Green's function gives the
    unique node ranks making every gauge node balanced.  Returns ranks for i=1..N-1.
    """
    part = P.normalize(part)
    N = sum(part)
    lam = [part[i] if i < len(part) else 0 for i in range(N + 1)]
    w = [lam[i] - lam[i + 1] for i in range(N)]         # w[0..N-1], indices 1..N used
    ranks = []
    for i in range(1, N):
        val = Fraction(0)
        for j in range(1, N):
            val += (min(i, j) - Fraction(i * j, N)) * w[j - 1]
        if val.denominator != 1:
            raise ValueError(f"non-integer rank N_{i} = {val} for partition {part}")
        ranks.append(int(val))
    return ranks, [w[i - 1] for i in range(1, N)]


def magnetic_quiver_from_partition(kind: str, part) -> Quiver:
    """Magnetic quiver whose Coulomb branch is the nilpotent orbit closure O_part.

    A-type (``su``) only in this unitary path: a linear quiver of N-1 U(N_i) nodes with
    w_i = lambda_i - lambda_{i+1} fundamental flavors on node i.  Orthosymplectic
    (``so``/``sp``) partitions are handled by ``magnetic_ortho`` (gated).
    """
    key = P._kind_key(kind)
    if key != "su":
        raise ValueError(
            "orthosymplectic partition quivers use the gated ortho path (Phase 5)"
        )
    part = P.normalize(part)
    if not part:
        raise ValueError("empty partition")
    ranks, flavors = _a_type_ranks(part)
    # keep the contiguous support; interior zeros would disconnect (shouldn't occur)
    keep = [i for i, r in enumerate(ranks) if r > 0]
    if not keep:
        # the zero orbit (1^N): the quiver is a point.
        return Quiver(gauge=[], flavors=[], edges={})
    lo, hi = keep[0], keep[-1]
    if any(ranks[i] == 0 for i in range(lo, hi + 1)):
        raise ValueError(f"disconnected orbit quiver for partition {part}")
    ranks = ranks[lo:hi + 1]
    flavors = flavors[lo:hi + 1]
    return linear_quiver(ranks, flavors, kind="U")


# ---------------------------------------------------------------------------
# parsing and formatting helpers
# ---------------------------------------------------------------------------

def _same_linear(a: "Quiver", b: "Quiver") -> bool:
    """Do two linear unitary quivers coincide (up to left-right reflection)?"""
    if a.n != b.n or not a.is_unitary() or not b.is_unitary():
        return False
    ar = [r for _k, r in a.gauge]
    br = [r for _k, r in b.gauge]
    af, bf = list(a.flavors), list(b.flavors)
    return (ar == br and af == bf) or (ar == br[::-1] and af == bf[::-1])


def recognize_orbit(quiver: "Quiver", cap: int = 24):
    """Recognise a linear unitary quiver as the magnetic quiver M(lambda) of an
    su-nilpotent orbit, returning ``('su', N, partition)`` or ``None``.

    A non-trivial orbit closure has quaternionic dimension >= N-1 (the minimal orbit),
    so only partitions of N <= dim_H + 1 can match -- this keeps the search tiny.
    """
    if quiver.n == 0 or not quiver.is_unitary():
        return None
    for (a, b), m in quiver.edges.items():
        if b != a + 1 or m != 1:
            return None                        # not a simple linear chain
    target_dim = quiver.coulomb_dim()
    for N in range(2, min(cap, target_dim + 1) + 1):
        for part in P.partitions_of(N):
            try:
                mq = magnetic_quiver_from_partition("su", part)
            except ValueError:
                continue
            if mq.coulomb_dim() != target_dim:
                continue
            if _same_linear(mq, quiver):
                return ("su", N, part)
    return None


def identify_theory(quiver: "Quiver") -> dict:
    """Cross-reference a magnetic quiver against known 3d N=4 moduli spaces.

    Returns ``{"name", "detail"}`` naming the Coulomb branch when recognized
    (nilpotent orbit closures, T[SU(N)], Kleinian singularities, free space)."""
    if quiver.n == 0:
        return {"name": "trivial", "detail": "the Coulomb branch is a point."}
    # single U(1) with f flavors -> C^2/Z_f (Kleinian A_{f-1}); f=1 is free C^2
    if (quiver.n == 1 and quiver.kind(0) == "U" and quiver.rank(0) == 1
            and not quiver.edges):
        f = quiver.flavors[0]
        if f == 1:
            return {"name": "free hypermultiplet (C^2)",
                    "detail": "U(1) with one flavor; Coulomb branch = C^2."}
        if f >= 2:
            return {"name": f"Kleinian singularity C^2/Z_{f} (A_{f - 1})",
                    "detail": f"U(1) with {f} flavors; the reduced moduli space of one "
                              f"SU(2) instanton on C^2/Z_{f} = the A_{f - 1} du Val slice."}
    rec = recognize_orbit(quiver)
    if rec is not None:
        _kind, N, part = rec
        label = P.orbit_label(part)
        if part == (N,):
            return {"name": f"T[SU({N})]  (nilcone of sl({N}))",
                    "detail": f"self-mirror Gaiotto-Witten theory; Coulomb branch = the "
                              f"closure of the regular nilpotent orbit {label} of sl({N})."}
        if part == (2,) + (1,) * (N - 2):
            return {"name": f"minimal nilpotent orbit of sl({N})",
                    "detail": f"reduced one-SU({N})-instanton moduli space; magnetic "
                              f"quiver [1]-(1)^{N - 2}-[1], Coulomb branch O-bar_{label}."}
        return {"name": f"nilpotent orbit O_{label} of sl({N})",
                "detail": f"Coulomb branch = closure of the nilpotent orbit {label} of "
                          f"sl({N}); quaternionic dimension {quiver.coulomb_dim()}."}
    return {"name": "3d N=4 quiver (unrecognized)",
            "detail": f"Coulomb branch is a {quiver.coulomb_dim()}-quaternionic-"
                      "dimensional symplectic singularity."}


def parse_partition(text: str, N: int | None = None) -> tuple:
    """Parse ``"2,1,1"`` / ``"[2,1,1]"`` / ``"2 1 1"`` into a partition tuple.

    If ``N`` is given, the parts must sum to ``N`` (else it is inferred).
    """
    if text is None:
        raise ValueError("empty partition")
    cleaned = text.strip().strip("[]() ")
    if not cleaned:
        raise ValueError("empty partition")
    tokens = [tok for tok in cleaned.replace(",", " ").split() if tok]
    try:
        parts = [int(tok) for tok in tokens]
    except ValueError:
        raise ValueError(f"partition must be integers; got {text!r}")
    if any(p <= 0 for p in parts):
        raise ValueError("partition parts must be positive integers")
    part = P.normalize(parts)
    if N is not None and sum(part) != N:
        raise ValueError(f"partition {list(part)} does not sum to N={N}")
    return part


_GAUGE_RE = None
_FLAVOR_RE = None


def parse_quiver(text: str) -> Quiver:
    """Parse a linear 3d N=4 quiver from text such as ``1-2-[3]`` or ``[1]-1-1-[1]``.

    Tokens are separated by ``-``.  A bare integer (optionally ``U``/``SU`` prefixed,
    parens allowed) is a gauge node U(n); ``[k]`` is a flavor node attaching ``k``
    fundamentals to the adjacent gauge node.  Consecutive gauge nodes are joined by a
    single bifundamental.  (Linear unitary quivers only; forks/loops are future work.)
    """
    import re
    global _GAUGE_RE, _FLAVOR_RE
    if _GAUGE_RE is None:
        _GAUGE_RE = re.compile(r"(U|SU|USp|SO|O)?\(?(\d+)\)?")
        _FLAVOR_RE = re.compile(r"\[(\d+)\]")
    if not text or not text.strip():
        raise ValueError("empty quiver")
    tokens = [t.strip() for t in text.strip().split("-") if t.strip()]
    if not tokens:
        raise ValueError("empty quiver")
    seq = []
    for tok in tokens:
        mf = _FLAVOR_RE.fullmatch(tok)
        mg = _GAUGE_RE.fullmatch(tok)
        if mf:
            seq.append(("f", int(mf.group(1))))
        elif mg:
            rank = int(mg.group(2))
            if rank <= 0:
                raise ValueError(f"gauge rank must be positive in token {tok!r}")
            seq.append(("g", mg.group(1) or "U", rank))
        else:
            raise ValueError(f"cannot parse quiver token {tok!r}")
    gauge = [(item[1], item[2]) for item in seq if item[0] == "g"]
    if not gauge:
        raise ValueError("quiver has no gauge nodes")
    if len(gauge) > MAX_NODES:
        raise ValueError(f"too many gauge nodes ({len(gauge)} > {MAX_NODES})")
    flavors = [0] * len(gauge)
    gi = -1
    for item in seq:
        if item[0] == "g":
            gi += 1
        else:
            flavors[max(gi, 0)] += item[1]
    edges = {(i, i + 1): 1 for i in range(len(gauge) - 1)}
    # a bifundamental between two orthosymplectic nodes is a real rep -> half-hyper
    osp = {"USp", "SO", "O"}
    half_edges = {(i, i + 1) for i in range(len(gauge) - 1)
                  if gauge[i][0] in osp and gauge[i + 1][0] in osp}
    return Quiver(gauge=gauge, flavors=flavors, edges=edges, half_edges=half_edges)


def monopole_series_json(quiver_text: str, order: int = 12) -> dict:
    """Coulomb-branch data for a hand-typed 3d N=4 quiver (custom-input mode)."""
    order = _clamp_order(order)
    quiver = parse_quiver(quiver_text)
    warnings = []
    coeffs, info = monopole_hilbert_series(quiver, order)
    warnings.extend(info.get("warnings", []))
    sym = coulomb_global_symmetry(quiver)
    return {
        "available": True,
        "mode": "custom",
        "quiver": _quiver_json(quiver),
        "physics": {
            "coulomb_dim": quiver.coulomb_dim(),
            "global_symmetry": sym,
            "classification": info["classification"],
            "hilbert_series": _hilbert_series_json(coeffs, info, order),
            "identified": identify_theory(quiver),
        },
        "orbit": None,
        "warnings": warnings,
    }


def _hilbert_series_json(coeffs, info, order: int) -> dict:
    """Package the monopole-formula result, or an 'omitted' notice when the lattice
    sum exceeded the work budget."""
    if coeffs is None:
        return {
            "available": False,
            "convention": "2Delta (scalar ~ t)",
            "order": order,
            "reason": info.get("reason", "Hilbert series unavailable."),
        }
    return {
        "available": True,
        "convention": "2Delta (scalar ~ t)",
        "order": order,
        "coeffs": coeffs,
        "string": hilbert_series_string(coeffs),
        "converged": info["converged"],
    }


def hilbert_series_string(coeffs, var: str = "t") -> str:
    """Render a coefficient list as ``1 + 8 t^2 + 27 t^4 + ... + O(t^n)``."""
    terms = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        elif k == 1:
            terms.append(f"{var}" if c == 1 else f"{c} {var}")
        else:
            coef = "" if c == 1 else f"{c} "
            terms.append(f"{coef}{var}^{k}")
    body = " + ".join(terms) if terms else "0"
    return f"{body} + O({var}^{len(coeffs)})"


def _quiver_json(quiver: Quiver) -> dict:
    return {
        "num_nodes": quiver.n,
        "gauge": [[k, r] for k, r in quiver.gauge],
        "flavors": list(quiver.flavors),
        "edges": [[a, b, m] for (a, b), m in sorted(quiver.edges.items())],
        "adjacency": quiver.to_adjacency(),
        "ranks": [r for _k, r in quiver.gauge],
        "dims": [1] * quiver.n,
        "labels": quiver.labels(),
        "balances": [quiver.balance(i) for i in range(quiver.n)],
        "coulomb_dim": quiver.coulomb_dim(),
    }


# ---------------------------------------------------------------------------
# JSON entry points (consumed by conformalmanifold.api)
# ---------------------------------------------------------------------------

def magnetic_from_partition_json(kind: str, N: int, partition_text: str,
                                 order: int = 12) -> dict:
    """Build M(lambda) for a nilpotent orbit and compute its Coulomb-branch data."""
    order = _clamp_order(order)
    part = parse_partition(partition_text, N)
    key = P._kind_key(kind)
    if not P.valid_partition(kind, part):
        raise ValueError(
            f"{list(part)} is not a valid {key} nilpotent-orbit partition"
        )
    warnings = []
    orbit = {
        "partition": list(part),
        "label": P.orbit_label(part),
        "transpose": list(P.transpose(part)),
        "complex_dim": P.orbit_dim(kind, part),
        "quaternionic_dim": P.orbit_dim_quaternionic(kind, part),
    }
    quiver_json = None
    physics = None
    if key == "su":
        quiver = magnetic_quiver_from_partition(kind, part)
        coeffs, info = monopole_hilbert_series(quiver, order)
        warnings.extend(info.get("warnings", []))
        sym = coulomb_global_symmetry(quiver)
        quiver_json = _quiver_json(quiver)
        physics = {
            "coulomb_dim": quiver.coulomb_dim(),
            "global_symmetry": sym,
            "classification": info["classification"],
            "hilbert_series": _hilbert_series_json(coeffs, info, order),
            "identified": identify_theory(quiver),
        }
    else:
        # orthosymplectic: the alternating O/USp magnetic quiver (Barbasch-Vogan +
        # collapse, Hanany-Kalveks) is experimental and not yet constructed here.
        warnings.append(
            "the orthosymplectic magnetic-quiver construction (alternating O/USp nodes) "
            "is experimental and not yet implemented; showing the orbit dimension and "
            "the Kraft-Procesi Hasse diagram only."
        )
    return {
        "available": True,
        "mode": "partition",
        "algebra": {"kind": key, "N": sum(part)},
        "quiver": quiver_json,
        "physics": physics,
        "orbit": orbit,
        "warnings": warnings,
    }


def _clamp_order(order) -> int:
    try:
        order = int(order)
    except (TypeError, ValueError):
        order = 12
    return max(2, min(order, 24))
