"""Seiberg duality of a quiver WITH superpotential (DWZ mutation).

A quiver-with-potential is

  arrows : dict  label -> (src, tgt)
  W      : list of (coeff, cycle: tuple[label, ...])   # closed cyclic words

Mutation mu_k at node k is the Derksen-Weyman-Zelevinsky mutation, which for
these N=1 quiver gauge theories *is* Seiberg duality on gauge node k:

  1. mesons:    for each a: i->k and b: k->j add M[ab]: i->j;
  2. reverse:   a: i->k  ->  a*: k->i ;   b: k->j  ->  b*: j->k;
  3. potential: [W] (each consecutive through-k pair (a, b) replaced by M[ab])
                +  sum_ab  M[ab] b* a*;
  4. reduce:    integrate out quadratic (mass) terms by the F-term equations.

Generalisations over the toolkit original: a W term may pass through k more
than once (all disjoint through-k pairs are replaced), coefficients are exact
rationals, and any structure the reduction cannot handle raises `WMutationError`
(the caller then falls back to adjacency-only tracking) instead of asserting.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction


class WMutationError(Exception):
    """The mutation/reduction hit a structure this implementation can't
    handle (e.g. an adjoint at the dualized node, or an ambiguous mass
    term); the caller should degrade to adjacency-only tracking."""


def canon_cycle(cycle):
    """Canonical rotation of a cyclic word (lexicographically minimal)."""
    cyc = tuple(cycle)
    rots = [cyc[i:] + cyc[:i] for i in range(len(cyc))]
    return min(rots)


def _composes(word, arrows):
    """True iff the cyclic word is a closed path: tgt(f_i) == src(f_{i+1})."""
    n = len(word)
    return all(arrows[word[i]][1] == arrows[word[(i + 1) % n]][0]
               for i in range(n))


def tiling_potential(t):
    """(arrows, W) of a BraneTiling: label -> (src, tgt) and the toric
    two-term superpotential as signed cyclic words in COMPOSITION order.

    Black-vertex (negative) terms come out of the dimer in the opposite
    cyclic orientation, so each word is reversed if it only composes
    backwards."""
    arrows = {f["label"]: (f["src"], f["tgt"]) for f in t.fields}
    W = []
    for term in t.superpotential:
        w = tuple(term["fields"])
        if not _composes(w, arrows):
            w = w[::-1]
            if not _composes(w, arrows):
                raise WMutationError(f"superpotential term {term['fields']} "
                                     "is not a closed path")
        W.append((Fraction(1 if term["sign"] > 0 else -1), w))
    return arrows, W


def mutate(arrows, W, k):
    """DWZ mutation / Seiberg duality at node k.  Returns (arrows', W')."""
    incoming = [l for l, (s, t) in arrows.items() if t == k and s != k]
    outgoing = [l for l, (s, t) in arrows.items() if s == k and t != k]
    if any(s == t == k for (s, t) in arrows.values()):
        raise WMutationError(f"adjoint at node {k}: mutation not defined here")
    other = {l: st for l, st in arrows.items() if k not in st}

    # --- step 1: mesons ---
    meson = {}
    new = dict(other)
    for a in incoming:
        for b in outgoing:
            lbl = f"M[{a}{b}]"
            meson[(a, b)] = lbl
            new[lbl] = (arrows[a][0], arrows[b][1])

    # --- step 2: reversed dual quarks ---
    rev = {}
    for a in incoming:
        rev[a] = a + "*"
        new[a + "*"] = (k, arrows[a][0])
    for b in outgoing:
        rev[b] = b + "*"
        new[b + "*"] = (arrows[b][1], k)

    # --- step 3: rewrite the potential ---
    newW = []
    inc, out = set(incoming), set(outgoing)
    for coeff, cyc in W:
        cyc = list(cyc)
        n = len(cyc)
        pairs = [m for m in range(n)
                 if cyc[m] in inc and cyc[(m + 1) % n] in out]
        # through-k pairs are disjoint (a ends at k, b starts at k; overlap
        # would need a label that is both -- an adjoint, excluded above)
        repl_at = {m: meson[(cyc[m], cyc[(m + 1) % n])] for m in pairs}
        skip = {(m + 1) % n for m in pairs}
        rebuilt = [repl_at.get(m, cyc[m]) for m in range(n) if
                   (m in repl_at or m not in skip)]
        if any(l in inc or l in out for l in rebuilt):
            raise WMutationError("a k-quark survived the rewrite (term "
                                 "passes through k non-consecutively)")
        newW.append((coeff, tuple(rebuilt)))
    for a in incoming:                       # Delta term: + M[ab] b* a*
        for b in outgoing:
            newW.append((Fraction(1), (meson[(a, b)], rev[b], rev[a])))

    return _reduce(new, newW)


def _reduce(arrows, W):
    """Integrate out quadratic (mass) terms by the F-term equations.

    For a quadratic term coeff*Tr(p q): eliminate the field of the pair whose
    only OTHER appearance is a single term (its F-term then solves for the
    partner as a path, which is substituted everywhere); both fields drop."""
    arrows = dict(arrows)
    W = [(Fraction(c), tuple(w)) for c, w in W]

    for _ in range(10000):
        quad = next(((i, c, w) for i, (c, w) in enumerate(W) if len(w) == 2),
                    None)
        if quad is None:
            break
        idx, coeff, (p, q) = quad
        if p == q:
            raise WMutationError("quadratic term with a repeated field")

        def other_terms(x):
            return [i for i, (c, w) in enumerate(W) if i != idx and x in w]

        # pick which of the pair to eliminate: need exactly one other term
        cand = [(x, y) for (x, y) in ((p, q), (q, p))
                if len(other_terms(x)) == 1
                and w_count(W[other_terms(x)[0]][1], x) == 1]
        if not cand:
            # a field appearing ONLY in the mass term: F-term sets the partner
            # to zero -- drop every term containing the partner.
            lone = [(x, y) for (x, y) in ((p, q), (q, p))
                    if len(other_terms(x)) == 0]
            if not lone:
                raise WMutationError("no integrable field in a mass term")
            x, y = lone[0]
            W = [(c, w) for i, (c, w) in enumerate(W)
                 if i != idx and y not in w]
            arrows.pop(x, None)
            arrows.pop(y, None)
            continue
        elim, keep = cand[0]
        oi = other_terms(elim)[0]
        oc, ow = W[oi]
        ow = list(ow)
        e = ow.index(elim)
        m = len(ow)
        repl = tuple(ow[(e + 1 + t) % m] for t in range(m - 1))
        factor = -oc / coeff                 # F-term of elim: keep = factor*repl
        out = []
        for i, (c, w) in enumerate(W):
            if i == idx or i == oi:
                continue
            if keep in w:
                w2, cc = [], c
                for sym in w:
                    if sym == keep:
                        cc *= factor
                        w2.extend(repl)
                    else:
                        w2.append(sym)
                out.append((cc, tuple(w2)))
            else:
                out.append((c, tuple(w)))
        W = out
        arrows.pop(elim, None)
        arrows.pop(keep, None)
    else:
        raise WMutationError("mass-term reduction did not terminate")

    combo = Counter()
    for c, w in W:
        combo[canon_cycle(w)] += c
    Wfinal = [(c, w) for w, c in combo.items() if c != 0]
    # every arrow in W must still exist; every dropped arrow must be gone
    for _, w in Wfinal:
        for l in w:
            if l not in arrows:
                raise WMutationError(f"reduced W references dropped field {l}")
    return arrows, Wfinal


def w_count(word, x):
    return sum(1 for l in word if l == x)


def adjacency_of(arrows, n):
    """Adjacency matrix (n x n) of a label->(src,tgt) arrow dict."""
    A = [[0] * n for _ in range(n)]
    for (s, t) in arrows.values():
        A[s][t] += 1
    return A


def w_json(W):
    """Serialise the potential for the web API: coefficients exact, integers
    shown as integers."""
    out = []
    for c, w in sorted(W, key=lambda cw: cw[1]):
        c = Fraction(c)
        coeff = str(c.numerator) if c.denominator == 1 else str(c)
        out.append({"coeff": coeff, "fields": list(w)})
    return out
