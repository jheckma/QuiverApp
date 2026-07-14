"""Hasse diagrams of the Higgs mechanism for magnetic quivers.

The Higgs branch of a 4d N=2 theory (= Coulomb branch of the magnetic quiver) is a
symplectic singularity; its symplectic leaves form a poset whose Hasse diagram has
edges labeled by the transverse elementary slice of each minimal degeneration
(Bourget-Cabrera-Grimminger-Hanany-Sperling-Zafrir-Zhong, arXiv:1908.04245).

Two routes are provided:

  * ``nilpotent_hasse`` -- for a nilpotent orbit closure the poset is exactly the
    dominance order on partitions, with edges labeled by Kraft-Procesi transitions
    (Cabrera-Hanany, arXiv:1609.07798).  Robust and fully combinatorial.
  * ``quiver_subtraction_hasse`` -- the general quiver-subtraction algorithm for an
    arbitrary unitary quiver (Phase 6, best-effort; imports ``magnetic``).

Import direction is strictly ``partitions <- magnetic <- hasse``; ``magnetic`` never
imports this module.
"""

from __future__ import annotations

from . import partitions as P


# ---------------------------------------------------------------------------
# nilpotent-orbit Hasse diagram (dominance poset + Kraft-Procesi slice labels)
# ---------------------------------------------------------------------------

def nilpotent_hasse(kind: str, N: int, top=None) -> dict:
    """Poset of nilpotent orbits of the classical algebra on the size-N defining rep.

    If ``top`` is a partition, the poset is restricted to the closure of that orbit
    (all orbits ``mu <= top`` in dominance) -- i.e. the Higgs-branch Hasse diagram of
    the theory whose Higgs branch is O-bar_top.  Otherwise the full nilcone is used.

    Returns nodes (one per valid orbit partition), covering edges labeled by the
    transverse slice, and a layering by quaternionic dimension.
    """
    key = P._kind_key(kind)
    universe = P.orbits(kind, N)
    if not universe:
        raise ValueError(f"no valid {key} orbits for N={N}")
    if top is not None:
        top = P.normalize(top)
        universe = [p for p in universe if P.dominates(top, p)]

    # id each orbit; order by quaternionic dimension (largest orbit on top)
    dims = {part: P.orbit_dim(kind, part) // 2 for part in universe}
    order = sorted(universe, key=lambda p: (dims[p], p))
    idx = {part: i for i, part in enumerate(order)}

    warnings = []
    nodes = []
    for part in order:
        nodes.append({
            "id": idx[part],
            "partition": list(part),
            "label": P.orbit_label(part),
            "dim_H": dims[part],
        })

    edges = []
    rels = [(a, b) for a in universe for b in universe
            if a != b and P.covers(a, b, universe)]
    for a, b in rels:
        edge = {"from": idx[a], "to": idx[b], "dim_H": dims[a] - dims[b]}
        if key == "su":
            slc = P.kp_slice_type(a, b)
            edge.update({"slice": slc["label"], "letter": slc["letter"], "k": slc["k"]})
        else:
            # orthosymplectic Kraft-Procesi slice types (a/A/b/c/d/D/...) are not yet
            # classified here; report the transverse dimension only.
            edge.update({"slice": f"[{dims[a] - dims[b]}]", "letter": "?",
                         "k": dims[a] - dims[b]})
        edges.append(edge)
    if key != "su":
        warnings.append(
            "orthosymplectic slice-type identification is not implemented; edges show "
            "the transverse quaternionic dimension only."
        )

    layers = _layers_by_dim(nodes)
    return {
        "nodes": nodes,
        "edges": edges,
        "layers": layers,
        "warnings": warnings,
    }


def _layers_by_dim(nodes) -> list:
    """Group node ids into layers of equal quaternionic dimension, ascending."""
    by_dim: dict = {}
    for nd in nodes:
        by_dim.setdefault(nd["dim_H"], []).append(nd["id"])
    return [by_dim[d] for d in sorted(by_dim)]


# ---------------------------------------------------------------------------
# JSON entry point (consumed by conformalmanifold.api)
# ---------------------------------------------------------------------------

def hasse_json(kind: str = "su", N: int = 0, quiver_text: str = "",
               method: str = "partition", partition: str = "") -> dict:
    """Higgs-branch Hasse diagram for the frontend.

    When ``partition`` is given the poset is restricted to that orbit's closure (the
    Higgs branch of the selected theory); otherwise the full nilcone is returned.
    """
    if method == "subtraction":
        from . import magnetic as Mg  # noqa: F401  (Phase 6)
        return quiver_subtraction_hasse(quiver_text)
    N = int(N)
    if N < 1:
        raise ValueError("N must be a positive integer")
    top = None
    if partition:
        from . import magnetic as Mg
        top = Mg.parse_partition(partition, N)
        if not P.valid_partition(kind, top):
            raise ValueError(f"{list(top)} is not a valid {P._kind_key(kind)} orbit")
    data = nilpotent_hasse(kind, N, top=top)
    return {
        "available": True,
        "method": "partition",
        "algebra": {"kind": P._kind_key(kind), "N": N},
        "top": list(top) if top else None,
        "nodes": data["nodes"],
        "edges": data["edges"],
        "layers": data["layers"],
        "warnings": data["warnings"],
    }


def quiver_subtraction_hasse(quiver_text: str) -> dict:
    """Higgs-branch Hasse diagram of a hand-typed quiver via quiver subtraction.

    Implemented for balanced linear A-type quivers by recognizing them as
    nilpotent-orbit quivers M(lambda) and returning the exact Kraft-Procesi poset of
    the orbit closure; other quivers are reported as not yet supported.
    """
    from . import magnetic as Mg
    quiver = Mg.parse_quiver(quiver_text)
    rec = Mg.recognize_orbit(quiver)
    if rec is None:
        raise ValueError(
            "Hasse via quiver subtraction currently supports balanced linear A-type "
            "(nilpotent-orbit) quivers; this quiver was not recognized. Non-orbit "
            "quivers (affine, orthosymplectic, forks) are future work."
        )
    kind, N, part = rec
    data = nilpotent_hasse(kind, N, top=part)
    return {
        "available": True,
        "method": "subtraction",
        "algebra": {"kind": kind, "N": N},
        "recognized_as": {"partition": list(part), "label": P.orbit_label(part)},
        "top": list(part),
        "nodes": data["nodes"],
        "edges": data["edges"],
        "layers": data["layers"],
        "warnings": data["warnings"] + [
            f"quiver recognized as M({P.orbit_label(part)}) of su({N}); showing the "
            "Kraft-Procesi orbit-closure Hasse diagram."
        ],
    }
