"""Drive the three steps sequentially and report the result."""

from __future__ import annotations

import numpy as np

from .chartable import build_character_table
from .conformal import conformal_manifold_dim
from .groups import MatrixGroup, make_group
from .quiver import build_quiver


def _fmt_c(z: complex) -> str:
    z = complex(z)
    re = 0.0 if abs(z.real) < 1e-4 else z.real
    im = 0.0 if abs(z.imag) < 1e-4 else z.imag
    if im == 0.0:
        return f"{int(round(re))}" if abs(re - round(re)) < 1e-4 else f"{re:.2f}"
    return f"{re:.2f}{im:+.2f}i"


def run(group: MatrixGroup | str, verbose: bool = True) -> dict:
    """Run Step 1 -> 2 -> 3 for `group` (a MatrixGroup or a library name)."""
    if isinstance(group, str):
        group = make_group(group)

    out = lambda *a: print(*a) if verbose else None

    # ---- Step 1: the group --------------------------------------------------
    out("=" * 68)
    out(f"STEP 1  Group  Gamma = {group.name}")
    out("-" * 68)
    out(f"  {group.description}")
    out(f"  |Gamma| = {group.order}    abelian = {group.is_abelian()}")
    out(f"  generators ({len(group.generators)}):")
    for k, gen in enumerate(group.generators):
        out(f"    g{k} =")
        for row in np.round(gen, 4):
            out("        [" + "  ".join(f"{z.real:+.3f}{z.imag:+.3f}i" for z in row) + "]")

    # ---- Step 2: the quiver -------------------------------------------------
    table = build_character_table(group)
    quiver = build_quiver(group, table)
    out("=" * 68)
    out("STEP 2  McKay quiver")
    out("-" * 68)
    out(f"  #irreps (= #gauge nodes) r = {table.num_irreps}")
    out(f"  irrep dimensions           = {table.dims}")
    out(f"  node gauge groups          = "
        + ", ".join(quiver.node_rank_label(i) for i in range(quiver.num_nodes)))
    out(f"  chi_Q on classes           = [" +
        ", ".join(_fmt_c(z).strip() for z in table.chi_Q) + "]")
    out(f"  total bifundamentals       = {quiver.num_arrows}")
    out(f"  connected                  = {quiver.is_connected()}")
    out("  adjacency  A[i->j] (#arrows):")
    for i, row in enumerate(quiver.adjacency.astype(int)):
        out(f"    R{i}: [" + ",".join(f"{x:2d}" for x in row) + "]")
    out(f"  cubic superpotential terms Tr(A^3) = {quiver.num_cubic_terms()}")

    # ---- Step 3: conformal manifold ----------------------------------------
    result = conformal_manifold_dim(group)
    out("=" * 68)
    out("STEP 3  Conformal manifold (Leigh-Strassler / GKSTW)")
    out("-" * 68)
    out("  " + str(result).replace("\n", "\n  "))
    out("=" * 68)

    return {
        "group": group,
        "table": table,
        "quiver": quiver,
        "result": result,
        "dim_conf": result.dim_conf,
    }
