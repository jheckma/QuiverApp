"""Step 2 -- the McKay quiver of the C^3/Gamma orbifold gauge theory.

Following Douglas-Moore / Lawrence-Nekrasov-Vafa / Kachru-Silverstein:

    * one gauge node per irrep R_i of Gamma; gauge group U(N * dim R_i)
      (or SU(...)); the ranks are N * dim R_i.
    * chiral bifundamentals: the number of arrows i -> j is

          a_{ij} = mult of R_j in (Q (x) R_i) = < chi_Q chi_i , chi_j >

      where Q is the defining 3-dim rep (the action on C^3).
    * the cubic superpotential descends from W_{N=4} = Tr Phi^1[Phi^2,Phi^3];
      its terms are the Gamma-invariant closed 3-loops of the quiver,
      counted (with orientation) by Tr(A^3) on the adjacency matrix A.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chartable import CharacterTable, build_character_table
from .groups import MatrixGroup


@dataclass
class McKayQuiver:
    table: CharacterTable
    adjacency: np.ndarray             # A[i, j] = #arrows from node i to node j
    dims: list[int]                   # dim R_i  (rank multiplier of node i)

    @property
    def num_nodes(self) -> int:
        return self.adjacency.shape[0]

    @property
    def num_arrows(self) -> int:
        return int(round(self.adjacency.sum().real))

    def node_rank_label(self, i: int) -> str:
        d = self.dims[i]
        return f"U({d}N)" if d > 1 else "U(N)"

    def num_cubic_terms(self) -> int:
        """Oriented closed 3-loops = Tr(A^3): the cubic superpotential terms."""
        A = self.adjacency
        return int(round(np.trace(A @ A @ A).real))

    def is_connected(self) -> bool:
        A = self.adjacency
        n = A.shape[0]
        reach = (A + A.T + np.eye(n)) > 0
        for _ in range(n):
            reach = (reach @ reach) > 0
        return bool(reach.all())


def build_quiver(group: MatrixGroup, table: CharacterTable | None = None) -> McKayQuiver:
    if table is None:
        table = build_character_table(group)

    r = table.num_irreps
    sizes = np.array(table.class_sizes, dtype=float)
    G = float(group.order)
    chars = table.chars
    chi_Q = table.chi_Q

    A = np.zeros((r, r))
    for i in range(r):
        for j in range(r):
            # < chi_Q chi_i , chi_j > = (1/|G|) sum_k |C_k| chi_Q(k) chi_i(k) conj(chi_j(k))
            val = np.sum(sizes * chi_Q * chars[i] * np.conj(chars[j])) / G
            A[i, j] = round(val.real)

    return McKayQuiver(table=table, adjacency=A, dims=table.dims)
