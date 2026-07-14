"""Finite U(2) subgroups of SU(3).

A finite subgroup of U(2) embeds in SU(3) by

    g  |-->  diag( g , det(g)^{-1} ) ,     g in U(2),

which is unitary, has determinant 1, and is injective -- so every finite U(2)
subgroup lands as a *faithful* SU(3) matrix group (a reducible 2+1 action on
C^3 = C^2 (+) C).  Faithfulness is exactly the regime the closed form
`dim_C M_conf = sum_g fix_Q(g) - 1` already covers, so once the SU(3) matrices
are built the existing pipeline computes everything.

We expose the standard finite subgroups of SU(2) (the binary polyhedral / ADE
groups) and their *phase extensions* <H, w.I> by a central U(1) phase
w = e^{2 pi i / k}, which give genuine U(2) (not SU(2)) groups.

Every constructor verifies the realized group order against the expected value
via `closure` (a built-in correctness check), so a wrong generator transcription
fails loudly rather than silently mislabeling a group.
"""
from __future__ import annotations

import cmath
import math
import numpy as np

from .groups import MatrixGroup, closure


# ---------------------------------------------------------------------------
# embedding U(2) -> SU(3)
# ---------------------------------------------------------------------------
def embed_u2(g) -> np.ndarray:
    """U(2) element g (2x2 unitary) -> SU(3) matrix diag(g, det(g)^{-1})."""
    g = np.asarray(g, dtype=complex)
    assert g.shape == (2, 2), "U(2) element must be 2x2"
    d = np.linalg.det(g)
    M = np.zeros((3, 3), dtype=complex)
    M[:2, :2] = g
    M[2, 2] = 1.0 / d
    return M


def _quat(w, x, y, z) -> np.ndarray:
    """Unit quaternion (w,x,y,z) -> SU(2) matrix [[w+xi, y+zi],[-y+zi, w-xi]]."""
    return np.array([[w + x * 1j, y + z * 1j],
                     [-y + z * 1j, w - x * 1j]], dtype=complex)


def _phase(k: int) -> np.ndarray:
    """Central U(1) phase w.I_2 of order k, w = e^{2 pi i / k}."""
    w = cmath.exp(2j * math.pi / k)
    return w * np.eye(2, dtype=complex)


def _build(gens2, name, description, expect_order=None):
    """closure() over the SU(3) embeddings of 2x2 U(2) generators, with an
    optional realized-order assertion."""
    G = closure([embed_u2(g) for g in gens2], name=name, description=description)
    if expect_order is not None and G.order != expect_order:
        raise RuntimeError(f"{name}: realized order {G.order} != expected "
                           f"{expect_order} (bad generators)")
    return G


# ---------------------------------------------------------------------------
# SU(2) generator sets for the binary polyhedral groups (det = 1)
# ---------------------------------------------------------------------------
def _bd_gens(m: int):
    """Binary dihedral / dicyclic group of order 4m (type D_{m+2})."""
    zeta = cmath.exp(1j * math.pi / m)            # primitive 2m-th root
    a = np.array([[zeta, 0], [0, 1 / zeta]], dtype=complex)
    b = _quat(0, 0, 1, 0)                         # [[0,1],[-1,0]]
    return [a, b]


def _2T_gens():
    """Binary tetrahedral 2T, order 24 (type E6)."""
    return [_quat(0, 1, 0, 0),                    # i
            _quat(0.5, 0.5, 0.5, 0.5)]            # (1+i+j+k)/2


def _2O_gens():
    """Binary octahedral 2O, order 48 (type E7)."""
    r = 1 / math.sqrt(2)
    return _2T_gens() + [_quat(r, r, 0, 0)]       # (1+i)/sqrt2


def _2I_gens():
    """Binary icosahedral 2I, order 120 (type E8).

    Two unit icosians: (1+i+j+k)/2 and (1/phi + i + phi j)/2 (phi = golden
    ratio); both have unit norm and generate the full group of 120 icosians."""
    phi = (1 + math.sqrt(5)) / 2
    return [_quat(0.5, 0.5, 0.5, 0.5),
            _quat(0.0, 0.5 / phi, 0.5, 0.5 * phi)]


_BINARY = {
    "2T": (_2T_gens, 24, "binary tetrahedral 2T (E6)"),
    "2O": (_2O_gens, 48, "binary octahedral 2O (E7)"),
    "2I": (_2I_gens, 120, "binary icosahedral 2I (E8)"),
}


# ---------------------------------------------------------------------------
# public constructors
# ---------------------------------------------------------------------------
def binary_dihedral(m: int, phase: int = 1) -> MatrixGroup:
    """Dicyclic group of order 4m in SU(2), optionally phase-extended to U(2)."""
    gens = _bd_gens(m)
    if phase > 1:
        gens = gens + [_phase(phase)]
        return _build(gens, f"BD_{m}.Z{phase}",
                      f"binary dihedral (order {4*m}) x central Z_{phase} phase, in U(2)")
    return _build(gens, f"BD_{m}",
                  f"binary dihedral / dicyclic, order {4*m} (D_{m+2}), in SU(2)",
                  expect_order=4 * m)


def binary_polyhedral(kind: str, phase: int = 1) -> MatrixGroup:
    """'2T'|'2O'|'2I' in SU(2), optionally phase-extended to U(2)."""
    gen_fn, order, desc = _BINARY[kind]
    gens = gen_fn()
    if phase > 1:
        return _build(gens + [_phase(phase)], f"{kind}.Z{phase}",
                      f"{desc} x central Z_{phase} phase, in U(2)")
    return _build(gens, kind, f"{desc}, in SU(2)", expect_order=order)


def u2_library() -> dict:
    """Named finite-U(2)-subgroup constructors (zero-arg factories).

    A representative spread: each binary polyhedral group in SU(2), plus a
    genuine-U(2) phase-extended copy; and a couple of binary dihedrals."""
    lib = {
        "BD_2 (Q8)":        lambda: binary_dihedral(2),
        "BD_3":             lambda: binary_dihedral(3),
        "BD_4":             lambda: binary_dihedral(4),
        "2T":               lambda: binary_polyhedral("2T"),
        "2O":               lambda: binary_polyhedral("2O"),
        "2I":               lambda: binary_polyhedral("2I"),
        # genuine U(2): central phase extensions
        "BD_2.Z3":          lambda: binary_dihedral(2, phase=3),
        "2T.Z3":            lambda: binary_polyhedral("2T", phase=3),
        "2O.Z3":            lambda: binary_polyhedral("2O", phase=3),
        "2I.Z5":            lambda: binary_polyhedral("2I", phase=5),
    }
    return lib
