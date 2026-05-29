"""Step 1 -- finite subgroups Gamma < SU(3) as explicit 3x3 matrix groups.

A subgroup of SU(3) *is* the action on C^3 (treated as a point set): every
element is a 3x3 complex matrix and the group is generated from a handful of
generators by taking the closure under multiplication.

The library follows the standard naming conventions for the classification of
finite subgroups of SU(3) (Miller-Blichfeldt-Dickson; see also Yau-Yu and the
physics reviews):

    Z_n(a,b,c)      cyclic, diag(w^a, w^b, w^c),  a+b+c = 0 mod n
    Delta(3 n^2)    (Z_n x Z_n) :| Z_3            -- "trihedral" series
    Delta(6 n^2)    (Z_n x Z_n) :| S_3            -- "trihedral" series
    A4 = Delta(12)  tetrahedral          (= Delta(3*2^2))
    A5 = Sigma(60)  icosahedral / Valentiner-type
    ...

New groups are added simply by supplying SU(3) generator matrices.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np

# matrices are compared up to this tolerance when building the closure
_TOL = 1e-7
_ROUND = 6


def _key(m: np.ndarray) -> tuple:
    """Hashable fingerprint of a matrix, robust to floating point noise."""
    flat = np.round(m, _ROUND).ravel()
    # normalise -0.0 to 0.0 so the key is canonical
    return tuple((float(z.real) + 0.0, float(z.imag) + 0.0) for z in flat)


def _is_su3(m: np.ndarray) -> bool:
    if m.shape != (3, 3):
        return False
    unitary = np.allclose(m.conj().T @ m, np.eye(3), atol=1e-6)
    sdet = np.isclose(np.linalg.det(m), 1.0, atol=1e-6)
    return unitary and sdet


@dataclass
class MatrixGroup:
    """A finite matrix group, generated from `generators` (3x3 complex)."""

    name: str
    generators: list[np.ndarray]
    elements: list[np.ndarray]
    description: str = ""

    @property
    def order(self) -> int:
        return len(self.elements)

    @property
    def Q(self) -> str:
        """Label of the defining (action on C^3) representation."""
        return "Q (defining 3-dim rep, = action on C^3)"

    def is_abelian(self) -> bool:
        gs = self.generators
        for a in gs:
            for b in gs:
                if not np.allclose(a @ b, b @ a, atol=1e-6):
                    return False
        return True

    def chi_Q(self, g: np.ndarray) -> complex:
        """Character of the defining rep = trace of the matrix."""
        return complex(np.trace(g))

    def fix_Q(self, g: np.ndarray) -> int:
        """Number of unit eigenvalues of g on C^3 (dim of its fixed subspace)."""
        eig = np.linalg.eigvals(g)
        return int(np.sum(np.isclose(eig, 1.0 + 0.0j, atol=1e-5)))


def closure(generators: list[np.ndarray], name: str = "", description: str = "",
            max_order: int = 5000) -> MatrixGroup:
    """Generate the finite group spanned by `generators` (breadth-first)."""
    for g in generators:
        if not _is_su3(g):
            raise ValueError(f"generator is not in SU(3):\n{np.round(g, 4)}")

    ident = np.eye(3, dtype=complex)
    elements = {_key(ident): ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for a in frontier:
            for gen in generators:
                p = a @ gen
                k = _key(p)
                if k not in elements:
                    elements[k] = p
                    nxt.append(p)
                    if len(elements) > max_order:
                        raise ValueError(
                            f"group exceeded max_order={max_order}; "
                            "is it really finite / are generators correct?")
        frontier = nxt
    return MatrixGroup(name=name or "Gamma",
                       generators=list(generators),
                       elements=list(elements.values()),
                       description=description)


# --------------------------------------------------------------------------
# generator helpers
# --------------------------------------------------------------------------
def _w(n: int) -> complex:
    return cmath.exp(2j * math.pi / n)


def cyclic(n: int, weights: tuple[int, int, int]) -> MatrixGroup:
    """Z_n acting as diag(w^a, w^b, w^c).  Requires a+b+c = 0 mod n (CY/SUSY)."""
    a, b, c = weights
    if (a + b + c) % n != 0:
        raise ValueError(
            f"weights {weights} violate the Calabi-Yau condition "
            f"a+b+c = 0 mod {n} (got {(a + b + c) % n}); "
            "Gamma would not lie in SU(3).")
    if math.gcd(math.gcd(math.gcd(a % n, b % n), c % n), n) != 1:
        raise ValueError(
            f"weights {weights} give a NON-faithful Z_{n} action "
            f"(gcd(a,b,c,n) = {math.gcd(math.gcd(math.gcd(a % n, b % n), c % n), n)}"
            f" != 1): the matrix image is the smaller faithful quotient "
            f"Z_{n // math.gcd(math.gcd(math.gcd(a % n, b % n), c % n), n)}. "
            "The closed form below assumes a faithful Gamma; reduce the "
            "weights/order accordingly, or handle the non-faithful case "
            "separately.")
    w = _w(n)
    gen = np.diag([w ** a, w ** b, w ** c]).astype(complex)
    return closure([gen], name=f"Z_{n}({a},{b},{c})",
                   description=f"cyclic orbifold C^3/Z_{n}, weights {weights}")


# permutation matrix cycling the three coordinates -- the Z_3 of the Delta series
_E = np.array([[0, 1, 0],
               [0, 0, 1],
               [1, 0, 0]], dtype=complex)


def delta_3n2(n: int) -> MatrixGroup:
    """Delta(3 n^2) = (Z_n x Z_n) :| Z_3,  order 3 n^2."""
    w = _w(n)
    A = np.diag([w, w ** (-1), 1]).astype(complex)      # first Z_n
    B = np.diag([1, w, w ** (-1)]).astype(complex)      # second Z_n
    return closure([A, B, _E], name=f"Delta(3*{n}^2)={3 * n * n}",
                   description=f"trihedral series Delta(3 n^2), n={n}, "
                               "(Z_n x Z_n) semidirect Z_3")


def delta_6n2(n: int) -> MatrixGroup:
    """Delta(6 n^2) = (Z_n x Z_n) :| S_3,  order 6 n^2."""
    w = _w(n)
    A = np.diag([w, w ** (-1), 1]).astype(complex)
    B = np.diag([1, w, w ** (-1)]).astype(complex)
    # the extra S_3 transposition, signed to keep det = +1 (stay in SU(3))
    F = np.array([[-1, 0, 0],
                  [0, 0, -1],
                  [0, -1, 0]], dtype=complex)
    return closure([A, B, _E, F], name=f"Delta(6*{n}^2)={6 * n * n}",
                   description=f"trihedral series Delta(6 n^2), n={n}, "
                               "(Z_n x Z_n) semidirect S_3")


def tetrahedral_A4() -> MatrixGroup:
    """A4 = Delta(12) = Delta(3*2^2), the tetrahedral group, order 12."""
    g = delta_3n2(2)
    g.name = "A4 = Delta(12)"
    g.description = "tetrahedral group A4 = Delta(3*2^2), order 12"
    return g


def octahedral_S4() -> MatrixGroup:
    """S4 = Delta(24) = Delta(6*2^2), the octahedral group, order 24."""
    g = delta_6n2(2)
    g.name = "S4 = Delta(24)"
    g.description = "octahedral group S4 = Delta(6*2^2), order 24"
    return g


def icosahedral_A5() -> MatrixGroup:
    """A5 = Sigma(60), the icosahedral group, order 60 (3-dim irrep in SU(3))."""
    # standard 3-dim icosahedral generators (golden-ratio realisation)
    phi = (1 + math.sqrt(5)) / 2
    # order-2 and order-5 generators of A5 in its 3-dimensional irrep
    s = np.array([[-1, 0, 0],
                  [0, -1, 0],
                  [0, 0, 1]], dtype=complex)
    half = 0.5
    t = np.array([[half, -phi / 2, (phi - 1) / 2],
                  [phi / 2, (phi - 1) / 2, -half],
                  [(phi - 1) / 2, half, phi / 2]], dtype=complex)
    return closure([s, t], name="A5 = Sigma(60)",
                   description="icosahedral group A5 = Sigma(60), order 60")


# --------------------------------------------------------------------------
# the user-facing library
# --------------------------------------------------------------------------
def library() -> dict[str, callable]:
    """Named constructors.  Values are zero-arg (or curried) factories."""
    return {
        # a few ready-made cyclic examples
        "Z3(1,1,1)": lambda: cyclic(3, (1, 1, 1)),
        "Z5(1,1,3)": lambda: cyclic(5, (1, 1, 3)),
        "Z6(1,2,3)": lambda: cyclic(6, (1, 2, 3)),
        "Z10(2,3,5)": lambda: cyclic(10, (2, 3, 5)),
        # non-abelian
        "A4 = Delta(12)": tetrahedral_A4,
        "S4 = Delta(24)": octahedral_S4,
        "Delta(27)": lambda: delta_3n2(3),
        "Delta(54)": lambda: delta_6n2(3),
        "Delta(75)": lambda: delta_3n2(5),
        "A5 = Sigma(60)": icosahedral_A5,
    }


def list_groups() -> list[str]:
    """Names available via make_group / the interactive selector."""
    return list(library().keys())


def make_group(name: str) -> MatrixGroup:
    """Build a group from the library by name."""
    lib = library()
    if name not in lib:
        raise KeyError(f"unknown group '{name}'. Available: {list_groups()}")
    return lib[name]()
