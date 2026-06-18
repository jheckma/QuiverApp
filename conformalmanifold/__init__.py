"""conformalmanifold -- conformal manifolds of C^3/Gamma orbifold quiver gauge theories.

Given a finite group Gamma < SU(3) (i.e. N D3-branes probing the Calabi-Yau
singularity C^3/Gamma), this package carries out the Douglas-Moore /
Lawrence-Nekrasov-Vafa / Kachru-Silverstein construction:

    Step 1  groups   : realise Gamma as an explicit 3x3 matrix group in SU(3).
    Step 2  quiver   : build the McKay quiver (nodes = irreps, arrows from Q).
    Step 3  conformal: count Leigh-Strassler / GKSTW exactly marginal couplings
                       and return dim_C of the conformal manifold.

The headline closed form (faithful Gamma, connected quiver) is

    dim_C M_conf = ( sum_{g in Gamma} fix_Q(g) ) - 1

where fix_Q(g) is the number of unit eigenvalues of g acting on C^3.  For a
cyclic group Z_K acting with weights (a,b,c) this reduces to the familiar

    dim_C M_conf = gcd(K,a) + gcd(K,b) + gcd(K,c) - 1.

The companion module `conformalmanifold.toric` extends this from orbifolds to the
broader class of *toric* Calabi-Yau singularities (the conifold, suspended pinch
point, Y^{p,q}, L^{a,b,c} and del Pezzo cones).  There the same count is written
directly on the toric quiver (Leigh-Strassler / NSVZ marginal-coupling counting)
and equals the geometric closed form  dim_C M_conf = B - 1, where B is the number
of boundary lattice points of the toric diagram (external legs of the (p,q) web).
"""

from .groups import MatrixGroup, library, list_groups, make_group
from .chartable import CharacterTable
from .quiver import McKayQuiver
from .conformal import conformal_manifold_dim, ConformalResult
from .pipeline import run
from .toric import ToricQuiver, make_toric, list_toric, default_toric_library

__all__ = [
    "MatrixGroup",
    "library",
    "list_groups",
    "make_group",
    "CharacterTable",
    "McKayQuiver",
    "conformal_manifold_dim",
    "ConformalResult",
    "run",
    "ToricQuiver",
    "make_toric",
    "list_toric",
    "default_toric_library",
]
