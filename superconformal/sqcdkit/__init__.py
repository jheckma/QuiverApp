"""sqcdkit: a calculator of 4d N=1 superconformal observables.

From the Lagrangian data of an SU(N) quiver gauge theory (gauge ranks, matter
arrows with R-charges, single-trace superpotential, optional SU(N) global
flavor nodes) compute the protected observables: the superconformal R via
a-maximization (exact, including irrational), the central charges a and c,
the 't Hooft anomalies (R, flavor SU(N)^3 / SU(N)^2-R, abelian flavor), the
Hofman-Maldacena and unitarity diagnostics, and the superconformal index.

A theory is a plain dict, e.g. SQCD = SU(2) with 3 flavors:

    sqcd = {
        "ranks": [2],
        "arrows": [
            {"label": "Q",  "source": 0, "target": 1, "r_charge": "1/2"},
            {"label": "Qb", "source": 2, "target": 0, "r_charge": "1/2"},
        ],
        "superpotential": [],
    }
    obs = scft_observables(sqcd, flavor_ranks=[3, 3])   # -> R, a, c, anomalies, ...
    I   = index_series(sqcd, 4, flavor_ranks=[3, 3], derive_r="amax")

Requires sympy (and mpmath for the algebraic-number recovery).
"""

from sqcdkit.a_maximization import (
    AMaxError,
    superconformal_central_charges,
    central_charges_match,
    with_superconformal_r,
    central_charge_scft_bounds,
    one_loop_beta_coefficients,
    asymptotic_freedom_report,
    mesonic_unitarity_scan,
    flavor_thooft_anomalies,
    abelian_flavor_anomalies,
    scft_observables,
)
from sqcdkit.r_repair import repair_r_charges, RRepairError
from sqcdkit.superconformal_index import (
    SuperconformalIndexError,
    index_series,
    index_pq,
    index_matches,
    index_symbolic,
    format_index,
)

__version__ = "0.1.0"

__all__ = [
    "AMaxError",
    "RRepairError",
    "SuperconformalIndexError",
    "superconformal_central_charges",
    "central_charges_match",
    "with_superconformal_r",
    "central_charge_scft_bounds",
    "one_loop_beta_coefficients",
    "asymptotic_freedom_report",
    "mesonic_unitarity_scan",
    "flavor_thooft_anomalies",
    "abelian_flavor_anomalies",
    "scft_observables",
    "repair_r_charges",
    "index_series",
    "index_pq",
    "index_matches",
    "index_symbolic",
    "format_index",
]
