# superconformal

A calculator of 4d N=1 **superconformal observables** for the SCFTs of D3-branes
probing C³/Γ (and SU(N) quiver gauge theories more generally). It complements
`conformalmanifold/`: where that returns the conformal-manifold dimension, this
returns the superconformal data of the same quiver.

From the quiver Lagrangian (gauge ranks, matter arrows, single-trace W) it computes:

- the superconformal R-symmetry by **a-maximization** (exact, including irrational);
- the **central charges** a, c and the **'t Hooft anomalies** (R, flavor, baryonic);
- **Hofman–Maldacena** and unitarity diagnostics;
- the **superconformal index** (rational R; and irrational R as a symbolic τ-series).

Main entry points (in `sqcdkit/`):

```python
from sqcdkit import scft_observables, index_series, index_symbolic
```

## Examples

`examples/` has self-contained runnable scripts — the C³/Γ quivers conifold,
dP0 (C³/Z₃), SPP, and C³/A4 (tetrahedral), plus flavored SQCD and Kutasov.
Run them with this folder on the path:

```
PYTHONPATH=superconformal python superconformal/examples/c3_a4.py
```

## Tests

```
pytest superconformal/tests/
```

Validated against textbook results, and against arXiv:2310.15792 for the
non-abelian orbifold C³/A4.

## Dependencies

`sympy>=1.10` and `mpmath` (added to the repo `requirements.txt`).
