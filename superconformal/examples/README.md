# Examples

Self-contained, runnable scripts. Each builds one theory as a plain dict and
prints its protected observables. Run them from the repo root with this folder
on the path (`PYTHONPATH=superconformal`).

| Script | Theory |
|---|---|
| `sqcd.py` | SQCD: SU(Nc) with Nf flavors. `R_Q = 1 - Nc/Nf`, central charges, flavor `SU(Nf)^3`, and the s-confinement index of SU(2), Nf=3 |
| `conifold.py` | Klebanov-Witten SU(2) x SU(2). Central charges and the lowest-order index (4 mesons + 6 baryons) |
| `kutasov.py` | SU(2) + adjoint + flavors, `W = Tr Phi^3`. Rational R fixed by the superpotential |
| `spp.py` | Suspended pinch point: non-orbifold toric quiver with an adjoint. Irrational R (`sqrt(97)`), exact a, c |
| `c3_a4.py` | C^3/A4 tetrahedral non-abelian orbifold (McKay quiver). All fields R=2/3, `a = 3N^2 - 19/24` |
| `irrational_index.py` | Toric quivers with irrational a-max R. The superconformal index as a symbolic `tau`-series (`index_symbolic`) |

Run, e.g.:

```bash
PYTHONPATH=superconformal python superconformal/examples/sqcd.py
```
