# conformalmanifold

Conformal manifolds of the 4D N=1 SCFTs obtained from **N D3-branes probing a
Calabi–Yau singularity C³/Γ**, with Γ a finite subgroup of SU(3).

Given Γ, the package runs the Douglas–Moore / Lawrence–Nekrasov–Vafa /
Kachru–Silverstein construction in three steps and returns the (complex)
dimension of the conformal manifold in closed form.

## Install / run

Only dependency is `numpy`.

```
pip install numpy
python -m conformalmanifold                 # interactive group menu
python -m conformalmanifold "A4 = Delta(12)"
python -m conformalmanifold --cyclic 10 2 3 5
python -m conformalmanifold --list
python examples.py
pytest tests/
```

### Web frontend

A dependency-free single-page UI (Python stdlib server + numpy only):

```
python -m conformalmanifold --serve          # http://localhost:8000
python -m conformalmanifold.webapp --port 9000
```

Pick a preset group from the dropdown, or choose **"Custom cyclic Z_n(a,b,c)…"**
to enter your own weights. The page shows the group data, the McKay quiver
(drawn as a node-and-arrow diagram plus the adjacency matrix), and the conformal
manifold dimension. API: `GET /api/groups`, `GET /api/compute?name=...`
(or `name=__cyclic__&n=&a=&b=&c=`).

## Database / sweep

`conformalmanifold.database` runs the full pipeline over a set of groups and
records **all** computed data (no filtering — that is left to downstream sweep
code) as one row per group in a SQLite table `quivers`:

```
python -m conformalmanifold.database quivers.db        # build the default sweep
```

```python
from conformalmanifold.database import build_database, default_entries, record_for
build_database(default_entries(), "quivers.db")        # SU(3) core + U(2) subgroups
```

Columns: `name`, `family`, `description`, `grp_order`, `is_abelian`,
`third_coord_fixed` (genuine SU(2) embedding vs U(2)), `num_nodes`,
`irrep_dims`, `num_arrows`, `connected`, `num_cubic_terms`, `fixed_sum` (`S`),
`per_direction`, `dim_conf`, `note`, `generators` (3×3 complex, JSON). The build
is idempotent (upsert on `name`), so a sweep can be extended and rebuilt without
duplicating rows.

## The pipeline

**Step 1 — the group.** Γ is realised as an explicit 3×3 matrix group inside
SU(3); the matrices *are* the action on C³ (treated as a point set). Built-in
families follow standard names from the classification of finite SU(3)
subgroups:

| name | group | order |
|------|-------|-------|
| `Z_n(a,b,c)` | cyclic, `diag(ωᵃ,ωᵇ,ωᶜ)`, `a+b+c≡0 (mod n)` | `n` |
| `Delta(3n²)` | `(Zₙ×Zₙ)⋊Z₃` | `3n²` |
| `Delta(6n²)` | `(Zₙ×Zₙ)⋊S₃` | `6n²` |
| `A4 = Delta(12)` | tetrahedral | 12 |
| `S4 = Delta(24)` | octahedral | 24 |
| `A5 = Sigma(60)` | icosahedral | 60 |

You can also pass any generators of your own (`groups.closure([...])`).
The cyclic action must be **faithful** (`gcd(a,b,c,n)=1`) and Calabi–Yau
(`a+b+c≡0 mod n`); both are checked.

**Finite U(2) subgroups of SU(3)** (`conformalmanifold.u2groups`). Every finite
subgroup of U(2) embeds in SU(3) by `g ↦ diag(g, det(g)⁻¹)` (unitary, det 1,
injective ⇒ *faithful*), a reducible `2+1` action on `C³ = C² ⊕ C`. Built-in:

| name | group | order |
|------|-------|-------|
| `BD_m` | binary dihedral / dicyclic (`D_{m+2}`) | `4m` |
| `2T` | binary tetrahedral (`E6`) | 24 |
| `2O` | binary octahedral (`E7`) | 48 |
| `2I` | binary icosahedral (`E8`) | 120 |
| `H.Z_k` | `H` × central `U(1)` phase `e^{2πi/k}` (genuine U(2)) | varies |

A pure-SU(2) embedding (`det = 1`) leaves the third coordinate fixed
(`C²/Γ × C`) and gives `dim_C M_conf = |Γ| + 1`; a central phase extension
(`H.Z_k`, `k` odd) makes the action genuinely 3-dimensional. Use
`binary_dihedral(m, phase=k)` / `binary_polyhedral('2I', phase=k)` for arbitrary
members, or `embed_u2(g)` to drop in your own U(2) generators.

**Step 2 — the McKay quiver.**
- one gauge node per irrep `R_i` of Γ, with gauge group `U(N·dim R_i)`;
- bifundamental chirals: the number of arrows `i→j` is the multiplicity
  `a_ij = ⟨χ_Q·χ_i, χ_j⟩`, where `Q` is the defining 3-dim rep;
- the cubic superpotential descends from
  `W_{N=4}=Tr Φ¹[Φ²,Φ³]`; its terms are the Γ-invariant closed 3-loops,
  counted by `Tr(A³)`.

The irreducible characters are obtained numerically from the class algebra
(Burnside / Dixon), so the quiver is built for any finite matrix group.

**Step 3 — the conformal manifold.** By Green–Komargodski–Seiberg–Tachikawa–
Wecht the conformal manifold is the space of exactly marginal couplings modulo
the complexified global symmetry. These marginal couplings are **field-theory
data, not invariants of the orbifold geometry**: the holomorphic gauge couplings
(one `τ` per quiver node) and the cubic superpotential couplings
(one per Γ-invariant closed 3-loop), with the broken global symmetry generated
by the U(1)'s acting on the bifundamentals. Performing the GKSTW count
`#couplings − dim_C G_F` for these theories collapses to a closed form that can
be written as a character sum over Γ:

```
fix_Q(g) = number of unit eigenvalues of g on C³
dim_C M_conf = ( Σ_{g∈Γ} fix_Q(g) ) − 1            (faithful Γ, connected quiver)
```

The `fix_Q` sum is a convenient bookkeeping of that field-theory count — *not* a
count of fixed loci or geometric invariants of C³/Γ. For a diagonal cyclic
action it evaluates to

```
dim_C M_conf = gcd(n,a) + gcd(n,b) + gcd(n,c) − 1.
```

The result is **independent of the number of branes N**.

## Example output

```
$ python -m conformalmanifold --cyclic 10 2 3 5
STEP 3  Conformal manifold (Leigh-Strassler / GKSTW)
  S = sum_g fix_Q(g)            = 8
  per-coordinate fixed counts  = [2, 1, 5] (sum 8)
  dim_C M_conf = S - 1          = 7   (independent of N)
```

(`gcd(10,2)+gcd(10,3)+gcd(10,5)−1 = 2+1+5−1 = 7`.)

## Toric Calabi–Yau quivers (`conformalmanifold.toric`)

The orbifold pipeline above covers `C³/Γ`. The companion module
`conformalmanifold.toric` extends the same conformal-manifold count to the
broader class of **toric** Calabi–Yau three-fold singularities — the conifold,
the suspended pinch point, the `Y^{p,q}` and `L^{a,b,c}` families, and the del
Pezzo cones. Each geometry carries a **quiver gauge theory** (gauge nodes,
oriented bifundamental arrows, a two-term toric superpotential) and a **toric
diagram** (a convex lattice polygon in `Z²`, the Newton polygon of the brane
tiling's Kasteleyn determinant, defined up to `GL(2,Z)` + translation).

### Labeling scheme

Geometries are named by canonical string labels (parsed by `make_toric`). There
are two tiers:

**Explicit quiver** (`ToricQuiver` — gauge nodes, arrows, superpotential, plus
the toric diagram; dimension cross-checked both ways):

| label | geometry | `dim_C M_conf` |
|-------|----------|----------------|
| `C3` | flat `C³` (N=4 SYM as N=1) | 2 |
| `conifold` | `C(T^{1,1})`, Klebanov–Witten `= Y^{1,0}` | 3 |
| `Y(p,q)` | `Y^{p,q}`, `0<q<p` | 3 |
| `L(1,5,2)` | `L^{1,5,2}` (explicit non-`Y` member) | 3 |
| `dP0` | del Pezzo 0 `= C³/Z₃(1,1,1) = C(P²)` | 2 |
| `dP1` | del Pezzo 1 `= Y^{2,1}` | 3 |

**Diagram only** (`ToricDiagram` — any toric CY3 by its toric diagram;
`dim_C M_conf = B − 1` and `#gauge = 2·area` from the polygon, no hand-built
tiling). There are infinitely many toric CY3 — one per lattice polygon — so the
library ships a representative spread and `from_diagram` accepts any polygon:

| label | family | `#gauge` | `dim` |
|-------|--------|----------|-------|
| `dP2`, `dP3` | del Pezzo cones (`dPₙ`: `n+3` gauge, `dim n+2`) | 5, 6 | 4, 5 |
| `F0` | cone over `P¹×P¹` (Hirzebruch F0) | 4 | 3 |
| `Z(n,m)` | `C³/(Z_n×Z_m)` orbifold (`dim = n+m+gcd−1`) | `n·m` | … |
| `L(a,b,c)` | general `L^{a,b,c}` quadrilateral | `a+b` | `B−1` |

```
python -m conformalmanifold --toric conifold       # explicit quiver
python -m conformalmanifold --toric "Y(5,2)"
python -m conformalmanifold --toric dP3             # diagram only
python -m conformalmanifold --toric "Z(3,4)"
python -m conformalmanifold --list-toric
```

```python
from conformalmanifold.toric import from_diagram
g = from_diagram("my-geom", [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)])
g.dim_conf(), g.num_gauge_groups          # (5, 6)   -- any toric diagram works
```

### Convention / scope

The dimension reported is the conformal-manifold dimension counted in the
**gauge + toric-superpotential coupling space** — the "β-deformation sector" —
the *same convention the orbifold pipeline uses*. For a **generic** toric CY3
this is the full conformal manifold (Benvenuti–Hanany, hep-th/0502043: a general
`Y^{p,q}` is 3-dimensional, matching `Y(p,q)`). At **symmetry-enhanced** special
points the true conformal manifold is larger, because extra non-toric marginal
operators open up — the two known cases are flagged with a `note`:

- `conifold = Y^{1,0}`: this count gives **3**, the true dimension is **5**
  (Benvenuti–Hanany) — enhanced `SU(2)×SU(2)` flavour symmetry;
- `C³` / N=4 SYM: this count gives **2**, the true dimension is **3**
  (Leigh–Strassler `τ, β, h`).

### The dimension, two ways

For each geometry the conformal-manifold dimension is given as **two
formulations of the same toric count** (the test-suite asserts they agree — they
are two views of one combinatorial quantity, not independent derivations):

1. **Field theory (authoritative).** Leigh–Strassler / NSVZ marginal-coupling
   counting written directly on the quiver,
   ```
   dim_C M_conf = (n_gauge + n_W) − rank(M),
   ```
   where `M` is the incidence matrix of the linearised β-function conditions
   (one row per gauge node + one per superpotential term; one column per field;
   entry `1` iff the field touches that node / appears in that term). This is the
   same count `conformal.py` performs for the orbifolds, now on a general quiver.

2. **Geometry (closed form).** The number of boundary lattice points of the
   toric diagram, minus one,
   ```
   dim_C M_conf = B − 1,
   ```
   where `B` is the number of external legs of the `(p,q)` 5-brane web dual to
   the toric diagram. For `C³/Z_K(a,b,c)` the triangle has
   `B = gcd(K,a) + gcd(K,b) + gcd(K,c)`, so `B − 1` reproduces the orbifold
   character formula `Σ_g fix_Q(g) − 1` exactly — the toric statement *contains*
   the orbifold one.

`Y^{p,q}` and the smooth `L^{a,b,c}` whose toric quadrilateral has primitive
edges (`B = 4`) give `dim = 3`; geometries with longer boundaries give more.

### Database

`python -m conformalmanifold.database quivers.db` builds **three** tables in one
shot: the orbifold `quivers` table, a `toric_quivers` table (explicit-quiver
geometries: quiver sizes, `dim_conf`, the geometric `B − 1` cross-check,
lattice-point counts, and the full arrows + superpotential as JSON), and a
`toric_diagrams` table (diagram-only geometries: `#gauge = 2·area`, `dim_conf`,
boundary/interior points, edge lengths, polygon as JSON).

```python
from conformalmanifold.database import (build_toric_database,
                                         build_toric_diagram_database)
from conformalmanifold.toric import (default_toric_library,
                                     default_toric_diagram_library)
build_toric_database(default_toric_library(), "quivers.db")
build_toric_diagram_database(default_toric_diagram_library(), "quivers.db")
```

### Non-isolated singularities

The **explicit-quiver** library is restricted to **isolated** toric CY3 (plus
smooth `C³`), where the field-theory count and the geometric `B − 1` closed form
agree — the tests assert `LS == B − 1` there. Two non-isolated cases are
deliberately kept out of the *explicit-quiver* library, because their `N=2`
symmetry enhancement breaks that identity:

- `C²/Z_n × C` (the `A_{n-1}` necklace, a line of `A`-singularities): `dim = n+1`,
  which disagrees with the `N=1` incidence count `n`.
- `SPP`, `xy = z w²` (a line of `A_1` singularities).

The **diagram-only** catalog is not so restricted: it does include non-isolated
geometries (e.g. some `C³/(Z_n×Z_m)`), where only the geometric `dim = B − 1` is
reported and no field-theory `LS` cross-check is asserted — so no false identity
is claimed. The `C²/Γ × C` orbifolds are also covered by the `C²/Γ × C` branch of
the orbifold pipeline, via the character formula.

## Scope / caveats

- The closed form is for a **faithful** Γ (a genuine SU(3) subgroup), which
  gives a connected quiver. Non-faithful actions — where the abstract group has
  a kernel on C³, so the matrix image is a proper quotient — change the overall
  `−1` (the single decoupled U(1)) and are rejected by the cyclic constructor
  rather than answered incorrectly.
- The McKay character table is computed numerically; group orders are capped
  (`closure(..., max_order=5000)`) as a finiteness guard.
```
