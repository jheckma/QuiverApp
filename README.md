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

**Step 2 — the McKay quiver.**
- one gauge node per irrep `R_i` of Γ, with gauge group `U(N·dim R_i)`;
- bifundamental chirals: the number of arrows `i→j` is the multiplicity
  `a_ij = ⟨χ_Q·χ_i, χ_j⟩`, where `Q` is the defining 3-dim rep;
- the cubic Leigh–Strassler superpotential descends from
  `W_{N=4}=Tr Φ¹[Φ²,Φ³]`; its terms are the Γ-invariant closed 3-loops,
  counted by `Tr(A³)`.

The irreducible characters are obtained numerically from the class algebra
(Burnside / Dixon), so the quiver is built for any finite matrix group.

**Step 3 — the conformal manifold.** By Green–Komargodski–Seiberg–Tachikawa–
Wecht the conformal manifold is the space of marginal couplings modulo the
complexified global symmetry. For these orbifold theories the count collapses
to a closed form controlled by how the three holomorphic coordinates of C³ are
left fixed:

```
fix_Q(g) = number of unit eigenvalues of g on C³
dim_C M_conf = ( Σ_{g∈Γ} fix_Q(g) ) − 1            (faithful Γ, connected quiver)
```

By Burnside this equals the sum, over the holomorphic directions, of the number
of group elements fixing that direction. For a diagonal cyclic action it is

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

## Scope / caveats

- The closed form is for a **faithful** Γ (a genuine SU(3) subgroup), which
  gives a connected quiver. Non-faithful actions — where the abstract group has
  a kernel on C³, so the matrix image is a proper quotient — change the overall
  `−1` (the single decoupled U(1)) and are rejected by the cyclic constructor
  rather than answered incorrectly.
- The McKay character table is computed numerically; group orders are capped
  (`closure(..., max_order=5000)`) as a finiteness guard.
```
