# Handoff — QuiverApp: SCFT data + inverse algorithm (2026-06-21)

> **Status update (2026-06-22):** everything below is now **committed** on
> `toric-resolutions`, plus the follow-up polish: orbifold-tab crash fix,
> Kasteleyn Newton-polygon check, per-field R-charges (Butti–Zaffaroni), and a
> planar brane-tiling rendering (convex faces). Suite now **248 passing**.
> This doc is kept for historical context; the "open" items below are done.

Repo: `C:\Users\jheckman\Desktop\Blobfish\conformal-manifold`
Branch: **`toric-resolutions`** (the active QuiverApp branch).

Run the app: `python -m conformalmanifold.webapp --port 8000` → http://localhost:8000
Tests: `python -m pytest tests/ -q` → **248 passing**.

---

## What got done today

### 1. Superconformal (SCFT) data on the display — DONE, both tabs
New module **`conformalmanifold/scft.py`** + **`tests/test_scft.py`** (12 tests).
Wired into `api.py` (`summarize` → `"scft"` key; `summarize_toric_web` → `"scft"`)
and rendered in `static/index.html` (new "Superconformal data" cards on both tabs;
JS `renderScft` for orbifold, `renderToricScft` for toric).

- **Orbifold C³/Γ tab** (exact, closed form): gauge group ∏ᵢU(Nᵢ), every chiral
  R = 2/3, R_W = 2; 't Hooft anomalies Tr R (=0) and Tr R³; central charges
  **a = c = (|Γ|/4)N²**, a/a_{N=4} = |Γ|. Key identity: McKay relation
  Σⱼ a_ij d_j = 3 dᵢ ⇒ Tr R = 0 ⇒ a = c exactly.
- **Toric (p,q)-web tab**: Martelli–Sparks–Yau volume (Z-)minimisation over the
  Reeb vector b=(3,b₂,b₃). a = c = π³N²/(4 Vol); reports a/N² = 3/(4·g_min)
  (π-free), the Reeb vector, Vol(X₅)/Vol(S⁵), and per-corner R-charges (Σ=2).
  Damped Newton on the strictly-convex volume; starts at b₂,b₃ = 3·centroid.
- **Verified**: C³ → ¼, conifold → 27/64, dP0 → ¾; dP0 (toric) cross-checks the
  orbifold route C³/Z₃ exactly. All in `test_scft.py`.

### 2. Inverse algorithm (toric diagram → quiver) — DONE, general polygons
New module **`conformalmanifold/inverse.py`** + **`tests/test_inverse.py`**
(23 tests). Wired into `api.py` (`summarize_toric_web` → `"inverse_quiver"` key,
`max_gauge=40`) and rendered in `static/index.html` (new "Reconstructed quiver —
inverse algorithm" card; JS `renderInverse` + `drawTiling`).

Method = **Gulotta "properly ordered dimer"** (arXiv:0807.3012), zig-zag / medial-
graph picture:
- zig-zag windings = primitive outward edge-normals of the toric diagram (= the
  (p,q) legs), one per primitive boundary segment;
- draw them as straight geodesics on T² = R²/Z²; the arrangement is the **medial
  graph** of the tiling. Crossings (|det(w_k,w_l)| per pair) = fields; faces with
  all-forward darts = white vertices, all-backward = black, mixed = gauge faces;
- combinatorial map: `alpha` (arc involution, from ordering each path's crossings
  by parameter) + `sigma` (rotation at a crossing, by direction angle); faces are
  orbits of `phi = sigma∘alpha`. **Face type from dart signs** (the clean
  criterion that made it work: all +1 → white, all −1 → black, mixed → gauge).
- "Properly ordered" placement found by a **deterministic seeded base-point
  search** (`inverse_quiver`, LCG, up to 400 attempts) that rejects degenerate
  placements and accepts the first consistent tiling (counts, anomaly-free,
  two-term W). Reproducible across runs (test_inverse_determinism).

Returns BraneTiling: quiver adjacency, fields (with src/tgt + white/black vertex),
two-term superpotential (white=+, black=−), tiling node/edge positions on T², and
a `checks` dict (gauge=2·area, fields=Σ|det|, white=black, anomaly-free, toric W,
Euler V−E+F=0).

**Verified** (test_inverse.py): C³(1,3), conifold(2,4), dP0(3,9), F0(4,8),
dP1(4,10), dP2(5,11), dP3(6,12), C³/(Z₂×Z₂)(4,12), SPP(3,7) — all textbook
counts; orbifold cases (dP0, Z₂×Z₂) cross-check the McKay quiver node degrees.

---

## Files touched
- NEW `conformalmanifold/scft.py`, `tests/test_scft.py`
- NEW `conformalmanifold/inverse.py`, `tests/test_inverse.py`
- `conformalmanifold/api.py` — 3 new imports; `"scft"` in `summarize`;
  `"scft"` + `"inverse_quiver"` in `summarize_toric_web`
- `conformalmanifold/static/index.html` — 3 new cards + JS (`renderScft`,
  `renderToricScft`, `renderInverse`, `drawTiling`, `fmtCC`)

Prototype scratch files in the system temp dir (`proto*.py`) are throwaway.

---

## Open / next steps
1. **Visual QA in a browser** — backend is fully tested and JS balances clean
   (`node` isn't installed here, so the script wasn't run through a real parser).
   Load the page and eyeball the three new cards, especially the dimer rendering
   (`drawTiling`: 2×2 cell tiling, nearest-image edges; node positions are the
   circular-mean of their crossings, may look slightly off near cell seams).
2. **Commit** — not yet committed. Suggested: one commit for SCFT, one for the
   inverse algorithm.
3. Possible polish: Kasteleyn-determinant Newton-polygon check as an explicit
   extra test (currently the geometry is self-certified because windings come
   from the input). Per-field (not just per-corner) R-charges on the toric tab.
   Isoradial node placement for a prettier tiling.
4. Known prior issue (separate, from memory): the older centroid-junction dual-web
   rendering violated local charge conservation; the tropical/Legendre dual in
   `resolution.py` already fixed that — not related to today's work.
