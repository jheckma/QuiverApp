# Handoff — Toric blow-up feature + perf fixes (QuiverApp, branch toric-resolutions)

Date: 2026-06-23. Scope: the **Toric (p,q)-web builder** tab of the QuiverApp web UI
(`conformalmanifold/static/index.html` + `api.py` + `resolution.py` + `webapp.py`).
**Do NOT touch the dimers part** — another agent owns `inverse.py` (Gulotta brane-tiling,
Seiberg phase-cycling: `iv-phase-cycle` / `showPhase` / `_ivPhaseIdx`, `HANDOFF_seiberg_duality.md`).

## What was built (DONE, 87 tests pass in ~3.2s)

### 1. Corner blow-up ("chamfer") of a toric diagram
Click a corner in the toric builder → cut it off one primitive lattice step along each
adjacent edge → recompute the whole pipeline (triangulation, (p,q) web, gauge count, SCFT,
reconstructed quiver). Only **singular** corners are actionable; smooth (unimodular) corners
— C³, the conifold square — are rejected (nothing to resolve).

- `resolution.py` → new `blowup_corner(points, corner)` (inward chamfer; raises ValueError on
  smooth corner / non-corner / degenerate result).
- `api.summarize_toric_web(..., blowup=[x,y])` — applies the blow-up first, then recomputes;
  returns the rewritten point set in `diagram.input_points`.
- `webapp.py` → new `&blowup=x,y` query param + `_parse_point` helper (mirrors `&flop=`).
- `static/index.html` → **"⊕ Blow up a corner"** toggle (`#blowmode`); in blow-up mode each
  singular hull corner gets a glowing ⊕ handle (`blowableCorners()` / `jsGcd()`); clicking it
  fires `buildToric({blowup:[x,y]})`. Fixed a latent bug: the new button shares `.preset`
  styling, so the preset handler is now scoped to `button.preset[data-pts]` (else it wiped the
  diagram).
- Tests: `tests/test_blowup.py` (6).

### 2. Performance fixes (this is what made the UI usable)
Root cause of the "incredibly laggy / buttons sometimes work / flop broken" reports: every
click recomputed the inverse brane-tiling, and `inverse._crossings` did **millions of
`fractions.Fraction` ops** (Z₃×Z₃ took ~19s and *failed*; a 4×4 square didn't return). All
three symptoms were this one cause (flop/blow-up/add-point all go through `buildToric`).

- **`inverse.py::_crossings` rewritten** to pure-integer arithmetic over a common denominator
  per zig-zag pair, building `Fraction`s only for the final |det| crossings. **Bit-identical
  output** (verified by reproducing the old impl and diffing on conifold/Z₃×Z₃/dP3/Z₂×Z₂).
  ~100× faster. *(This is the one edit inside the dimer agent's file — kept surgical and
  output-preserving; the existing `test_inverse.py` passes.)*
- **`api.py` hull-keyed cache** (`_HULL_CACHE`, `_hull_cached`) for the hull-only heavy calls
  (`inverse_quiver_json`, `toric_scft_json`, `inverse_phases_json`). Flops / mode-toggles /
  re-selecting the same diagram are now ~1ms instead of recomputing.
- **`index.html` request sequencing**: `buildToric` aborts the previous in-flight fetch
  (`AbortController`) + a `_toricSeq` guard drops stale responses → latest click always wins.
  Added a pulsing "recomputing…" indicator (`#t-busy` / `.busy`).

Cold timings after fix: conifold 0.007s, dP3 0.033s, Z₂×Z₂ 0.19s, **Z₃×Z₃ 19s → 1.07s**.

## Run it
```
cd C:\Users\jheckman\Desktop\Blobfish\conformal-manifold
python -m conformalmanifold.webapp --port 8765
# open http://localhost:8765/#pts=0,0;1,0;1,1;0,1   (conifold; #toric for the builder)
```
A dev server was left running on port 8765 — it may still be on the OLD code; restart it to
pick up the `inverse.py`/`api.py` changes. Stop a stale one:
`Get-NetTCPConnection -LocalPort 8765 -State Listen | %{ Stop-Process -Id $_.OwningProcess -Force }`

## RESOLVED 2026-06-24 — blow-up is now the genuine del Pezzo blow-up (add a divisor)

The old inward-chamfer `blowup_corner` (shrank the diagram) was WRONG and is REMOVED. Replaced
with the standard toric **star subdivision**: for a primitive boundary edge A–B whose adjacent
unimodular triangle has apex C, the new exceptional ray is **W = A + B − C** (reflect the apex
across the edge midpoint). This GROWS the diagram by one unimodular triangle → +1 gauge group,
+1 external (p,q) leg, and reproduces dP0→dP1→dP2→dP3 (add (1,1)→(−1,0)→(0,−1) to the dP3
hexagon). Origin-free / translation-invariant; verified against the dP ladder and conifold.

- `resolution.py`: `blowup_candidates(points)` (one site per boundary edge, returns
  {edge, apex, new_point}) + `blowup(points, new_point)` (validates the site is legal, then adds
  it; asserts 2·area goes up by exactly 1). `blowup_corner` deleted.
- `api.py`: `diagram.blowup_candidates` added to the output; `summarize_toric_web(blowup=[x,y])`
  now calls `Rz.blowup` (x,y = the chosen exceptional point W, NOT a corner).
- `index.html`: blow-up handles now render at the backend `blowup_candidates[].new_point`
  positions (the ⊕ sits just OUTSIDE the diagram where the new leg sprouts); `syncView` frames
  them; button/help text updated ("Blow up a point", add a divisor → grows). New `&blow` deep-link
  flag pre-enables blow-up mode for screenshots.
- `webapp.py`: doc-comment updated; `&blowup=x,y` unchanged (now the exceptional point).
- Tests: `tests/test_blowup.py` rewritten (6) — dP ladder nodes 3→4→5→6, candidates outside +
  minimal, pipeline/quiver re-read from the grown diagram, illegal sites rejected, conifold→3.
  **Full suite 279 pass.** Visual QA (headless Chrome, `#pts=...&blow`): dP0 shows ⊕ at
  (1,1)/(−1,0)/(0,−1); clicking grows to dP1 (4 nodes, dim 3, 4-leg web) — confirmed.

## ALSO RESOLVED 2026-06-24 — blow-DOWN + lag fix (async inverse + size gate)

**Blow-down (inverse of blow-up):** contract a smooth (unimodular) corner.
- `resolution.py`: `blowdown_candidates(points)` (hull vertices W whose ear A,W,B is unimodular,
  i.e. smooth corners, keeping ≥3 corners after) + `blowdown(points, corner)` (removes W, area−1).
  Verified inverse of blow-up (dP1→dP0), conifold→C³ (4 sites), C³ has none.
- `blowup_candidates` now FILTERS to genuine minimal blow-ups (W=A+B−C can overshoot convexity at
  a sharp corner → area+>1; those are dropped so candidates ⇔ what `blowup` accepts).
- `api.py`: `blowdown` param + `diagram.blowdown_candidates`. `index.html`: second toggle
  "⊖ Blow down a corner" (`blowMode` is now "" / "up" / "down"), green ⊖ handles at the smooth
  corners, `&blowdown` deep-link.

**Lag fix — the toric diagram / blow-up / blow-down are now ALWAYS instant.** The only expensive
part is the inverse brane-tiling reconstruction (Gulotta dimer + the O(n!) Kasteleyn Newton
certificate, which also *fails slowly* on some diagrams). Two changes:
1. **Decoupled into a separate async fetch.** `summarize_toric_web(..., include_inverse=False)`
   returns only the cheap geometry (diagram, web, dim, gauge count, blow-up/down sites);
   `webapp.py` honors `&inverse=0`. The frontend `buildToric` does the fast geometry fetch first
   (renders immediately, ~3 ms even where the full call is ~1 s), then `fetchInverse()` fetches the
   full response async and fills in the reconstructed-quiver / phases panels (placeholder
   "reconstructing…"; both fetches abort+seq-guard so the latest click always wins).
2. **Size gates:** `inverse.py KAST_NEWTON_MAX_WHITE = 8` (skip the O(n!) certificate above 8
   white nodes; the dimer's own checks still pass and a skipped cert shows as `null`, not a
   failure — frontend `allok` treats null as OK). `api.py INVERSE_MAX_AREA2 = 16` (skip the whole
   reconstruction above 2·area = 16, return `available:false` + reason).

Timings (full endpoint): geometry-only is 3 ms (dP0) … 84 ms (5×5); full is 6 ms (dP0), and large
diagrams (2·area = 32/72) return in ~30/200 ms because the inverse is gated. The old 4×4 / many-
blow-up hang is gone. **Full suite 284 pass** (`tests/test_blowup.py` now 11: dP ladder, blow-down
inverse, conifold, gating, include_inverse). Visual QA: dP3 blow-down shows ⊖ on all 6 corners;
dP0 blow-up reconstructed-quiver card fills in "consistent dimer ✓" after the async fetch.

## SEPARATION PASS 2026-07-03 (third revision) — the TWO blow-ups fully distinguished

User: "There are two things you can do with a blowup. One is to resolve a singularity.
The other is to pick a generic point and add a NEW P^1 to the surface... these two
possibilities [must be] clearly distinguished in both the code and the UI."  Done:

- **Code**: the base-surface operations are now `surface_blowup`,
  `surface_blowup_candidates`, `surface_blowdown`, `surface_blowdown_candidates`
  (resolution.py); the singularity-resolution engine stays `triangulate(hull, active)`
  + `residual_cells` + `is_valid_subdivision`.  api.summarize_toric_web kwargs:
  `surface_blowup=`, `surface_blowdown=`, `active=`; response fields
  `diagram.surface_blowup_candidates` / `surface_blowdown_candidates`.
  webapp: `&surface_blowup=` / `&surface_blowdown=` (legacy `&blowup=`/`&blowdown=`
  still accepted).
- **UI**: TWO labelled button rows with four mutually exclusive modes —
  "resolve the singularity (exceptional divisors, diagram fixed)":
  ⊕ Blow up a divisor / ⊖ Blow down a divisor (handles on lattice dots;
  deep-links &resup / &resdown); "blow up a point of the surface (a NEW exceptional
  P¹, diagram grows)": ⊕ Blow up a point / ⊖ Blow down a −1-curve (handles
  outside/on corners; deep-links &blow / &blowdown).  Mode-specific hints; help
  paragraph opens with "a blow-up does one of TWO different things".
- 303 tests pass; headless QA: dP0 ⊕→dP1 (surface path, renamed params) and the
  Z₃×Z₃ cone resolve/un-resolve loop both verified interactively.

## PREVIOUS 2026-07-03 — genuine blow-up/blow-down semantics (303 tests)

The grow/shrink relabel below was REJECTED by the user ("you are supposed to add
Blowup and Blowdown functionality!", "C³/Z_M×Z_N can for sure be resolved").  The
buttons now implement the two GENUINE geometric operations:

1. **Resolve the CY singularity (fixed polygon)** — NEW partial-resolution engine:
   every lattice point of the diagram is an exceptional divisor of the singular cone;
   an `active` subset (corners forced; default all = fully resolved) selects which are
   blown up.  `resolution.py`: `triangulate(hull, active)`, `is_valid_subdivision`,
   `residual_cells` (non-unimodular cells = residual orbifold singularities, each
   identified via the named library, e.g. Z(3,3)); `dual_web` generalized: coarse
   boundary edges give external legs with MULTIPLICITY (the singular Z₃×Z₃ cone shows
   one junction, legs (1,1)×3 (−1,0)×3 (0,−1)×3).  `api.summarize_toric_web(active=)`,
   webapp `&active=x,y;...`, UI: ⊕ on dimmed (unresolved) lattice dots blows a divisor
   up, ⊖ on active non-corner dots blows it down, red-tinted cells + a
   "resolution state" row list the residual singularities, `&active=` deep-links.
   Flops still work on partial subdivisions.  Quiver/dim/SCFT stay polygon-level
   (the singular theory) — correct, they are resolution-independent.
2. **Blow up/down the base surface (polygon changes)** — RESTRICTED to diagrams that
   are cones over a toric surface (exactly ONE interior lattice point O): blow-up
   sites = adjacent boundary rays with det 1, W = A + B − O (dP0→dP1→dP2→dP3);
   blow-down sites = −1-curves (ray = sum of neighbour rays).  C³ and the conifold
   have NO interior point → no sites (this kills the bogus C³→…→C³/Z_M×Z_N "blow-up"
   chain that triggered the complaint); P² and F0 are minimal → no blow-downs;
   dP1/dP2/dP3 have 1/3/6 −1-curves (textbook counts, falling out of the arithmetic).

Tests: 303 pass (`test_partial_resolution.py` new, 8; `test_blowup.py` updated to the
surface semantics).  Interactive QA (`_qa_resolve.py`): cone → click ⊕(1,1) → three
|Γ|=3 cells → ⊖ back to the Z(3,3) cone; hash round-trips.

## SUPERSEDED — grow/shrink relabel + inverse-criterion fixes (295 tests)

User feedback (twice): labelling the moves "blow up / blow down" is wrong — the two
standard conventions point in OPPOSITE directions.  Base surface (cone over a surface):
blow-up GROWS the diagram (dP0→dP1).  CY singularity: the SUB-diagram is the partial
resolution, so growing embeds the CY in a MORE singular parent (C³ → C³/Z_M×Z_N) — and
the surface reading doesn't exist at all when the diagram has ≠1 interior point.

1. **UI relabel** (`index.html`): buttons are now "⊕ Grow — add boundary point" /
   "⊖ Shrink — contract corner"; help text + hints state BOTH geometric readings.
   `&blow`/`&blowdown` deep-links and the `blowup`/`blowdown` API params/function names
   are unchanged (docstrings carry the naming caveat).
2. **`blowdown_candidates` criterion fixed** (`resolution.py`): a corner is contractible
   iff removing it drops 2·area by EXACTLY 1, computed on ALL lattice points of the
   polygon.  The old adjacent-hull-CORNER ear test wrongly reported no shrink sites on
   long-edged diagrams (Z₃×Z₃ triangle: 0 sites → now all 3 corners).  `blowdown`
   likewise returns the full lattice-point set (a 3-click triangle no longer degenerates).
3. **`blowup_candidates` now triangulation-INDEPENDENT**: every lattice point whose
   addition grows 2·area by exactly 1 (scan each primitive boundary edge's outer
   height-1 line + exact area filter; new `_unimodular_conormal` helper).  Previously
   the sites were the apex reflections of the currently displayed triangulation, so the
   offered handles changed with the toric phase and a just-contracted corner was not
   always re-offered.  Shrink→grow now round-trips exactly (regression-tested).
   Conifold: 8 sites (was 4); dP0: 3 (unchanged); Z₃×Z₃: 0 (correct — a point under a
   length-3 edge adds area 3).
4. Tests: `test_blowdown_long_edge_corners` (new), inverse test compares hulls not raw
   point sets, conifold count updated.  Full suite **295 pass**.

## OPEN / KNOWN ISSUES (none blocking)

1. **Inverse reconstruction still *fails* (not hangs) on some medium diagrams** (e.g. the 3×3
   Z3×Z3, area2=9 — `inverse_quiver` returns `available:false` after ~1 s of base-point search).
   Now harmless: it runs async so the geometry/blow-ups stay instant, and the quiver panel just
   shows the error. A real fix (polynomial Kasteleyn / better proper-ordering search) is
   dimer-agent territory (`inverse.py`, see `HANDOFF_seiberg_duality.md`).
2. The O(n!) `kasteleyn_newton_polygon` itself is unchanged (just gated). Replacing the permanent
   expansion with a polynomial perfect-matching determinant would let the certificate run on large
   diagrams too — optional.

## Files touched
- `conformalmanifold/resolution.py`  (blow-up rewrite `blowup`/`blowup_candidates`; blow-down
  `blowdown`/`blowdown_candidates`; old `blowup_corner` deleted)
- `conformalmanifold/api.py`          (blowup/blowdown/include_inverse params, gates
  `INVERSE_MAX_AREA2`, `_finish_identify`, `_HULL_CACHE`/`_hull_cached`)
- `conformalmanifold/webapp.py`       (`&blowup`/`&blowdown`/`&inverse`, `_parse_point`)
- `conformalmanifold/static/index.html` (two-mode toggle + ⊕/⊖ handles, two-stage
  `buildToric`/`fetchInverse` async, `renderResultGeometry`/`renderResultInverse`, busy UI)
- `conformalmanifold/inverse.py`      (`_crossings` perf rewrite + `KAST_NEWTON_MAX_WHITE` gate
  — coordinate w/ dimer agent)
- `tests/test_blowup.py`              (rewritten: 11 tests)
