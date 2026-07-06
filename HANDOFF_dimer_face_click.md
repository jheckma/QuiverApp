# Handoff — Seiberg duality by clicking dimer faces

Repo: `C:\Users\jheckman\Desktop\QuiverApp\`  ·  Branch: **`toric-resolutions`**
Tests: `python -m pytest tests/ -q` → **326 passing** (316 before this feature + 10 new).
Status: **DONE, verified headless, UNCOMMITTED** (user has not asked to commit).

## What shipped (2026-07-06)

Square gauge faces of the drawn brane tiling are now clickable Seiberg-duality
(urban renewal) moves — the dimer display is an input surface, not just an
output, on equal footing with the quiver-node clicks that already existed.
Every face also carries a faint node-id label, so the face ↔ quiver-node
correspondence is finally visible on screen.

### Backend — `conformalmanifold/inverse.py`

- **`face_polygons(t)`** (right after `square_gauge_faces`): per gauge face,
  returns `{"node": displayed_id, "square": bool, "poly": [[x,y],...]}`,
  sorted by node id. The polygon walks the face's dart orbit (the same
  `_trace_faces(dimer, _FACE_HAND)` trace `square_gauge_faces`/`urban_renewal`
  use), taking each dart's START vertex in the **same universal-cover
  coordinates `drawTiling` draws with**: white vertex at `white_glob[w] + T`,
  black at `black_glob[b] − h_e + T`. The per-edge integer translate `T`
  propagates by matching the shared vertex of consecutive darts:
  - across a shared **black** vertex (prev dart had s=+1): `T += h_next − h_prev`
  - across a shared **white** vertex: `T` unchanged.

  The walk closes because the face-cocycle condition (Σ signed homologies
  around a gauge face = 0 — the same rows `solve_homology` enforces)
  telescopes exactly to the accumulated shift. Works on mutated phases too
  (schematic `_spanning_layout` embeddings): the polygon is whatever the
  drawing shows, which is what a click target must be.
- `_tiling_json` now emits `tiling.faces = face_polygons(t)` alongside the
  existing `edges`/`edge_h`. The `square_faces` top-level field is unchanged
  (tests and the quiver-node glow depend on it).

### Frontend — `conformalmanifold/static/index.html`

- CSS: `.dface` (transparent polygon + faint centroid label, labels are
  `pointer-events:none`) and `.dface.sq` (dashed warm outline, warm bold
  label, hover fill, cursor pointer) — visual language matches `.node.sq`.
- `drawTiling`: a face-overlay pass between the edge pass and the vertex-dot
  pass. Each face polygon is replicated over the same lattice-translate
  window as the edges (drawn when its centroid passes `inview`). Square faces
  get `data-qnode="<node>"` on the polygon; non-square faces are label-only.
- Click wiring: `$("iv-tiling").addEventListener("click", _quiverNodeClick)` —
  the polygon's `data-qnode` rides the exact same path as a quiver-node click
  (`_dualPath` push → `syncDualHash` → `buildToric` → `&dualize=` deep-link,
  path badge, reset button). No new state machine.
- Help text in the `iv-dualize` row now says "click any quiver node — or a
  glowing square face of the brane tiling — to dualize it."

### Tests — `tests/test_inverse.py` (10 new)

- `test_face_polygons_cover_close_and_flag_squares` (parametrized over all 9
  CASES): one polygon per gauge node in displayed order; square flags ==
  `square_gauge_faces(t)`; sizes even, ≥4, summing to 2E (every dart is one
  corner); translate closure around every face orbit **including the
  wraparound transition** (which `face_polygons` itself never emits — that's
  the independent check).
- `test_face_polygons_on_mutated_phase`: F0 phase II via
  `inverse_phases_json` — faces still emitted on the urban-renewal product,
  sizes {4,4,8,8}, squares == `square_faces`.

## Verification (headless Chrome, `_qa_dimerclick.py`)

Same CDP pattern as `_qa_driver.py`; screenshots in `_qa/dc_0*.png`. Start the
server first: `python -m conformalmanifold.webapp --port 8000` — **the port is
a `--port` flag; a bare positional port is silently ignored** (bit me: server
came up on 8000 while the driver aimed at 8011, everything read as null).

- F0 seed: 4 labelled faces (49 replicas), all square/clickable, 8 fields.
- Click face 0 **in the tiling** → 12-field phase II, badge `seed → n0`,
  hash `&dualize=0`, squares now {0,2} (heritage labels), quiver glow matches.
- Click face 0 again → involution back to 8 fields (`seed → n0 → n0`).
- C3: hexagonal faces labelled, zero clickable (no square faces).
- Deep-link `#pts=...&dualize=0` cold-load restores the dualized state.
- Zero console errors/exceptions across all runs.

## Notes / conventions for next time

- **Repo location**: moved 2026-07-06 to `Desktop\QuiverApp` (was
  `Blobfish\conformal-manifold`). Old handoffs cite the old path.
- The polygon coordinate convention is COUPLED to `drawTiling`'s
  (`white_g` / `black_g − edge_h`, replicated over Z²). If the drawing
  convention ever changes, `face_polygons` must change with it — the
  closure test will not catch a uniform convention drift, only the
  screenshots will.
- Non-square faces are deliberately not clickable in the dimer: a non-square
  move leaves the toric regime, and the tiling it would be clicked on ceases
  to exist. The quiver nodes remain the input surface for those (unchanged).
- `_qa_dimerclick.py` is a throwaway QA driver in the repo root like the
  other `_qa_*.py` scripts; `_qa/` output stays untracked.

## Follow-up fix (same day): pts= hash desync — the "wrong quiver" report

User report: "the quiver after clicking a face looks wrong." Root cause was NOT
the duality (audited clean, see below) but a pre-existing deep-link bug the new
face clicks made easy to hit: **`#pts=` was only read at page load and never
rewritten when the diagram changed** (dots/presets reset `_dualPath` but not the
pts token), while `&dualize=`/`&phase=` update live. Flow that bit: window
opened on `#pts=F0` → user switched to dP3 (hexagon) → clicked 3 faces → hash
became `pts=F0&dualize=2.1.2` — a chimera. The DISPLAY was the correct 14-field
dP3 dual (squares {0,1,2,3,5}, matches the screenshot), but any reload applied
the dP3 path to F0 → a 4-node ranks-(1,3,5,1) non-toric quiver that looks like
a botched duality.

Fix: in `buildToric` stage-1 success, the hash rebuild now drops+unshifts a
fresh `pts=` token (same place `active=` is synced). Verified headless
(`_qa_ptssync.py`): preset click syncs pts; face click appends dualize; cold
reload reproduces the identical state; dot edits sync too. Known edge: builds
that fail the hull<3 guard (tEmpty) don't sync — harmless, no meaningful state.

Audits written for the report (keep, they're fast):
- `_qa_faceaudit.py` — for every square face of all 9 seeds: reversal at the
  dualized node, spectator chirality == meson composition, anomaly freedom,
  involution; plus label/click-polygon containment (each face centroid inside
  exactly its own polygon) for every phase AND every post-click tiling. All clean.

## Second follow-up (same day): repo-wide mojibake from the Codex session

A concurrent **Codex** session building the BPS-quiver tab (bps.py, test_bps.py,
api.py/webapp.py wiring, index.html tab — its own handoff: HANDOFF_bps_quiver.md)
re-saved files via PowerShell 5.1, double-encoding UTF-8: `index.html` was hit
TWICE ('✓'→'âœ“'→'Ã¢Å“â€œ') and flipped to CRLF; `api.py` once. Repaired with
`_qa_fixenc.py` (iterative cp1252-roundtrip peel, passthrough bytes 81 8D 8F 90
9D, only runs that repair SHORTER; CRLF→LF restored). Verified: non-ASCII run
inventory now identical to HEAD in both directions; page renders ✓/—/⊕
correctly; 334 tests pass (326 + 8 BPS); face-click QA green, 0 console errors.
Codex's one out-of-scope edit — `drawQuiver` gained an optional `labels` opt —
is backward-compatible (defaults to node index), left in place.
**If Codex keeps working here, tell it: no PS Get-Content/Set-Content/-replace
on repo files — Edit-tool/Python only** (same rule we already follow), else the
mojibake returns.

## Open (optional) items

- Inherited niggles from `HANDOFF_seiberg_duality.md` +
  `HANDOFF_toric_blowup.md` (unchanged by this feature).
- Cosmetic: on busy mutated-phase layouts the schematic face polygons can be
  non-convex; hover-fill looks slightly ragged there. Purely visual.
- Possible follow-up: hover on a dimer face could highlight the matching
  quiver node (and vice versa) — the `data-qnode` plumbing already supports it.
