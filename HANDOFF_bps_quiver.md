# Handoff: BPS Quiver Tab

Date: 2026-07-06

## What changed

Added a new **BPS quivers** tab to the existing single-page QuiverApp UI. The tab interprets an input quiver as a BPS quiver and reports stable objects for a chamber encoded by a mutation sequence.

Main files:

- `conformalmanifold/bps.py`
  - New finite-chamber BPS backend.
  - Parses exchange matrices or arrow lists.
  - Tracks c-vectors through quiver mutation.
  - Reports stable hypermultiplet charges for green mutations.
  - Detects green vs non-green sequences and maximal-green sequences.
  - Optionally computes central charges from phases and masses and checks phase ordering.
  - Preserves zero-based labels when input uses node `0`.
  - Rejects self-arrows, nonpositive multiplicities, arrow-list 2-cycles, invalid skew matrices, and nonpositive masses.

- `conformalmanifold/api.py`
  - Added `summarize_bps_quiver(...)` wrapper.

- `conformalmanifold/webapp.py`
  - Added `GET /api/bps_quiver?kind=...&quiver=...&sequence=...&phases=...&masses=...`.

- `conformalmanifold/static/index.html`
  - Added third tab: `BPS quivers`.
  - Added matrix/arrow input, mutation sequence input, optional phases/masses, examples, chamber status, quiver diagram, exchange/adjacency matrices, final exchange matrix, stable-object table, and mutation log.
  - Reuses existing `drawQuiver` / `drawAdj`; `drawQuiver` now accepts optional display labels.
  - Deep link `#bps` opens and computes the BPS tab.

- `tests/test_bps.py`
  - Covers A2 maximal-green sequence `1 2 1` with stable charges `[1,0]`, `[1,1]`, `[0,1]`.
  - Covers partial green and non-green warning semantics.
  - Covers arrow multiplicity parsing, zero-based labels, 2-cycle rejection, mass validation, and invalid skew matrices.

## Current API contract

A chamber is represented by a finite mutation sequence. For each green mutation, the positive c-vector is emitted as a stable hypermultiplet charge/dimension vector in that chamber.

Example:

```text
quiver:   1->2
sequence: 1 2 1
```

returns stable objects:

```text
gamma_1
gamma_1 + gamma_2
gamma_2
```

and marks the chamber as maximal green.

## Verification performed

- `python -m py_compile conformalmanifold\bps.py conformalmanifold\api.py conformalmanifold\webapp.py`
- `python -m pytest tests\test_bps.py tests\test_pipeline.py -q` -> 34 passed
- `python -m pytest -q` -> 334 passed
- Direct API smoke test through `conformalmanifold.api.summarize_bps_quiver(...)`
- Live HTTP smoke tests against `http://127.0.0.1:8017`:
  - `/` serves BPS tab markup and script
  - `/api/bps_quiver` returns A2 stable spectrum
  - invalid non-skew matrix returns HTTP 400 with a useful error
- Frontend script imported successfully with DOM/fetch stubs. The in-app browser runtime had no available browser instance, so screenshot/click QA could not be performed here.

## Running locally

Current dev server was started on:

```text
http://127.0.0.1:8017/#bps
```

Normal command:

```powershell
python -m conformalmanifold.webapp --port 8017
```

## Notes / next steps

- The current central-charge phase check is intentionally lightweight: it checks non-increasing phases on the chosen branch after charges are computed. It does not yet fully solve branch-cut/wall-crossing chamber validation.
- The BPS tab reports finite-chamber c-vector spectra. It does not yet implement a general stability solver from arbitrary central charges independent of a mutation sequence.
- Browser-level visual QA should be run when an in-app browser or Chrome surface is available.
- Existing unrelated dirty work was present before/alongside this task, notably Seiberg-duality/dimer-click work in `static/index.html`, `inverse.py`, `tests/test_inverse.py`, and QA/handoff files. Do not revert those while continuing BPS work.

## Geometry-source correction

The BPS tab was corrected to use quivers already built from geometry:

- `/api/bps_quiver?source=toric&pts=...&sequence=...` now calls `api.summarize_toric_web(...)`, selects the existing geometry-built quiver (`inverse_quiver` when available, otherwise named-library `quiver`), and runs BPS chamber analysis on that quiver.
- The BPS exchange matrix is `B = A - A^T`, where `A` is the geometry quiver adjacency. Self-arrows are ignored in `B`; reciprocal arrows are netted and reported as a warning.
- The BPS tab defaults to dP0 and offers dP1, F0, C3, conifold, custom points, and "current toric tab diagram". Preset buttons now compute immediately.
- Geometry quiver nodes are treated as 0-based, matching the toric quiver drawings.

Verification after this correction:

- `python -m pytest tests\test_bps.py -q` -> 10 passed
- `python -m pytest -q` -> 336 passed
- Fresh dev server started on `http://127.0.0.1:8018/#bps`
- HTTP smoke test for dP0 geometry returned adjacency `[[0,3,0],[0,0,3],[3,0,0]]` from the inverse brane-tiling quiver and stable charges `gamma_0`, `3 gamma_0 + gamma_1`, `gamma_2` for sequence `0 1 2`.

## 2026-07-06 update: toric diagram, SW curve, clickable BPS mutations

User requested that the BPS tab show the toric diagram, expose the underlying Seiberg-Witten curve, and use the same clickable-node interface for mutations.

Changes made:

- `conformalmanifold/api.py`
  - Added Newton-polynomial helpers for the toric diagram.
  - `summarize_bps_toric_quiver(...)` now includes `source.lattice_points` from the toric resolution and `source.sw_curve`, e.g. `H(x,y) = c[-1,-1] x^-1 y^-1 + c[0,0] + c[0,1] y + c[1,0] x = 0` for dP0.

- `conformalmanifold/static/index.html`
  - BPS tab now has a `Toric diagram and SW curve` panel.
  - `bpsDrawToric(...)` renders the toric polygon and lattice points from the same geometry payload used to build the BPS quiver.
  - `renderBpsGeometry(...)` fills the SW/Newton curve panel from `source.sw_curve`, with a frontend fallback if the backend field is absent.
  - The BPS quiver now calls `drawQuiver(..., {clickable:true})` and registers `handleBpsQuiverClick(...)` on `#bps-quiver-svg`.
  - Clicking a BPS quiver node appends that node id to `#bps-sequence` and recomputes the chamber.

- `tests/test_bps.py`
  - Added assertions that geometry-sourced BPS output carries lattice points and the SW curve.

Verification performed after this update:

- `python -m py_compile conformalmanifold\api.py conformalmanifold\bps.py conformalmanifold\webapp.py` -> passed
- `python -m pytest tests\test_bps.py -q` -> 10 passed
- `python -m pytest -q` -> 336 passed
- API smoke through `api.summarize_bps_toric_quiver(...)` returned the dP0 curve and lattice points.
- Live route smoke against `http://127.0.0.1:8019/api/bps_quiver?source=toric&pts=1,0;0,1;-1,-1&sequence=0%201%202` returned stable count `3`, lattice points, and the SW curve.
- Static frontend smoke confirmed `bps-toric-svg`, `bps-sw-curve`, `bpsDrawToric`, `renderBpsGeometry`, `handleBpsQuiverClick`, clickable `drawQuiver`, and the BPS quiver click listener are present.

Browser QA note:

- Browser plugin instructions were followed. The browser runtime loaded, but `agent.browsers.list()` returned `[]`, so no in-app browser/Chrome instance was available for screenshot/click QA from Codex.
- Node is also not installed (`node` command not found), so a JS DOM execution smoke could not be run. Verification used Python tests, HTTP checks, and static HTML checks instead.

Current fresh dev server for this update:

```text
http://127.0.0.1:8019/#bps
```

Next likely BPS work:

- Add a clear way to remove/reset the clicked mutation sequence without leaving the BPS tab.
- If richer SW data is desired, promote `source.sw_curve` from a symbolic Newton polynomial to a structured list of moduli/coefficients and known mass constraints for each named geometry.
- Run visual browser QA once a browser instance is actually registered with the Browser plugin.

