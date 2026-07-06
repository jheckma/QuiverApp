# Resume handoff — 2026-07-06 session (dimer face clicks, state-sync fix, encoding repair)

Repo: `C:\Users\jheckman\Desktop\QuiverApp\` (moved out of Blobfish; GitHub `jheckma/QuiverApp`)
Branch: **`toric-resolutions`**, merged to **`main`** at the end of this session (user request).
Tests at commit time: **334 passing** (`python -m pytest tests/ -q`).
Run the app: `python -m conformalmanifold.webapp --port 8000` (NOTE: `--port` flag; a bare
positional port is silently ignored and you get 8000 anyway).

## What shipped today (this session — Claude)

1. **Seiberg duality by clicking dimer faces** (`HANDOFF_dimer_face_click.md` has full detail):
   - `inverse.py::face_polygons(t)` — per-gauge-face boundary polygons in the exact
     universal-cover coordinates `drawTiling` uses; exported as `tiling.faces`.
   - `index.html` — `.dface` overlay: every face shows its node id; SQUARE faces glow
     (dashed warm) and are clickable, feeding the same `_dualPath`/`&dualize=` machinery
     as quiver-node clicks. 10 new tests.
2. **pts= hash-sync fix** — `#pts=` was read-only-at-load while `&dualize=` updated live,
   so switching diagrams then clicking faces built chimera URLs (old diagram + new path)
   that reproduced as garbage quivers on reload. `buildToric` now rewrites `pts=` on every
   successful build. This was the root cause of the user's "the quiver looks wrong" report —
   the duality engine itself audited clean (`_qa_faceaudit.py`, all seeds, all phases).
3. **Mojibake repair** — the concurrent Codex session PS-rewrote files: `index.html`
   double-encoded + CRLF'd, `api.py` once. Repaired via `_qa_fixenc.py` (iterative
   cp1252 peel); non-ASCII inventory verified identical to pre-damage HEAD.

## Concurrent-agent situation (IMPORTANT for next session)

- **Codex** built a **BPS-quiver tab** in this working tree: `conformalmanifold/bps.py`,
  `tests/test_bps.py` (8 tests), `/api/bps_quiver` endpoint, `#bps` tab in index.html,
  `HANDOFF_bps_quiver.md`. Committed this session (separate commit, credited to Codex).
  Its one shared-code edit (`drawQuiver` optional `labels` opt) is backward-compatible.
- **RULE for any agent in this repo: never edit files with PowerShell
  Get-Content/Set-Content/-replace** (PS5.1 double-encodes UTF-8, silent mojibake) —
  Edit-tool or Python I/O only. If mojibake reappears, run `python _qa_fixenc.py`
  (dry run), then `--apply`.
- **origin/main also advanced in parallel** (AdS6/CFT5 defect groups + anomalies,
  in-app bibliography, `superconformal/` package, `add-4d-superconformal` branch) —
  other agents push to this repo; always `git fetch` and check both branch tips
  before assuming local state.

## QA drivers in repo root (all headless-Chrome CDP, server must be on :8000)

- `_qa_dimerclick.py` — face-click feature end-to-end (F0 flip/involution, C3 negative,
  deep-link restore).
- `_qa_faceaudit.py` — duality-correctness audit (reversal/chirality/anomaly/involution)
  + click-polygon containment; backend-only, fast.
- `_qa_ptssync.py` — pts= hash sync replay (preset switch → face click → cold reload).
- `_qa_fixenc.py` — mojibake diagnose/repair.
- Screenshots land in `_qa/` (untracked).

## Open / next steps

- Optional: hover cross-highlight dimer face ↔ quiver node (plumbing exists via data-qnode).
- Codex's BPS tab: only smoke-tested (tests pass, tab loads); deeper physics review pending.
- Inherited niggles in `HANDOFF_seiberg_duality.md` + `HANDOFF_toric_blowup.md`.
- Cosmetic: schematic-layout face polygons can be non-convex on mutated phases.
- Pre-existing engine gate: phase enumeration capped at 2·area ≤ 6.

## Unrelated same-day context (other sessions/repos)

- Taiga QA response written + compiled in Blobfish: `blobfish_qa_response_july2026.tex/pdf`
  (byline "Team Blobfish"); ZMxZN answer-collision findings accepted, remediation plan in doc.
- Blobfish `zmzn-quiver/batch.jsonl` + `batch2.jsonl` are FROZEN (answer collisions) —
  regenerate with unique canonical answers before ANY submission.
