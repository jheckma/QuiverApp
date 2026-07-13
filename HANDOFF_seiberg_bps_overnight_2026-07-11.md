# Handoff — Seiberg-duality fix + BPS-quiver / spectral-network work (overnight 2026-07-11 → 07-12)

Repo: `C:\Users\jheckman\Desktop\QuiverApp`  ·  Branch: **`main`** (working tree, UNCOMMITTED)
Tests: `python -m pytest tests/ -q` → **372 passing** (was 356; +8 BPS, +8 inverse).
All work verified with the full suite AND a live headless-Chrome smoke test.

## 1. Seiberg duality: the "invalid dual phases when you dualize the wrong node" bug — FIXED

### Root cause (found by fuzzing multi-duality paths over many toric polygons)
Through a cascade of dualities the meson step legitimately **creates adjoints**
(self-loops) at gauge nodes.  Dualizing a node that *has* an adjoint is a
Kutasov-type duality, not ordinary SU(N) Seiberg duality.  `wmutation.mutate`
already refused it (raised `WMutationError`, so the tracked W dropped to None),
but `quiver_seiberg` had **no adjoint guard** and ran its naive meson/mass
arithmetic on the self-loop, emitting a nonsensical quiver with **negative arrow
multiplicities**.  This reproduced on C³/Z4, rect(2×3), L^{1,3,1}, the
trapezoid, etc. after 3–4 dualities.

Confirmed physically correct by the literature: SPP's adjoint node and F0
phase-II's 4N_c nodes behave exactly as the fix now dictates (Closset/FHH/FFHH).

### The fix (`conformalmanifold/inverse.py`, `quiver_seiberg`)
- Guard: if `A[k][k] != 0` raise `ValueError` ("node k has an adjoint … Kutasov-type").
- Flavor count now excludes the self-loop (defensive; moot given the guard).
The UI already turns a `ValueError` into a refused move.

### Frontend robustness (`static/index.html`) — 5 fixes from the frontend audit
1. **BUG 1 (was high):** a refused move used to leave the inverse card stuck on an
   error out of sync with the still-valid path.  Now it pops the failing move,
   flashes the reason, and **re-renders the last valid state** (`renderResultInverse`).
2. Top quiver no longer flickers back to the undualized seed mid-path
   (skip the stage-1 library redraw while `_dualPath` is active).
3. Stale `&phase=` token stripped on diagram change (no wrong-phase display).
4. Dimer face overlay: visibility uses "any polygon vertex in view" (matches the
   edge rule) so straddling square cells stay clickable.
5. `_quiverNodeClick` uses `closest("[data-qnode]")`; vertex dots get
   `pointer-events:none` so corner clicks aren't swallowed.

### Regression tests (`tests/test_inverse.py`)
- `test_quiver_seiberg_refuses_adjoint_node`
- `test_deep_dualities_stay_physical` (parametrized L131/C3Z4/rect23/trap/F0/dP3):
  exhaustive len-2/3 + fixed deep reproducers; asserts non-negative adjacency,
  positive ranks, rank-weighted anomaly freedom at every node.

Node-identity / face↔node mapping was audited (agent + empirical) and is SOUND;
the backend was already anomaly-free + involutive; dP0's rank growth is the
correct Markov triples. The negative-adjacency adjoint bug was the real defect.

## 2. BPS quivers for the 5d SCFTs — central charges + chambers (Closset–Del Zotto)

New in `conformalmanifold/bps.py`:
- **`orbifold_chamber(B)`** — the canonical finite chamber of every toric 5d
  SCFT: the n fractional branes as hypermultiplets (charges = basis e_i),
  central charges from the ε-narrow decreasing-phase recipe.  State counts
  verified: dP0→3, F0→4, dP1→4, dP2→5, dP3→6.
- **`bps_invariants(B)`** — KK/D0 charge δ=(1,…,1) (verified null vector of B),
  flavor directions = rest of ker B, flavor rank (E_n: 0,1,1,2,3 ✓).
- **`find_maximal_green_sequence(B)`** — DFS finder; correctly returns None for
  the mutation-infinite local-P² (Markov) quiver, finds one for the conifold.
- **`default_central_charges` / `orbifold_central_charges`** — default (mass,phase)
  so central charges always display in decreasing-phase order.
- `bps_quiver_from_adjacency_json`: empty sequence ⇒ canonical orbifold chamber;
  always attaches `invariants` + `canonical_chamber`; fills default Z_i.

Frontend: BPS tab defaults each geometry to the empty sequence (canonical
chamber as headline), plus new cards "BPS charge invariants" (δ, flavor, MGS)
and "Spectral network".  Nodes still clickable to explore mutation chambers.

## 3. Spectral networks (Gaiotto–Moore–Neitzke) — NEW

`bps.spectral_network_description(newton_points)` builds, from the toric/dimer
data: the SW/mirror curve = Newton polynomial; genus = Coulomb rank; sheets;
the asymptotic **(p,q) legs = outward polygon-edge normals = dimer zig-zags =
network ends** (5-brane charge conservation Σ length·leg = 0 holds); charge
lattice; wall-crossing = quiver mutation.  Wired into `api.summarize_bps_toric_quiver`
as `spectral_network`; displayed in the new BPS card.  Refs 0907.3987, 1204.4824.

## Verification
- `python -m pytest tests/ -q` → 372 passed.
- Live headless Chrome (`_qa_seiberg_bps.py <port>`): F0 four successive dualities
  (12→20→24→36 fields, path badge accumulates, 0 console errors); L131 node-2
  adjoint move refused with state preserved; BPS tab dP0 shows 3-hyper chamber,
  δ=(1,1,1), flavor 0, MGS none, spectral legs+rank.  **This caught a `const sn`
  redeclaration I introduced that had broken the whole page** — now fixed.

## 1b. Backend physics-audit fix — `quiver_seiberg` naive fallback (INCORPORATED)
A backend audit found the adjacency-only fallback (used when the tracked W
degrades to None — which is BROAD: C3, dP0, dP1, dP2, Z2xZ2, not just SPP) had
three latent defects that would give WRONG duals on non-chiral quivers (your
Y^{p,q} / L^{a,b,c} classes):
1. blanket `min(B[i][j],B[j][i])` over ALL pairs cancelled pre-existing
   vector-like bifundamentals with no mass term (e.g. spectator conifold pairs;
   gave the EMPTY C^3/(Z2xZ2) dual);
2. the same loop ran over `i==j`, negating spectator adjoints (negative mult.);
3. the `i!=j` clause dropped adjoint mesons M_ii.
Fix (`quiver_seiberg`, validated to reproduce DWZ on C^3/(Z2xZ2) node 0 and dP0):
keep adjoint mesons; only a NEW meson can be massive, only against an opposite
existing arrow; never touch the diagonal in the mass loop.  DWZ
(`wmutation.mutate`) remains authoritative when W is tracked; this is the safer
fallback.  Pinned by `test_quiver_seiberg_meson_mass_rule_matches_dwz`.
The audit also VERIFIED correct: `wmutation._reduce` (F-term coeff/sign/lone-
field/through-k-twice), `seiberg_path` arrows<->A sync, ranks through
urban_renewal.

## Open / next
- Longer (non-orbifold) finite BPS chambers for F0/dP1/dP2/dP3 are "to be found
  numerically" (Closset–Del Zotto); only the canonical orbifold chamber is shipped.
- Changes are UNCOMMITTED on `main`; review the diff (`git diff --stat`: api,
  bps, inverse, index.html, test_bps, test_inverse; +570/−31) and commit when ready.
- `_qa_seiberg_bps.py` is a throwaway CDP driver like the other `_qa_*.py`.
