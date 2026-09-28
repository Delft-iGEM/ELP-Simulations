# Session summary — September 2026

Covers a long working session (Fable 5.1, continued on Opus 5, wrapped up on
Sonnet 5) that started as a bug audit of the CALVADOS-based ELP brush pipeline
and grew into the RGD placement study. Written so the session can be cleared
without losing the conclusions. Cross-references `docs/BUGS.md`,
`docs/GAPS.md`, `docs/LITERATURE.md`, `docs/RERUN.md`, `docs/V3-BATCH.md`,
`docs/DESIGN-crosslinker-models.md`, `docs/UNFINISHED.md` for detail; this file
is the index of what was concluded and why it matters.

## 1. The audit and the batch-destroying bug

The pipeline runs CALVADOS 2 (git commit `dc5beac158b72262f585947487c1dd8ef1602a1`)
on OpenMM 8.2.0, one bead per residue, Ashbaugh-Hatch + Debye-Hückel, Langevin
dynamics at dt = 10 fs, friction 0.01 ps⁻¹. Surface-grafted brushes fix residue 0
(tagged "Z", zero mass) to a wall; a reactive crosslinker turns dormant
zero-strength `HarmonicBondForce` terms on between nearby lysines.

**BUG-XL-7 (critical, fixed):** `Crosslinker.add_force` in `tools/crosslink.py`
never called `force.setUsesPeriodicBoundaryConditions(True)`. The reaction
criterion itself used minimum-image distances, so two lysines 0.6 nm apart
across a periodic boundary were correctly identified as a reactive pair — but
the bond OpenMM then added used the *raw* (non-minimum-image) distance. Every
crosslink formed across a boundary was stretched across the whole box, wrecking
the run. This destroyed **all 12 of the crosslinking runs** in the first
block-runs batch (SLURM jobs 265624–265643); the 8 non-crosslinking runs in
that batch were unaffected. Fixed, and a regression test was added that fails
if the fix is reverted. Full list of the 17 bugs found is in `docs/BUGS.md`;
other notable ones:

- Wall spring constant was implicitly k=1 (way too soft) — set to `WALL_K =
  5000.0`, plus `Z_ANCHOR = 2.0` so the anchor sits above the wall rather than
  inside it. Sub-wall occupancy went from 26.3% to 0.23%.
- `new_simulation.py` was only recentring x/y on `xinit`, not z.
- `surface.py` used a 10 nm build height, limiting footprint to 5.32 nm and
  capping achievable surface density; raised to 20 nm (footprint 3.80 nm),
  which is what made the 0.04/0.05 concentration sweeps possible.
- Various off-by-one / frame-indexing bugs in `crosslink.py`'s state file and
  `post_hoc_events`.
- A race in `template/job.sh` (`uv run python` vs an already-built venv) killed
  3 of an early 125-job batch on missing matplotlib style files; switched to
  invoking `.venv/bin/python` directly.

**"K-K crosslinking in the first few picoseconds"** — the user's original
worry — was addressed by the `min_span` parameter (default 3) in the reaction
criterion, which excludes lysine pairs too close in sequence to be a real
crosslink, and is independent of the PBC bug above.

All of this was done on the `robust-audit` branch (separately, a worktree also
exists at `~/ELP-Simulations-robust` — see project memory) and merged into
`main` afterward, preserving main's block-run outputs and the CSV analysis
functionality that had landed there in parallel.

## 2. Model validation against literature

Free-solution ELP4 runs (no source changes needed — CALVADOS's own
`topol: 'center'`, `ext_force: false` defaults already give free solution) were
compared against published hydrodynamic radii. `tools/hydrodynamic.py`
implements Kirkwood-Riseman Rh (pre-averaged ⟨1/r⟩, not per-configuration) with
proper chain unwrapping across periodic images before computing distances.

**Result: CALVADOS over-compacts ELPs by 4–17%** relative to measured Rh
(commit `133593f`). This is a property of the force field, not a bug in this
repo — noted in `docs/LITERATURE.md` (651 lines of comparison) as a known
systematic and something to caveat in any quantitative claim about absolute
chain size.

## 3. Lysine arrangement — the dominant design lever

20 distinct lysine arrangements were run on one fixed scaffold (same VPGIG /
RGD / VPGMG spacing, His-tag terminal, only the 4-lysine block's position
varied), 400 ns each, preattached mode at fraction 0.3, concentration 0.04.

**Result: lysine arrangement changes the inter-chain contact share by ~136×**
(0.6% → 87% of all K-K contacts being inter-chain rather than intra-chain),
depending purely on where the lysine block sits in the sequence — commit
`4d5e0dd`. This is the single largest design effect found all session, and
the practical conclusion is that **crosslinking/bridging behaviour is
controlled almost entirely by lysine placement**, not by concentration or
grafting fraction.

Contact analysis conventions settled during this work, applied to all
subsequent contact CSVs:
- 12 Å cutoff (not 10 Å — a locale bug briefly made an earlier CSV misleadingly
  show "100" instead of "10" because Dutch-locale Excel reads `.` as a
  thousands separator; fixed by writing CSVs with `;` delimiter, `,` decimal,
  UTF-8 BOM, in `tools/analysis_workbook.py`'s `write_contacts_csv`).
- Tethered (surface-pinned) K-K contacts are excluded from crosslinking
  contact counts — including them was found to inflate intra-chain contacts
  by 57–93% as a frozen-pin artifact (visible in the 200 ns vs 600 ns
  comparison of the concentration sweep; the "more pinning lowers inter share"
  conclusion drawn at 200 ns was contradicted once the 600 ns data came in and
  was retracted).
- Equilibration is checked in production analyses via a Chodera-style burn-in
  test (statistical inefficiency g, effective sample count N_eff) before
  contact statistics are trusted.

## 4. RGD domain placement — the main topic of the back half of the session

Question: does moving the 4 copies of the 23-residue motif
`LDASTVYAVTGRGDSPASSAASA` (which contains the RGD tripeptide) around the
sequence change how exposed it is at the brush surface?

**Method:** `tools/sasa.py` implements Shrake-Rupley SASA directly with
CALVADOS's own σ/2 bead radii (mdtraj's built-in version uses element radii —
carbon 0.17 nm — which is meaningless for coarse-grained beads that represent
whole residues at 0.45–0.68 nm). Probe radius 0.14 nm (water convention).
Periodic images are respected so a bead can be correctly occluded by a
neighbouring chain wrapped to the far side of the box.

24 distinct arrangements of the motif were run (preattached mode, fraction
0.3, concentration 0.04), and analysed for: motif SASA, RGD-only SASA, lysine
contacts, and average z-height of the RGD blocks in the layer (commit
`1e3f5a2`, results in `d83f382`).

**Headline result: placement barely moves exposure.** Across all 24 designs,
exposure ranges roughly 39.6%–42.1% of the free-domain reference — a ~2.5-point
spread — and it *trades off against bridging* (designs that expose the domain
most tend to reduce inter-chain lysine bridging, and vice versa).

**The one placement effect that is real:** the C-terminal slot. A paired
within-run test — comparing a C-terminal copy against its own siblings inside
the *same* run, which cancels seed effects exactly — found the C-terminal copy
significantly more exposed in **8 runs out of 8** (p = 0.004). This is the only
causally-established placement effect from this study.

**Best vs worst design (this session's final analysis, commit `f5786b2`,
`tools/sasa_two_sample.py`):** comparing `rgdarr08-all-cterm` (best, all 4
copies near C-terminal-favourable positions) against `rgdarr23-staircase`
(worst):

| | best (`all-cterm`) | worst (`staircase`) | difference |
|---|---|---|---|
| whole 23-mer, 4 copies | 79.20 ± 0.93 nm² (n=16 chains) | 77.46 ± 0.73 nm² | +1.74 nm² = +0.89 pts, Welch p = 2.5e-06 |
| RGD only, 4 copies | 11.73 ± 0.15 nm² | 11.57 ± 0.12 nm² | +0.17 nm² = +0.63 pts, p = 0.0013 |

Statistically significant, but:
- Only **2.2% relative** effect on the 23-mer, 1.5% on RGD alone.
- Individual chains overlap heavily (77.5–80.6 vs 75.8–78.5 nm² range) — the
  small p-value comes from averaging 16 chains, not from real separation
  between the distributions.
- The two runs are different seeds as well as different designs, so this test
  says "these two runs differ," not "the design causes it" — design and seed
  are confounded at the between-run level. Only the paired within-run test
  (§ above) is causal. About half the best-vs-worst gap is attributable to
  that one C-terminal effect (1.7 pts on one copy ≈ 0.42 pts spread over four,
  against 0.89 observed).

**Grafting cost:** comparing the same 24 sequences run as single free-floating
chains (no surface) vs grafted in the brush, grafting costs only **~4%** of
domain exposure (`docs/` SASA workbook, `sasa-grafted-vs-free.xlsx` from
commit `65047f9`). Same qualitative trend (best/worst ranking of designs)
holds free vs grafted, so the placement effect is intrinsic to the sequence,
not an artifact of the surface.

### Bottom line for design choices

**Lysine arrangement (136× effect on bridging) dominates over RGD/domain
placement (2.2% effect on exposure) by roughly two orders of magnitude.**
Domain placement is close to a free design parameter — put the C-terminal
copy at the C-terminus because it's free to do, but don't expect it to matter
much — and design effort should go into lysine placement if bridging /
crosslinking behaviour is the target property.

## 5. Statistics housekeeping (why the "noise" looked large)

Three different spreads were being conflated when the 24-design sweep first
came back:
- **Frame-to-frame scatter within one chain's trajectory** (~0.77 nm²) — not
  an error bar, just thermal fluctuation of the domain's exposure over time.
- **Chain-to-chain spread within one run** (SEM ≈ 0.06 nm² once averaged over
  16 chains and the production window) — this is the real precision on a
  single design's mean.
- **Design-to-design spread across the 24 runs** (SD ≈ 0.25–0.52 nm²) — this
  is what actually separates or fails to separate designs.

Resolving the "why is the noise so large" question meant keeping the
per-chain replication rather than collapsing it, which led to writing the
paired within-run test (causal, cancels seed) and the best-vs-worst Welch test
(§4, correlational, confounded with seed) as two separate tools with two
separate claims each is licensed to make.

## 6. Tooling added this session (all in `tools/`, all tested)

- `tools/sasa.py` — CG-correct Shrake-Rupley SASA, `find_motif`, `chain_sequence`.
- `tools/hydrodynamic.py` — Kirkwood-Riseman Rh with chain unwrapping, Rg,
  burn-in-aware SEM via statistical inefficiency.
- `tools/analysis_workbook.py` — per-run collection, `lysine_architecture`
  parsing (had an aliasing bug — `runs.append(cur); cur.clear()` was emptying
  the just-appended list because it stored a reference, not a copy; fixed by
  appending `list(cur)`), locale-correct CSV writer.
- `tools/build_sasa_workbook.py` — 3-sheet grafted-vs-free Excel workbook with
  a Statistics sheet documenting the three spreads above.
- `tools/sasa_two_sample.py` — per-chain paired/unpaired Welch comparison
  between two designs' total domain SASA, with the seed-confound caveat
  written into the module docstring so it travels with the code.
- `tools/analysis_batch.py` gained conditional step guards
  (`_Need(step, unless=...)`) so an analysis step can be skipped when its
  precondition wasn't run — added after the user caught an unnecessary
  equilibration re-run being triggered by a tag that should have been
  conditional.

All 63 tests pass as of the last commit (`f5786b2`) on `main`.

## 7. Known unfinished / suggested-but-not-run follow-ups

- **Seed vs design confound is not resolved.** The natural next step (offered,
  not requested) is 3 seeds × 4 contrasting designs (~2 GPU-hours) to properly
  separate design effect from seed effect in the between-run comparisons.
- **Crosslink modulus (Miller-Macosko / phantom network / RENT) cannot be
  computed yet.** The current crosslinker is single-bond, not the
  f-functional junction model that a modulus calculation needs. Design is
  written up in `docs/DESIGN-crosslinker-models.md` but not implemented.
- **`sim network` CLI command is missing its entry point** — the module and
  14 tests for it exist per `docs/UNFINISHED.md`, but nothing wires it into
  the CLI.
- Other bugs/gaps not yet actioned are catalogued in `docs/BUGS.md` and
  `docs/GAPS.md` — worth a re-read before starting new physics work.

## 8. Where things stand on git

- Branch: `main`, everything above is merged and pushed (through commit
  `f5786b2`).
- Batch 265624–265643 (block-runs, submitted 2026-09-19): the 12 crosslinking
  runs in that batch are corrupted by BUG-XL-7 and should not be trusted for
  crosslinking conclusions (see project memory
  `project-batch-265624-corrupted.md`); the 8 non-crosslinking runs are fine.
  The 12 were re-run post-fix (`docs/RERUN.md`).
- v3 batch (24 runs: 4 families × 3 modes × crosslink on/off, concentration
  0.015, crosslink prob 0.1, 20M steps) was created and submitted; see
  `docs/V3-BATCH.md`.
