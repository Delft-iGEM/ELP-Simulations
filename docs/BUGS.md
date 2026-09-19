# BUGS

Branch `robust-audit`. One entry per confirmed bug. Every entry says how it was
*measured*, not just what looked wrong — if you are an AI picking this up later,
re-measure before you trust it.

Status key: **fixed** (changed on this branch, test covers it) · **open** ·
**wontfix (documented)**.

The audit was cut short by a credit limit partway through; `docs/UNFINISHED.md`
lists what was still in flight, so do not read a quiet area as a clean one.

---

## BUG-XL-7 — crosslink bonds ignored the periodic boundary, and destroyed every reactive run

**File:** `tools/crosslink.py` `Crosslinker.add_force`.
**Severity:** CRITICAL. This silently corrupted **all 12 crosslinking runs** of
the 2026-09-19 block batch (SLURM 265624–265643). It is one missing line.

`add_force` builds a plain `HarmonicBondForce()` and never calls
`setUsesPeriodicBoundaryConditions(True)`. Every force CALVADOS builds itself
does (`calvados/interactions.py`: bonds, angles, restraints, scaled LJ/YU,
FENE), and the surface tethers get periodicity via `periodicdistance()` inside
their expression. The crosslink force was the only one without it.

The reaction criterion, correctly, uses the **minimum image**
(`min_image_distances`). So two lysines 0.6 nm apart *across* a box edge are a
legitimate pair and do get bonded. OpenMM then evaluates that bond on the **raw
coordinate difference** — one box length, 28.28 nm in the block runs — against
r0 = 0.6 nm. With k = 2000 kJ/mol/nm² the bond pulls with

    F = k·Δr = 2000 × 27.7 ≈ 5.5e4 kJ/mol/nm,  E ≈ 7.7e5 kJ/mol

on two beads of ~128 Da. The system is destroyed within a few thousand steps.

**Evidence.** `monoblock-brush-xl` (720 residues, 16 chains, 28.28 nm box,
`crosslink_start_step` = 1 000 000):

| frame | MD step | max abs coordinate |
| --- | --- | --- |
| 142 | 1 011 010 | 46 nm — healthy |
| 143 | 1 018 080 | **1.1e11 nm** |
| 999 | 7 070 000 | 1.6e13 nm |

The third bond ever formed (step 1 005 500, chains 2→15) has a raw separation of
24.1 nm against a minimum-image separation of 9.1 nm — i.e. it spans the
boundary. Divergence follows within two checks. The same signature appears in
every crosslinking run checked (`tetrablock-brush-xl` frame 143,
`triblock-64-pre-xl` frame 150), while all eight non-crosslinking runs
completed cleanly.

**Why nobody noticed.** Nothing reports it. OpenMM does not error on huge
coordinates, the `StateDataReporter` is configured without
`potentialEnergy`, and `metadata.csv` records "completed". The five runs that
*did* fail crashed only at the very last line, writing the final PDB
(`ValueError: coordinate "1670550345991.72" could not be represented in a
width-8 field`); the six that hit the 2 h walltime were never flagged at all
and would have been read as valid results. The earlier reactive test runs
(`does-link-work2`, `xl-dense-reactive`) and every audit smoke run produced
**only intra-chain events in small boxes**, so no bond ever spanned a boundary
and the bug stayed invisible — 114 of the 123 events in `monoblock-brush-xl` are
inter-chain.

**Fix:** `force.setUsesPeriodicBoundaryConditions(True)`, with a comment
explaining why. Regression test
`tests/test_crosslink.py::test_crosslink_force_uses_periodic_boundaries` builds
a two-bead system 0.5 nm apart through a 28 nm boundary and asserts the bond
energy is 10 kJ/mol (minimum image) and not 7.2e5 kJ/mol (raw). Verified to
fail with the line removed and pass with it restored. **Status: fixed.**

**Consequence for existing data:** every crosslinked trajectory produced before
this fix is physically meaningless from the first boundary-spanning bond
onwards, and so is every crosslink event recorded after that point (on
coordinates of 1e13 nm, minimum-image distances are noise). The non-crosslinked
runs are unaffected by this particular bug. Re-run the reactive batch.

**Related:** the same runs were 6–9.5x slower than their non-reactive
counterparts (1350 and 841 it/s against 7988 it/s; 1.02e3 vs 7.43e3 ns/day),
which is why six of them hit the 2 h walltime. The reaction loop calls
`context.getState(getPositions=True)` every `check_every` = 1000 steps, forcing
a full GPU→CPU position sync, and `updateParametersInContext` re-uploads all
32 640 bond parameters on every change. A 70 ns reactive run of these systems
needs ~2.5–3 h, not 2 h. See `docs/GAPS.md`.

## BUG-WALL-1 — the surface wall is ~5000x too soft (it barely exists)

**File:** `tools/new_simulation.py` (`PREPARE_TEMPLATE`, `ext_force_expr`), every
generated `simulations/*/prepare.py`, every `runtime/config.yaml`.
**Severity:** high — affects every run ever done, including the 20 queued jobs.

The surface is CALVADOS's one-sided external potential:

```
ext_force_expr = f'step({z_wall}-z)*0.5*({z_wall}-z)^2'
```

There is no force constant in the expression, so OpenMM reads k = 1 kJ/mol/nm².
At 293 K, kT = 2.44 kJ/mol:

| penetration below z_wall | energy cost | in kT |
| --- | --- | --- |
| 0.1 nm | 0.005 kJ/mol | 0.002 |
| 1.0 nm | 0.5 kJ/mol | 0.2 |
| 1.9 nm (reaching z = 0) | 1.8 kJ/mol | 0.74 |

A bead crosses the entire "surface" for less than one kT. The box is periodic in
z and the wall is not a lid, so a bead that keeps going reappears at the top.

Measured on the finished brush runs (`tetrablock`, `monoblock-test`,
`triblock-64-test`): 1.4–1.8 % of all bead samples sat below the wall onset,
0.5 % below z = 1.0 nm, deepest bead z = −0.66 nm. For the lysine-containing
runs, 6.6 % of all *lysine* samples were below the onset — which no
surface-binding or near-surface crosslink statistic survives.

**Fix:** `WALL_K = 5000.0` kJ/mol/nm², written into the expression
(`step(z_wall-z)*0.5*wall_k*(z_wall-z)^2`). 0.1 nm of penetration now costs
10 kT and the thermal penetration depth sqrt(kT/k) is 0.02 nm. It stays below
CALVADOS's own bond stiffness (kb = 8033 kJ/mol/nm², same bead masses), so it
adds no timescale the 10 fs Langevin step does not already resolve, and the form
stays harmonic so `tools.metadata._z_wall_from_expr` still parses it.
**Status: fixed.** Independently flagged by the literature review as LIT-2
(which recommends ≥ 2500).

**Verified before/after** on a 4-chain, 40-residue, 6000-step CPU run
(`zz-audit-wall-before` / `zz-audit-wall-after`, identical except for the two
wall fixes). Both fixes are in play together, which is the point — neither is
safe alone:

| | anchor z at t=0 | beads below wall at t=0 | bead samples below wall over the run | lowest bead ever |
| --- | --- | --- | --- | --- |
| before | 1.430 nm | 46.7 % | 26.3 % | 1.107 nm |
| after | 2.000 nm | 0 % | 0.23 % | 1.851 nm |

A 114x reduction in sub-wall occupancy, the anchor exactly where it is meant to
be, the deepest excursion now 0.05 nm (thermal, sqrt(kT/k) = 0.02 nm), and no
NaN or instability from the stiffer wall.

## BUG-WALL-2 — the brush anchor starts *below* the wall, and half the chain with it

**File:** `tools/new_simulation.py` `PREPARE_TEMPLATE.build_sim`.
**Severity:** high — affects every brush run, including the 20 queued jobs.

CALVADOS's `build_compact` lays a chain out as a serpentine through a cube of
side ≈ cbrt(N) bonds **centred on the origin**, so bead 0 is the cube's lowest
corner, at −0.5·cbrt(N)·0.38 nm in x, y *and* z. The template recentred only
x and y, then added a flat 2.0 nm in z:

```python
xy_centre = 0.5 * (comp.xinit[:, :2].min(axis=0) + comp.xinit[:, :2].max(axis=0))
xinit_centred = comp.xinit - np.array([xy_centre[0], xy_centre[1], 0.0])
pos = xinit_centred + np.array([x0, y0, 2.0])
```

Measured on the prepared `runtime/top.pdb` of the **currently queued** batch:

| run | residues/chain | anchor z | min bead z | beads below the 1.9 nm wall |
| --- | --- | --- | --- | --- |
| `tetrablock-brush-xl` | 1040 | 0.100 nm | 0.100 nm | 45.7 % |
| `triblock-64-brush-xl` | 1040 | 0.100 nm | 0.100 nm | 45.7 % |
| `monoblock-brush-xl` | 720 | 0.480 nm | 0.480 nm | 44.4 % |

The anchor is a zero-mass bead held **fixed** by OpenMM, so it sits 1.4–1.8 nm
*inside* the surface for the whole run and cannot relax out. Roughly 45 % of
every chain starts under the wall. With BUG-WALL-1's k = 1 wall this merely
leaked; with a stiff wall it would launch the system, so the two must be fixed
together.

The x/y recentring was wrong in its own right too: the bounding-box centre of
the cube is not bead 0, so the anchor sat ~2 nm off its lattice point laterally
and the grafting pattern did not match the lattice the box was built around.

**Fix:** shift by bead 0's own position, `xinit_anchored = comp.xinit - comp.xinit[0]`,
then place at `[x0, y0, Z_ANCHOR]` with `Z_ANCHOR = 2.0` nm. The anchor lands
exactly on its lattice point just above the wall onset and, because bead 0 is
the minimum corner, every other bead starts above it. `build_sim` now raises if
any bead is below `z_wall`, which also pins down the `build_compact` assumption
if CALVADOS ever changes it. **Status: fixed.**

## BUG-BOX-1 — `box_height` from the notebook's `Shared` was never passed on

**File:** `analyze.ipynb` (batch submit cell and local smoke-test cell).
**Severity:** medium.

`plan_run` showed a per-run box height in the planning table, but the call to
`create_simulation` omitted `box_height=`, so the run was built with the default
rule instead. The table and the simulation disagreed whenever `Shared` set a
height. **Fix:** pass `box_height=shared.box_height` in both cells;
`plan_box_height` now also warns when an explicit height is below the default.
**Status: fixed.**

## BUG-XL-1 — adjacent lysines "crosslink" through the backbone

**File:** `tools/crosslink.py` (pair construction).
**Severity:** high for any K-rich sequence — directly corrupts the primary output.

Every lysine pair was a candidate, including neighbours along the same chain.
Bonded CA beads are 0.38 nm apart and second neighbours 0.5–1.2 nm, i.e. always
inside any sensible reaction distance (0.8 nm, let alone the 1.2 nm the
literature review recommends). A `...KK...` motif therefore "reacted" on the
first check, at step 0, forever — and the resulting bond constrains nothing,
because it is parallel to a covalent backbone bond that is already there. These
fake events inflate the intra-chain count, which is *the* headline number (the
primary-loop estimate feeding elastic network theory).

**Fix:** new setting `crosslink_min_span` (default 3): intra-chain pairs closer
than this many residues along the backbone are dropped at construction, so they
are not even dormant bonds. Inter-chain pairs are never affected; `min_span=0`
restores the old behaviour. The summary reports how many pairs were excluded.
**Status: fixed.**

## BUG-XL-2 — a restart could reapply bond indices to the wrong pairs

**File:** `tools/crosslink.py` `save_state` / `load_state`.
**Severity:** high (silent wrong network), latent before BUG-XL-1.

The checkpointed network is stored as *bond indices* into the all-pairs list.
That list is only meaningful together with the settings that built it — and
`min_span` changes it. Restarting with a different `min_span`, or a different
site list, would restore `reacted` flags onto entirely different pairs.

**Fix:** `save_state` records the pair count and a hash of the pair list;
`load_state` refuses a mismatch, naming both values, instead of continuing into
a silently different network. A state file written before the hash existed is
accepted as "every pair kept". Also added: a consistency check that `used[]`
agrees with the bond counts implied by `reacted[]`, so a corrupt or hand-edited
state file fails loudly rather than mis-counting valence. **Status: fixed.**

## BUG-XL-3 — the event `frame` column was off by one

**Files:** `tools/crosslink.py` and `tools/surface.py` (`event_rows`),
`analyze.ipynb` (time axis).
**Severity:** medium — every event was attributed to the wrong frame.

The code used `frame = round(step / n_save)`. OpenMM's `DCDReporter` writes its
**first frame after `n_save` steps, not at step 0**, so saved frame k holds MD
step (k+1)·n_save and the frame nearest step s is `round(s / n_save) - 1`. Every
event pointed one frame late; a movie flashing a bond, or any per-frame join
against the trajectory, was consistently out. A bond present at initialisation
(step 0) predates every frame and is put on frame 0.

Same off-by-one in the notebook's time axis: `np.arange(n_frames) * unit_per_frame`
labels the first frame t = 0, which does not exist → `(np.arange(n_frames) + 1)`.

Verified against real smoke output: events CSVs written before the fix satisfy
the old convention, those written after satisfy the new one. **Status: fixed.**

## BUG-XL-4 — post-hoc replay silently ignored unmatched surface bonds

**File:** `tools/crosslink.py` `post_hoc_events`.
**Severity:** medium — breaks the "never both" invariant exactly when it matters.

Surface-bonded lysines are kept out of the post-hoc pair pool by matching each
`kind=surface` row to a site via `(chain, resid)`. An unmatched row was skipped
by `if slot is not None:` — silently. The replay would then crosslink a lysine
the live run had bonded to the surface, producing precisely the "bonded to the
ground *and* to another lysine" state the masking exists to prevent, with no
warning. Triggered by any events-file/trajectory mismatch: regenerated topology,
hand-edited CSV, a run analysed against the wrong simulation.

**Fix:** raise, naming the offending key and the keys the topology actually has.
**Status: fixed** (found while triaging a failing test; now covered by one).

## BUG-XL-5 — replay compared `start_step` against the wrong step

**File:** `tools/crosslink.py` `post_hoc_events`.
**Severity:** medium.

Same root cause as BUG-XL-3, on the other side: the loop used
`step = frame * n_save`, so frame 0 was treated as step 0. With
`start_step = n_save` the first real frame was skipped entirely, and every
surface-bond time comparison sat one frame out of step with the live run's.
**Fix:** `step = (frame + 1) * n_save`. **Status: fixed.**

## BUG-XL-6 — conversion counted lysines that cannot react

**File:** `tools/crosslink.py` `summary_lines`.
**Severity:** low (reporting only).

`conversion = 2 x events / (n_sites x valence)` counts surface-pinned lysines in
the denominator, but a pinned lysine can never crosslink. In `preattached` mode
with 40 % of lysines pinned this understates conversion by that fraction.
**Fix:** report capacity excluding pinned sites alongside the raw number.
**Status: fixed.**

## BUG-RUN-1 — a crashed or killed run left stale metadata

**File:** `template/run.py`.
**Severity:** medium (you lose a run's provenance).

`write_timing` / `write_metadata` ran only after `simulate()` returned. A NaN, a
SLURM walltime SIGTERM, or a bug in the reaction loop left `metadata.csv` still
saying what `prepare.py` planned — "not started", or the planned frame count —
for a run that had actually produced most of a trajectory.
**Fix:** `try: ... finally: _record_run()`, so status, frames written and
elapsed time are recorded on the way out and the exception is re-raised
untouched. **Status: fixed.**

## BUG-RUN-2 — the CPU fallback never triggered, then ran single-threaded

**File:** `template/run.py`.
**Severity:** medium.

The guard compared the requested platform against the *names* OpenMM lists. An
OpenMM built with CUDA support lists `"CUDA"` on every machine, GPU or not — it
only fails when a `Context` is created. So a CUDA run on a CPU-only node did not
fall back; it crashed. **Fix:** probe by actually creating a one-particle
Context on that platform. Additionally CALVADOS passes `config['threads']` to
the CPU platform (default 1), so a fallback inside a SLURM job would have left
every other allocated core idle — the fallback now takes threads from
`SLURM_CPUS_PER_TASK`. **Status: fixed** (smoke run `zz-audit-smoke-cuda-fallback`).

## BUG-BATCH-1 — SLURM walltimes were parsed with the wrong grammar

**Files:** `tools/run_batch.py` `parse_walltime`, `tools/metadata.py` `parse_slurm_time`.
**Severity:** medium — silently wrong job limits and wrong cost totals.

Both parsers read a two-field `H:MM` as hours:minutes. To SLURM, two
colon-separated fields are **minutes:seconds** (`sbatch(1)`), so a CSV cell
saying `2:00` buys two *minutes*, not two hours, and the notebook's walltime
totals and `walltime_used_fraction` were out by 60x for such rows. The forms
`D-HH` and `D-HH:MM` were not accepted at all.
**Fix:** both now implement SLURM's full grammar with SLURM's meaning.
`block-runs.csv` uses `2:00:00` throughout, so the queued batch is unaffected.
**Status: fixed.**

## BUG-BATCH-2 — the ELP shorthand's `G1 A1` is not a GAGA block

**Files:** `tools/run_batch.py` (reporting), `tools/elibpy.py` (behaviour).
**Severity:** high as a *scientific* issue; the code does what it documents.

In the shorthand, a token `X<n>` expands to n copies of the pentapeptide
`VPGXG`. So `G1 A1 G1 A1 ...` is **`VPGGGVPGAGVPGGGVPGAG...`**, not the
alanine-rich `GAGAGAGA` elastic block it reads as. Only a parenthesised
fragment, `(GAGAGAGA)`, is inserted literally — and `Example-csv.csv` uses
exactly that form, while `block-runs.csv` uses `G1 A1`.

Verified by expansion: the "triblock" mid-block in `block-runs.csv` comes out as
`...VPGGGVPGAGVPGGGVPGAG...`, i.e. more VPGxG repeats, not a GAGA block. The 20
queued jobs are therefore running **pentapeptide-only sequences** — a legitimate
ELP construct, but not the block copolymer the run names imply. The residue
counts (720/880/1040) are consistent either way, so nothing looked wrong.

**Fix (reporting only — the expansion itself is documented behaviour and is
unchanged):** `print_summary` now states that one letter = one VPGxG
pentapeptide, says only `(...)` fragments are literal, and prints the first 60
residues of the real sequence. A bare letter with no repeat count (`K V3`) is
now refused outright, since elibpy silently inserts it as a single literal
residue — almost never what a shorthand cell means and invisible in the residue
count. **Status: fixed (reporting).** Whether the queued batch is the chemistry
you meant is a scientific decision, and it is yours.

## BUG-BATCH-3 — a leading lysine is silently deleted in brush mode

**File:** `tools/new_simulation.py` `create_simulation`.
**Severity:** medium.

Brush mode rewrites residue 0 as the `Z` anchor tag, which CALVADOS treats as a
valine. If the sequence starts with `K`, that lysine is simply not in the
simulated chain: the run has one fewer reactive site than the sequence says, and
nothing said so. Several `block-runs.csv` rows begin with `K1`, so those chains
have 15 reactive lysines, not 16.
**Fix:** `sim new` warns, naming both counts, and suggests prepending a spacer.
**Status: fixed (warning).** Deliberately not auto-corrected — changing the
sequence behind your back would be worse.

## BUG-AN-1 — notebook start-frames were hard-coded past the equilibration check

**File:** `analyze.ipynb` (z-distribution cell, contact-counting cell).
**Severity:** medium.

`start_frame = 600` and `contact_start_frame = 0` were literals sitting below a
cell that computes the burn-in properly. Analyses silently used a burn-in
unrelated to the measured one — 600 frames for one run, none for another.
**Fix:** both default to `production_start_frame` from the equilibration cell,
overridable by typing a number. **Status: fixed.**

## BUG-AN-2 — a comment inviting someone to delete the `box=` argument

**File:** `analyze.ipynb` (contact-counting cell).
**Severity:** low (comment), but it guards a real invariant.

The comment claimed chains are "stored unwrapped". They are not: the DCD stores
whole molecules wrapped into the box. The `box=` argument is exactly what makes
a contact across the periodic edge count. **Fix:** comment corrected.
**Status: fixed.**

## BUG-AN-3 — `_Tee` broke any library that probes stdout

**File:** `tools/analysis_io.py`.
**Severity:** low, but it aborts whole notebook cells.

`save_text` replaces `sys.stdout` with a `_Tee` exposing only `write`/`flush`.
typer/click, tqdm and rich probe `isatty`, `encoding`, `fileno`; an
`AttributeError` from inside a `with save_text(...)` block killed the cell.
**Fix:** `__getattr__` delegates to the real stdout. **Status: fixed.**

## BUG-META-1 — `z_wall_nm` was documented as a well depth

**File:** `tools/metadata.py`.
**Severity:** low (documentation), but it is how BUG-WALL-1 stayed hidden.

The field was described as "Depth of the surface-attachment potential well". It
is not a depth and it is not a well: it is the *onset height* of a purely
repulsive wall, and its stiffness — the actual depth-like quantity — was not
recorded at all. **Fix:** description corrected; `ext_force_expr` was already
recorded and now carries `WALL_K`. **Status: fixed.**

---

## Not bugs, but worth knowing (checked, no change needed)

* **`System.addForce()` takes ownership of the force.** A temporary
  `openmm.System()` is garbage-collected immediately and deletes the force with
  it, so every later `getParticleParameters` call hits a dangling pointer and
  raises `std::bad_array_new_length`. Production always passes the long-lived
  `sim.system`, so this only bites tests — both test files now keep the System
  alive and say why.
* **The ramp start step is recorded only when a force object exists.** Post-hoc
  replay has no force and nothing to ramp, so this is correct; but a
  `Crosslinker`/`SurfaceBinder` used as a stand-in for a live run *without*
  calling `add_force` will silently report every bond as fully ramped.
* **The all-pairs dormant bond list is fine at this scale.** 16 chains x 16
  lysines = 256 sites = 32 640 bonds; `updateParametersInContext` re-uploads all
  of them, a few hundred kB, a handful of times per reaction — negligible beside
  a 1000-step chunk of MD. The 500-site warning is still the right threshold.
* **x/y wrapping is bookkeeping, not physics**, and minimum-image contact
  counting is correct across the boundary (already in `notes_titus`; re-verified).

## Verified invariants (crosslink / surface)

Covered by `tests/test_crosslink.py` and `tests/test_surface.py` — 47 tests,
7 s, no GPU:

* a surface-pinned lysine never forms a K–K bond, and a K–K bonded lysine never
  pins to the surface — in the same check step, at valence ≥ 2 where `used > 0`
  but the site is not saturated, across a save/load round trip, and in post-hoc
  replay;
* no reaction before `start_step`, for either reactor, including across chunk
  boundaries and when the two reactors have different `check_every`;
* bond stiffness rises monotonically 0 → k over `ramp_steps`, never jumps, and
  resumes mid-ramp at the fraction it had reached;
* restart determinism: same seed and checkpoint step → same subsequent draws;
* `used[]` stays consistent with `reacted[]`;
* the reaction rule is pure geometry, so the live run and the replay apply an
  identical rule.
