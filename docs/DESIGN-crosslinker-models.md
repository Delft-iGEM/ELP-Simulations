# Design: crosslinker models beyond "two lysines touch → bond"

Status: spec written 2026-09-19 on branch `robust-audit`, before implementation.
Read `docs/LITERATURE.md` for the numbers behind the defaults.

## Why the current model is wrong about topology

`tools/crosslink.py` bonds two lysine beads that come within `distance` and
lets each lysine take part in `valence` bonds. That is right for a
bifunctional crosslinker at exactly 1:1 stoichiometry (glutaraldehyde-like,
`valence = 1`). It is wrong for anything else:

* THPP and THPC — the crosslinkers actually used on lysine-ELPs — are 3- and
  4-functional. One molecule sits at a *junction* and bridges up to f lysines
  that are all mutually connected. `valence = 2` does not reproduce that: it
  lets lysine B bond A and then C, giving a linear A–B–C–D chain of
  crosslinks, which is a different topology (it has no branch points!) and
  leads to a different modulus estimate (Miller–Macosko counts junctions of
  degree ≥ 3).
* Nothing limits the *number* of crosslinker molecules. Experiments use a fixed
  crosslinker : amine ratio (often 1:1 reactive groups); the model can consume
  every lysine.
* The crosslinker is added from solution above the coating. In a dense brush
  it reacts top-down and can be depleted before it reaches the surface. The
  implicit model crosslinks the bottom of the brush exactly as readily as the
  top.

Two models fix this in increasing order of realism and cost.

## Model A — implicit f-functional junctions (cheap, do first)

Setting `functionality` f ≥ 2 (default 2, which is exactly today's
`valence = 1` behaviour). Optional `crosslinker_ratio` = crosslinker molecules
per lysine (None = unlimited; 1/f means one reactive arm per lysine).

State: `junction_of[slot]` (−1 = free), `members[j]` list of slots, plus the
existing `reacted[bond]` flags on the all-pairs dormant bond list.

Rule at each check, nearest pair first:

1. two free lysines within `distance` → new junction {a, b} (if the junction
   budget allows); bond a–b switched on (ramped).
2. free lysine c within `distance` of any member m of a junction J with
   |J| < f → c joins J; bonds c–m for **every** m in J are switched on (they are
   all within ~distance of each other so no energy spike). Each new bond is an
   event row; a new column `junction` gives J's id so the network tool can
   rebuild the star.
3. two junctions never merge, and a lysine never leaves its junction.

Surface-pinned lysines are excluded exactly as now (`blocked`). The legacy
`valence > 1` path stays for old runs but warns that it is not a real
crosslinker functionality.

Restart: `junction_of` and `members` go into `crosslink_state.json`; the pair
list hash must match (already enforced after the min_span change).

Post-hoc replay works unchanged because `react()` is pure geometry.

## Model B — explicit crosslinker beads (the candidate "4th mode")

A crosslinker molecule becomes a single CG bead of a new residue type
(`B` in `residues_CALVADOS2.csv`: σ ≈ 0.5 nm, λ ≈ 0.3, q = 0, MW of the
reagent), added as its own single-residue CALVADOS component with `nmol =
n_crosslinker`. `build_sim` scatters them uniformly in the empty solvent
above the brush (z between the brush top and the box top minus the cutoff).

Reaction: lysine K and crosslinker bead B within `distance` → harmonic bond
K–B (r0 ≈ 0.63 nm side chain + half the spacer; see LITERATURE.md), ramped as
now. Each B has f arms; each K binds at most one B. The dormant pair list is
K × B (256 × 128 = 32k bonds for a 16-chain block ELP run; fine).

What it buys, qualitatively:

* top-down crosslinking and crosslinker depletion — the architecture
  gradient through the coating;
* stoichiometry as a real knob (`n_crosslinker`), conversion vs. ratio curves;
* junctions with the right functionality for free (a B with 3 arms *is* an
  f = 3 junction);
* the network tool sees B nodes directly.

Costs: an extra component in prepare (template change), the crosslinker beads
add to the neighbour list (+ ~10–20 % beads), the K×B pair list, and one more
reactor in `run_reactive`'s loop. Dynamics of B in CALVADOS are as unphysical
as everything else (rates are not kinetics), but *ordering* is meaningful.

Mode name: `explicit` is not an attachment mode, it is a crosslinking mode, so
it is a `crosslinker: implicit | junction | explicit` setting orthogonal to
`mode: brush | free | preattached`. It should be usable with all three.

## Also fix, independent of the model

* `min_span` (done in the audit): adjacent lysines are always 0.38 nm apart
  and must never count as a crosslink.
* `start_step` must not default to 0. Default None → 10 % of the run, at least
  the equilibration the Rg/RMSD check reports for comparable runs; refuse to
  silently start at 0 in `sim new` unless `--crosslink-start-step 0` is given
  explicitly.
* `distance` and `r0` defaults follow LITERATURE.md (CA–CA capture radius for
  a lysine-lysine bridge is ~1.2–2.0 nm, not 0.8 nm; r0 ≈ 1.0–1.5 nm; k soft,
  a few hundred kJ/mol/nm²).
