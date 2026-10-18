# Opt-in models

Reactive crosslinking, CALVADOS-TEA, ionic strength and lysine attachment. All are off by default.

## Reactive crosslinking (opt-in)

Off by default, and the default is not merely "the feature does nothing": with
no `--crosslink-distance` the run never writes `runtime/crosslink.yaml`, `run.py`
never imports the crosslinking module, no force is added and no random numbers
are drawn. A non-reactive run is bit-for-bit the trajectory it was before this
existed — verified by running the same seeded system through both the old and
new pipelines and comparing coordinates.

Turned on, lysine pairs that come within the reaction distance bond permanently
*while the run is going*, so a formed crosslink tethers its chains and changes
what happens next:

```bash
sim new ... --crosslink-distance 0.8 --crosslink-start-step 2000000
sim new ... --crosslink-distance 0.8 --crosslink-valence 2   # trifunctional (THPP-like)
sim new ... --crosslink-distance 0.8 --crosslink-prob 0.01   # reaction- not diffusion-limited
```

`--crosslink-distance` is in **nanometres**, like the rest of the box maths and
OpenMM — not the ångströms an MDAnalysis-based post-hoc script works in. Every
startup line prints both so the two cannot be confused.

The knobs: `--crosslink-valence` (bonds per lysine; 1 = a bifunctional
crosslinker at 1:1, 2 = trifunctional junctions), `--crosslink-prob` (P(react |
within the cutoff at a check) — below 1 decouples the reaction rate from the
diffusion rate, which is the Damköhler sweep), `--crosslink-check-every`,
`--crosslink-k` / `--crosslink-r0` / `--crosslink-ramp-steps` (the bond itself,
switched on gradually so it doesn't spike the temperature), and
`--crosslink-start-step`, which you should always set — see below. The same
names are inputs in the notebook's planning cell.

Each reactive run writes, beside its trajectory: `crosslink_events.csv` (one row
per bond, columns `frame,time_ps,resid_i,chain_i,resid_j,chain_j,kind,span,x,y,z,
pymol_i,pymol_j,origin` — `frame` is the saved-trajectory frame the event rounds to,
`chain_*` are 0-based chain indices, `resid_*` are 1-based within their chain,
`x,y,z` is the midpoint of the pair, `pymol_*` are 1-based for PyMOL's
`index` selector, and `origin` is `run` for a crosslink — `initial` only ever
appears on the lysine–surface rows of the lysine-attached modes) and
`crosslink_summary.txt` (conversion, the
intra/inter split, loop spans, the settings used). The **intra fraction is the
primary output** — it is the primary-loop estimate that feeds real elastic
network theory. The network is checkpointed too (`crosslink_state.json`), so a
restarted run resumes with the bonds it had already made instead of silently
building a different network.

`sim crosslink <sim> --distance 0.8` applies the same rule *post-hoc* to a
finished trajectory, where a formed bond changes nothing. Run it on a
non-reactive trajectory and compare: it should count **at least as many**
crosslinks as the reactive run of the same system, and more wherever lysines
have competing partners, because it counts encounters a real network would have
prevented by tethering the chains at the first bond. The gap between the two is
the interesting number, and a post-hoc conversion near 100% mostly says the
cutoff was generous and nothing was ever restrained.

**Read before believing any of it** (the full list is in `tools/crosslink.py`):

- CALVADOS has no crosslinker, no solvent and no activation barrier. The
  reaction is a pure distance criterion on coarse-grained beads. The *sequence
  and topology* of events is meaningful; the *rate* is not — do not report event
  times as kinetics.
- Stock CALVADOS2 has no temperature-dependent hydrophobicity, so with TEA off
  (the default) it cannot reproduce the LCST transition that drives
  coacervation in real ELP crosslinking. If the transition matters, turn on
  [CALVADOS-TEA](#temperature-and-calvados-tea-opt-in).
- The reaction distance dominates every result. Sweep it (0.5, 0.8, 1.2 nm) and
  report the sensitivity rather than one number.
- Contacts in the first frames are artefacts of the initial placement, not
  encounters the dynamics produced. Set `--crosslink-start-step` past
  equilibration (the notebook's Rg/RMSD check says where that is) and state what
  you used.
- Valence changes the answer: Kawamoto et al. (Macromolecules 2015, 48, 8980)
  found A2+B4 networks far more prone to cyclic defects than A2+B3.

## Temperature and CALVADOS-TEA (opt-in)

**Temperature** is the thermostat temperature in **kelvin** (default 293.15 K,
i.e. 20 °C). Anything outside 150–500 K is refused, and a value that looks like
Celsius ("20") is named as such instead of queueing a 20 K run. On its own it
only changes kT: stock CALVADOS2 gives every residue a constant stickiness
(lambda), so a warmer run does **not** move the LCST.

**CALVADOS-TEA** (`tools/tea.py`) makes lambda a function of temperature, so
collapse on heating comes out of the force field. It implements Chen & Zeng
(2026), *JACS Au* 6(8), 4478–4491,
[doi:10.1021/jacsau.6c00523](https://doi.org/10.1021/jacsau.6c00523): hydration
free energies are extrapolated from 300 K, mapped to an excess pair free energy,
and scaled by **gamma** into a lambda shift. It is off by default, and with it off
a run is bit-for-bit the stock pipeline. With it on, `prepare.py` writes a
temperature-specific residues table into `runtime/` and never touches
`residues_CALVADOS2.csv`.

```bash
sim new ... --temperature 350 --tea                  # gamma 6 (the default)
sim new ... --temperature 285 --tea --tea-gamma 3
```

In a batch CSV the same settings are columns, with fallbacks in `Shared(...)`
(`temperature`, `tea`, `tea_gamma`):

| Column | Also accepted as | Values |
| --- | --- | --- |
| `temperature` | `temp`, `T`, `temperature (K)`, `thermostat` | kelvin |
| `tea` | `calvados tea`, `temperature dependent` | `yes`/`no`, `1`/`0`, `on`/`off` |
| `tea gamma` | `gamma` | positive number |

Read before using it:

- TEA only works between **280 and 380 K**, the validated range of the
  underlying hydration fits. A TEA row outside it is refused.
- **gamma = 6 is our own choice**, outside the {2, 3} the paper validated on
  CALVADOS. Runs before 2026-10-05 used 3, and `metadata.csv` records which one
  a run used.
- The paper's printed equation 6 has a sign that contradicts the rest of the
  paper. The code uses the opposite sign, which reproduces the check values in
  `tests/test_tea.py`. This has not been cross-checked against the authors' own
  code. `tools/tea.py` explains the reasoning.
- Some residues leave lambda's usual [0, 1] range (Trp, Tyr and Phe go above 1
  when warm; Glu goes slightly below 0). CALVADOS does not clamp lambda, so the
  module docstring lists where this happens.

The TEA series behind our results (single chains, brushes, lysine and RGD
arrangements, temperature sweeps) have their own scripts. See
[Analysis scripts](ANALYSIS.md#analysis-scripts).

## Ionic strength

CALVADOS screens charges with Debye–Hückel, set by the ionic strength in
**molar**. The default is 0.19 M, which every run used until the column existed.
In a batch CSV, set it per row with an `ionic strength (M)` column (also accepted
as `ionic` or `salt`). Blank cells fall back to `ionic` in `Shared(...)`, and
negative values are refused. `batches/Example-csv-ionic.csv` is the template with this
column. There is no `sim new` flag for it.

## Attaching through lysines (`free` and `preattached` modes)

By default (**brush** mode) residue 0 of every chain is tagged `Z`, given zero
mass and held fixed on its lattice point. Two more modes attach chains through
their **lysines** instead, with residue 0 an ordinary residue and nothing but
the surface-bonded lysines holding a chain:

| mode | at t = 0 | during the run |
| --- | --- | --- |
| `brush` | residue 0 fixed on the lattice | — |
| `free` | every chain has **one** lysine, picked at random, bonded to the surface | free lysines that come within `surface_distance` of the surface bond, with probability `surface_prob` per check, until `surface_max_fraction` of *all* lysines are bonded |
| `preattached` | a random `preattached_fraction` of **all** lysines is bonded, every chain getting at least one | nothing, unless `surface_dynamic=True` |

A surface bond is a stiff harmonic tether from the lysine bead to a fixed point
on the *tether plane* `z = z_wall + surface_tether`, so it holds x, y **and** z
(a bound lysine does not slide), and it never breaks. Like the K–K crosslinks,
every lysine gets a dormant tether up front and binding means ramping its
stiffness in over `surface_ramp_steps`. A lysine is in exactly one of three
states — free, surface-bonded, or K–K-bonded — and the two reactors exclude
each other's lysines. `surface_distance` is its own knob, separate from
`crosslink_distance`, and like it dominates the result: sweep it.

The cap is global. It has to leave room for one bond per chain, so
`surface_max_fraction` (and `preattached_fraction`) must be at least
`chains / lysines`; the planning cell and `sim new` refuse anything less and
say what the minimum is.

Starting geometries are **constructions**: each chain is a stack of vertical
hairpins on a cubic grid of the bond length, standing on the tether plane in
its own lattice cell, every pinned lysine exactly on its pin, no bead below the
plane, no two beads closer than a bond. `run.py` prints the potential energy
before and after minimisation, the energy in the tethers (≈ 0) and how far the
pins moved, and refuses to step if anything is off. Equilibrate before
sampling, and set `surface_start_step` past that point — the construction
leaves beads on the tether plane that would otherwise bind at the first check.

Both modes share the brush lattice, so the pins are spread over the surface
the way the chains are; within one chain the pins sit inside that chain's
footprint.

Everything the modes need lives in the notebook's `Shared` block as
`surface_*` settings, siblings of `crosslink_*` (all in nanometres; the
startup line prints ångström too):

```bash
sim new ... --mode free --surface-distance 0.8 --surface-tether 0.6 --surface-max-fraction 0.5
sim new ... --mode preattached --surface-preattached-fraction 0.4
```

The event record is shared with the K–K crosslinks: `crosslink_events.csv` gains
rows with `kind = surface` (`resid_j`, `chain_j`, `span`, `pymol_j` blank;
`x,y,z` is the pin) and a new last column `origin` — `initial` for a bond the
run started with, `run` for one it formed. `crosslink_summary.txt` and
`metadata.csv` report how many lysines are surface-bonded (count and fraction),
how many chains have one, two or three-plus bonds — the number that says
whether the result is a coating or loosely tethered chains — and, for `free`
mode, when the cap was reached. The notebook's last cell plots the same. The
per-lysine state is checkpointed in `surface_state.json`, so a restarted run
resumes with its bonds; a restart without that file is refused. `sim crosslink`
(post-hoc) keeps residue 0 as a site in these modes and keeps surface-bonded
lysines out of the pair pool from the step they bound.

Read `tools/surface.py` before believing any of it: no crosslinker, no
solvent, no activation barrier (the sequence of events is meaningful, the rate
is not); no LCST transition in CALVADOS2; the wall is purely repulsive unless
`surface_attraction` is set, so dynamic binding in `free` mode is set by
diffusion and `surface_distance`; and the tether height is a modelling choice
with no default that is right by default.

**Brush mode is untouched.** With `mode` blank or `brush` the run takes the
same code path as before: no `surface.yaml`, no extra force, the trajectory is
coordinate-for-coordinate what it was — verified by running the same seeded
system (with and without K–K crosslinking) through the old and new pipelines.
