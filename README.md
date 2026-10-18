# iGEM 2026 TU Delft — ELP Simulations

CALVADOS / OpenMM simulations of elastin-like polypeptides (ELPs).

## Setup

This project is managed with [`uv`](https://docs.astral.sh/uv/). Install it
by following the [official instructions](https://docs.astral.sh/uv/getting-started/installation/),
then clone the repo and run `uv sync`.

### Personal machine

```bash
git clone git@github.com:Delft-iGEM/ELP-Simulations.git
cd ELP-Simulations
uv sync
```

In VSCode, select `.venv/bin/python` as your interpreter. From a terminal,
activate the env with:

```bash
source .venv/bin/activate
```

### DelftBlue

Install `uv` once on a login node by following the
[official instructions](https://docs.astral.sh/uv/getting-started/installation/).

Then add the following to your `~/.bashrc` to redirect `uv`'s caches off the
small `/home` quota:

```bash
export UV_CACHE_DIR=/scratch/$USER/.uv-cache
export UV_PYTHON_INSTALL_DIR=/scratch/$USER/.uv-python
mkdir -p "$UV_CACHE_DIR" "$UV_PYTHON_INSTALL_DIR"
export PATH="$HOME/.local/bin:$PATH"
```

Reload your shell (`source ~/.bashrc`), then clone the repo to `/scratch`
and sync. This must be done **from a login node** — compute nodes have no
internet access:

```bash
cd /scratch/$USER
git clone git@github.com:Delft-iGEM/ELP-Simulations.git
cd ELP-Simulations
uv sync --frozen
```

### Where simulation data goes

Trajectories are big — a few hundred MB per simulation — and DelftBlue caps
`/home` at 30 GB while `/scratch` gives you 5 TB. So on DelftBlue the
`runtime/` folder is a **symlink**: the data itself is written to

```
/scratch/$USER/elp-data/oefeningen/<sim_name>/
```

`sim prepare` and `sim new` create the link for you, and everything else —
`run.py`, `job.sh`, `analyze.ipynb` — keeps using the familiar
`simulations/<sim_name>/runtime/` path, so nothing changes in day-to-day use.
`sim clean` deletes the data on `/scratch` along with the link.

Set `ELP_DATA_DIR` to put the data somewhere else. On a machine with no
`/scratch/$USER` (your laptop), `runtime/` stays a plain folder inside the repo
and nothing is symlinked.

Simulations created before this was set up keep their data in the repo until
you move it:

```bash
sim migrate <simulation>   # or --all
```

## Writing simulations

The simulations use [CALVADOS](https://github.com/KULL-Centre/CALVADOS) on
top of OpenMM. If you haven't worked with CALVADOS before, skim the examples
in their repo first.

To create a new simulation, copy an existing folder under `simulations/`
and rename it. The folder name should match `sim_name` and be descriptive —
e.g. `i70-surfaceattached` for 70 `VPGIG` repeats, surface-attached. You can
use `tools.elibpy.humanize_seq` for naming inspiration.

Inside the simulation folder, `prepare.py` defines the usual CALVADOS
`Config` and `Components` setup; edit parameters there. ELP sequences can
be built with:

```python
from tools.elibpy import build_sequence_with_features
```

You don't have to use `elibpy` — plain strings in the `sequences` dict
work fine too.

For anything that needs to happen on the OpenMM `System` itself (custom
restraints, extra forces, etc.), use the `build_sim` hook in `prepare.py`.
It runs after `sim.build_system()` and before `sim.simulate()`, so the
system exists and is ready to be modified.

The `runtime/` folder is generated and not committed. You **should NOT**
need to edit anything inside it by hand: `prepare.py` regenerates it
idempotently, so running prepare on any machine produces the exact same result.

The Analyze section of `analyze.ipynb` writes its results — the equilibration
verdict, the crosslinking counts, every plot — into `analysis/`, next to
`prepare.py` (`simulations/<sim_name>/analysis/`). Unlike `runtime/`, `analysis/`
**is** committed: it's small (text + PNGs), so pushing it is how a DelftBlue
checkout gets the same numbers and figures without re-running the notebook there.

## A batch of runs from one CSV

`analyze.ipynb`'s first section plans a whole **set** of simulations at once:
one CSV row per run, one job per row. The CSV lives in the project root
(`Example-csv.csv` is the template) and holds only what differs between runs:

```
Name ;Sequence;concentration;nmolecules;steps;crosslink distance;crosslink prob;mode;preattached fraction;walltime
example-v20-dilute;V20;0.02;16;7070000;;;;;0:40:00
example-block-xl;(GAGAGAGA) K1 V3 (VPGFG) V4;0.02;25;7070000;0.8;1.0;;;2:00:00
example-k4-free;K1 V4 K1 V4 K1 V4 K1 V4;0.02;16;7070000;;;free;;1:00:00
example-k4-pre;K1 V4 K1 V4 K1 V4 K1 V4;0.02;16;7070000;;;preattached;0.4;1:00:00
```

* **Any column may be missing, and any cell may be blank** — it then falls back
  to the value set once in the planning cell (`shared = Shared(...)`), which is
  also where everything that never varies within a batch lives: the margin, the
  length measure, the crosslink bond parameters, the SLURM partition.
* **`concentration` is always chains/nm²** here. The other two spacing modes
  (fraction of the chain's own length, grafted mass per area) are still in `sim
  new` and `create_simulation`, but a batch mixing modes wouldn't be comparable
  row to row.
* **`Sequence`** is either a raw single-letter sequence or the `tools/elibpy.py`
  shorthand: `V20` is 20 × `VPGVG`, `(VPGFG)` is inserted literally, so
  `K1 V3 (VPGFG) V4` is one lysine block, 3 × `VPGVG`, one `VPGFG`, 4 × `VPGVG`.
* **`walltime`** is per run (SLURM `HH:MM:SS`), so a long sequence and a short
  one can share a CSV without either waiting for the other's limit.
* **`mode`** is how the chains attach to the surface — blank or `brush` (the
  default, unchanged pipeline), `free` or `preattached` — see
  [Attaching through lysines](#attaching-through-lysines-free-and-preattached-modes).
  `preattached fraction` is only read for `preattached` rows (blank falls back to
  `surface_preattached_fraction` in `Shared`), so one CSV can sweep it.

The planning cell prints a comparison table of every run — box, spacing,
density, steps, frames, walltime — and, as 0/1 opt-ins, the full stats block
(`SHOW_STATS`) and the to-scale space-filling plot of the grafting lattice
(`SHOW_PLOTS`) for each one. Put a `0` behind either to skip it for a long CSV.

The table's footer adds up the walltimes two ways: the sum (what one job running
the rows back-to-back would have to ask for) and the maximum (the longest of
several separate jobs). Compare those against what earlier runs actually cost —
see `run_wall_seconds` / `ns_per_hour` under [Run metadata](#run-metadata).

### Names, and the list of used ones

`used-runs.csv` in the project root is the running list of every name and
sequence that has been submitted. The submit cell checks the whole CSV against
it **before preparing or queueing anything** — a name is refused if it repeats
inside the CSV, already exists as a folder under `simulations/`, or is already
in `used-runs.csv` — and appends the runs to it once `sbatch` has accepted them.
Nothing is recorded for a batch you only prepared, or for the local smoke-test
cell: the list is jobs that actually went out.

A reused *sequence* under a new name is reported as a note, not an error —
running the same ELP again at another concentration is a legitimate thing to do.

```bash
sim used              # show the list
sim used --backfill   # add any simulation on disk that isn't in it yet
```

`sim submit <name>` appends to the same list, so a one-off submission outside the
notebook is just as taken as one from a CSV.

## Preparing and running simulations

After `uv sync`, the `sim` CLI is available inside the venv. Activate it once:

```bash
source .venv/bin/activate   # Linux / macOS
.venv\Scripts\activate      # Windows (PowerShell)
```

Then use `sim` for all simulation tasks:

```bash
sim list                          # list available simulations
sim prepare i70-surfaceattached   # generate runtime/ folder
sim run     i70-surfaceattached   # run locally
sim submit  i70-surfaceattached   # sbatch on DelftBlue
sim clean   i70-surfaceattached   # delete its runtime data
sim migrate i70-surfaceattached   # move existing runtime data onto /scratch
sim metadata i70-surfaceattached  # write runtime/metadata.csv (--all for every simulation)
sim crosslink i70-surfaceattached --distance 0.8  # post-hoc crosslink detection
sim used                          # names/sequences already submitted (--backfill to fill gaps)
```

### Spacing the chains

`sim new` (and the notebook's planning cell) sets the distance between grafting
points in exactly one of three ways:

```bash
sim new ... --concentration 0.02       # chains/nm^2 — same chain count per area
sim new ... --spacing-fraction 0.5     # fraction of the chain's own length
sim new ... --mass-concentration 0.2   # ug/cm^2 — same mass of polymer per area
```

The last one is the one to use when comparing ELPs of different lengths against
an experiment that measures adsorbed mass: a chain twice as long weighs twice as
much, so it gets grafted sqrt(2) further apart and the surface still carries the
same µg/cm². Chain molar mass is summed from the residue table (the anchor bead
counts as the valine it is a tagged copy of, not as its MW = -2 simulation
placeholder). Whichever mode is used, all three numbers — chains/nm², fraction
of the chain length, and µg/cm² — are reported and recorded.

### Space around the lattice

By default the lattice fills the box exactly (`nx × spacing` by `ny × spacing`),
so the periodic images continue it without a seam. A rim chain is half a spacing
from the wall and half a spacing from the next box on the other side — exactly
one spacing from its neighbour, as crowded as a chain in the middle. That is a
uniform infinite grafted surface, and it is the right default.

The price is that chains genuinely sit at the boundary, lean across it, and are
written to the trajectory wrapped to the opposite side. The dynamics are correct
(and the analysis unwraps them), but in a rendered trajectory it looks like the
molecule teleporting. `--margin-fraction` leaves empty surface around the
lattice so that happens less:

```bash
sim new ...                         # default: no rim, lattice tiles the box exactly
sim new ... --margin-fraction 0.5   # half a spacing of rim on each side
sim new ... --margin-fraction 1.5   # more room still, at the cost of a bigger box
```

The rim is `margin_fraction × spacing` wide on every side, so the box becomes
`nx × spacing + 2 × margin`. At 0.5 an edge chain starts a full spacing from the
wall instead of half of one, which stops most of the wrapping.

The trade-off is exactly the uniformity above: with a rim, neighbours across the
boundary are `spacing + 2 × margin` apart instead of `spacing`, so the rim chains
really are less crowded than the ones in the middle and the run is a finite
grafted patch repeated periodically. Worth it for per-chain observables (Rg,
height, RMSD), where a dilute rim barely matters and a readable trajectory does;
not worth it for anything that depends on uniform lateral crowding. `spacing`
means the nearest-neighbour distance either way, so the three spacing knobs keep
their meaning; what a margin dilutes is the *box-averaged* density, reported
separately as `box_concentration_chains_per_nm2` and
`box_mass_concentration_ug_cm2`.

### Box height

The box is periodic in **z as well as x and y**. CALVADOS builds every
nonbonded force with the periodic cutoff method, and the "surface" is a
one-sided potential at low z, not a lid — so a bead that comes within the
electrostatics cutoff (4 nm) of the box top interacts with the periodic image
of the grafted layer at the bottom, as if it were next to the surface again.
The box must stay taller than the brush ever reaches, plus that cutoff.

The old fixed 50 nm was not enough: the 720–1040-residue block ELPs at
0.02 chains/nm² reached 51–63 nm (highest bead over the run, 0.16–0.19 × the
contour length). The height is now planned per run:

* default = 0.3 × the contour length, rounded up to 5 nm, never below 50 nm
  (`plan_box_height` in `tools/new_simulation.py`);
* override with `box_height` in the notebook's `Shared`, or `sim new ... --box-height 120`.

A taller box costs nothing to speak of — the neighbour list scales with the
beads, not the empty volume above them. After a run, `metadata.csv` records
`z_max_nm`, `box_headroom_nm` and `box_tall_enough`; the notebook's load cell
prints the same check. Headroom under the cutoff means the run touched its
own image: raise the height and rerun.

Wrapping in x/y is a different thing and harmless: a chain that leans across a
lateral edge is written translated by one box length, but the forces never see
that jump — with periodic boundaries the chain's image *is* its neighbour, at
the same grafting density. Contact counts are unaffected too, because the
notebook and the live crosslinker both measure minimum-image distances.

### Reactive crosslinking (opt-in)

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

### Temperature and CALVADOS-TEA (opt-in)

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
[Analysis scripts](#analysis-scripts).

### Ionic strength

CALVADOS screens charges with Debye–Hückel, set by the ionic strength in
**molar**. The default is 0.19 M, which every run used until the column existed.
In a batch CSV, set it per row with an `ionic strength (M)` column (also accepted
as `ionic` or `salt`). Blank cells fall back to `ionic` in `Shared(...)`, and
negative values are refused. `Example-csv-ionic.csv` is the template with this
column. There is no `sim new` flag for it.

### Attaching through lysines (`free` and `preattached` modes)

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

### Run metadata

Every prepared run gets a `runtime/metadata.csv` next to its `.dcd` — the
settings that define it (sequence, chain length, nmol, spacing and the fraction
of the chain length it corresponds to, the margin around the lattice, surface
concentration, box, the attachment `mode`, any crosslinking or surface-attachment
settings and how many bonds formed, total steps,
steps per frame, simulated time, temperature, ionic strength, pH, wall depth,
cutoffs, force field, platform, package versions) plus how far the run actually
got and **how long it took**. `prepare.py` writes it, `run.py` refreshes it when the run finishes, and
`sim metadata <sim>` (or `--all`) regenerates it — including for simulations
that predate the file, since everything is recovered from `config.yaml`,
`components.yaml`, `molecules.fasta`, `lattice.yaml`, `crosslink.yaml`,
`surface.yaml` and the DCD header. (Only the lattice has to be stated rather than derived: once there is a margin around
it, the box is one number holding two unknowns. Runs from before margins existed
have no `lattice.yaml` and are read back as the seamless lattices they are.)

#### How long a run takes

Deciding between one long job and several short ones needs the wall-clock cost
of a simulation, which nothing CALVADOS writes records — so `run.py` times
`simulate()` itself into `runtime/timing.yaml`, and `metadata.csv` reports:

| key | what it says |
| --- | --- |
| `walltime_requested` (+`_hours`) | the SLURM `--time` the job was given, read back out of `job.sh` |
| `run_wall_seconds`, `run_wall_clock` | what the MD actually took, summed over every leg |
| `run_legs` | how many jobs it took — >1 means it hit the limit and was resumed from its checkpoint |
| `walltime_used_fraction` | measured / requested, i.e. how much of the request was slack (>1 = killed at the limit) |
| `ns_per_hour`, `steps_per_second` | the rate, from the frames actually written, so it is meaningful mid-run too |
| `run_started_utc`, `run_finished_utc`, `partition`, `cpus_per_task` | when, and on what |

Only `simulate()` is inside the clock — building the system, and waiting in
DelftBlue's queue, are not what an extra sequence would cost. `walltime_requested`
and `partition` come from `job.sh`, so they are filled in for every run that has
one, including older ones; the measured fields stay empty for runs that finished
before they were recorded.

One row per setting (`key,value,unit,description`), so a batch of runs is one
concat away from a comparison table:

```python
import pandas as pd
from pathlib import Path
runs = {p.parent.parent.name: pd.read_csv(p, index_col="key")["value"]
        for p in Path("simulations").glob("*/runtime/metadata.csv")}
pd.DataFrame(runs).T
```

### Is it equilibrated?

The first cell of the Analyze section in `analyze.ipynb` plots the radius of
gyration and the RMSD of the grafted chains, each per frame and cumulatively,
works out how much burn-in to discard and how many independent samples the rest
is worth, and from that estimates how many steps the next run of a comparable
sequence actually needs. Both observables are there because they fail
differently: a chain reaches its final size long before it has forgotten the
conformation it started in.

### What the run actually crosslinked

Two different questions, two cells, and they do not give the same answer:

* **Contacts** — how often two crosslinkable residues came within
  `contact_cutoff` of each other, split intra-chain vs inter-chain, per pair and
  in total. That is the *opportunity* to react, and it is defined even for a run
  with no chemistry switched on.
* **Bonds formed** — what the run's reactor actually committed to, read from
  `runtime/crosslink_events.csv`. K–K crosslinks split into **intra-chain**
  (loops inside one polymer) and **inter-chain** (bridges, the ones that build a
  network), lysine–surface bonds split into the ones present at t = 0 and the
  ones formed during the run, how many of the system's lysines that consumed and
  how many are still free, plus how many chains ended up bridged to another and
  the largest connected group. Four plots: cumulative bonds over time, formation
  rate, the lysine budget, and network growth. A run with no
  `crosslink_events.csv` says so and plots nothing.

Bond *times* are not kinetics — see the caveats in `tools/crosslink.py` and
`tools/surface.py`. The counts and the intra/inter split are the meaningful part.

### Analyzing a whole CSV at once

The Analyze section is pointed at its target by `ANALYZE`, at the top of the
section. It takes one simulation name, a runs CSV (`block-runs.csv`,
`used-runs.csv`, ...), a pattern like `triblock-64-*`, or a list mixing them:

```python
ANALYZE = "block-runs.csv"
sim_names = resolve_targets(ANALYZE)
```

Names out of a CSV are slugified the same way `create_simulation` slugifies them
when it makes the folder, so `My_ELP 2` finds `simulations/my-elp-2/`. A name
that still doesn't match — a row renamed after it was planned, one that was
never submitted, one whose job hasn't written a trajectory yet — is printed with
the closest folder that *does* exist and dropped:

```
19 of 20 run(s) from block-runs.csv can be analyzed
   !  skipping monoblock-free-xl2: no simulations/monoblock-free-xl2/ folder — closest folder that does exist: monoblock-free-xl
   !  skipping tetrablock-pre-xl: no .dcd in runtime/ yet — still queued or running?
```

One missing ELP never costs you the rest of the batch. Work through the section
once on the first run, then the **last cell of the section** re-runs every cell
tagged `analysis` for each of the others in turn, each writing into its own
`simulations/<name>/analysis/`. The cells are read from the notebook *on disk*,
so save it after editing one. A run whose cells raise is reported and abandoned
— the batch carries on with the next name rather than analyzing it with the
previous run's trajectory still loaded — and a run that simply has nothing to
analyze (a sequence with no lysines, for the crosslink cells) is noted, not
counted as a failure. `tools/analysis_batch.py` has the machinery.

### Where one residue sits in z (`sim distribution`)

`sim distribution` plots the distribution of one amino acid's height above the
surface (z), pooled over every frame of a trajectory. For example, it shows where
the lysines end up in a brush:

```bash
sim distribution --pdb simulations/<sim>/runtime/top.pdb \
                 --dcd simulations/<sim>/runtime/<sim>.dcd --residue K
```

`--residue` takes a one-letter (`K`) or three-letter (`LYS`) code. `--bins` sets
the histogram bins (default 60), `-o` / `--output` sets where the plot goes (by
default next to the DCD), and `--show` also opens the plot interactively.

## Analysis scripts

The notebook covers one run or one CSV of runs. The series behind our results
have their own scripts in `tools/`, grouped below by the result they produce.
Run each one from the project root:

```bash
uv run python tools/<script>.py           # most also accept --help
```

The scripts follow two patterns:

- **Run generators** (`plan`, `init`, `submit --yes`): `plan` prints the runs,
  `init` writes `simulations/<name>/prepare.py` for each run and prepares it,
  and `submit --yes` sends them to SLURM.
- **Analysis scripts**: run with no arguments to analyse every finished run of
  the series, which writes a CSV and figures into `figures/<series>/`.
  `--plot-only` redraws the figures from that CSV without touching
  trajectories.

`tea_sweep.py` is the exception: `init` writes a signac-style
`workspace/<job_id>/` instead of `simulations/`.

`figures/` is not part of this repository, so run an analysis once before using
`--plot-only`. `plot_sasa_24.py` and several of the RGD scripts also read
`figures/free-rgdarr/sasa-24-arrangements.csv`, which no script in `tools/`
writes and which is not in this repository.

Some generators take their sequences from earlier runs, so `simulations/` keeps
just the `prepare.py` of those source runs (no results):

- `simulations/0r0g0d1-2-4-7-focused/`: the original 357-aa RGD construct,
  read by `tea_sequences.py`. That module is used by `tea_temps`, `tea_g6`,
  `tea_lysramp`, `tea_random` and `tea_tsweep`.
- `simulations/rgdarr01-*` … `rgdarr24-*`: read by `rgdarr_graft.py`.
- `simulations/free-rgdarr01-*` … `free-rgdarr24-*`: read by `tea_temps.py`,
  `tea_g6.py` and `tea_pairs.py`.
- `simulations/lysarr01/02/03/09/12/13/19/20-*`: the best and worst lysine
  arrangements, read by `tea_g6.py` and `tea_pairs.py`.
- `simulations/tea-{lys,rgd}-{best,worst}-free-350k/`: `tea_temps.py` takes
  the seeds for its no-TEA controls from these runs' `runtime/config.yaml`.
  Run `tea_pairs.py init` before `tea_temps.py init`.

Every other run folder is in the
[full repository](https://github.com/Delft-iGEM/ELP-Simulations).

### Temperature response (CALVADOS-TEA)

| Series | Generate runs | Analyse / plot |
| --- | --- | --- |
| (VPGXG)₈₄ sequences × temperatures | `tea_sweep.py init` / `submit` / `status` | `tea_analysis.py` (Rg, collapse fraction, scaling exponent ν) |
| 270 / 350 K set and no-TEA controls | `tea_temps.py` (sequences from `tea_sequences.py`) | `plot_tea_temps.py --rg-csv <csv>` → `figures/tt/` |
| γ = 6, 285 / 370 K single chains | `tea_g6.py` | `tea_g6_analysis.py`, `plot_tea_g6.py`, `tea_g6_states.py` (collapsed vs extended) → `figures/g6/` |
| Temperature sweep, 275–380 K in 5 K steps | `tea_tsweep.py` | `tea_tsweep_analysis.py` (thermoresponse), `tea_tsweep_sasa.py` (exposure) → `figures/g6-tsw/` |
| Random arrangements, 285 / 370 K | `tea_random.py` | `tea_random_analysis.py` (Rg, collapsed-state dwells) → `figures/g6-random/` |
| Best/worst pairs, brush vs free chain, 280 / 350 K | `tea_pairs.py` | `tea_pairs_sasa.py` → `figures/tea-pairs/` |

### Lysine placement and amount

| Series | Generate runs | Analyse / plot |
| --- | --- | --- |
| Lysine fraction, (VPGIG)₉₆ with 0–100 % VPGKG | `tea_lysfrac.py` | `tea_lysfrac_analysis.py` → `figures/g6-kfrac/` |
| Lysine ramp, 4 → 49 VPGKG | `tea_lysramp.py` | `tea_lysramp_analysis.py`, `tea_lysramp_dwells.py`, `plot_lysramp_290k.py` → `figures/g6-kspc-kend/` |
| 8-lysine arrangements (g6-k8-*) | `tea_g6.py` | `tea_k8_dwells.py` (collapsed-state dwells) → `figures/g6/` |
| Submodule 4: lysine contacts, Rg and SASA on the surface | `sub4_mains.py`, or the `*sub4*` CSVs through the notebook | `sub4_analysis.py` (`--runs 'sub4-extra2-*'` for a subset) → `figures/sub4/` |
| Lysine-arrangement runs (`lysarr*`) | `lysine-arrangement-sweep.csv` through the notebook | `build_lysine_workbook.py` → `lysine-arrangement-results.xlsx` |

### RGD exposure (functional motif accessibility)

| Series | Generate runs | Analyse / plot |
| --- | --- | --- |
| 24 RGD arrangements (`rgdarr*`, preattached) | `rgd-arrangement-sweep.csv` through the notebook | the notebook's Analyze section |
| The same 24 as free single chains (`free-rgdarr*`) | — (no generator script in this repository) | `plot_sasa_24.py`, `rgd_distance.py` (mean distance between the 4 RGD sites) → `figures/free-rgdarr/` |
| The same 24, grafted at 0.03 chains/nm² for 500 ns (`rgdg*`) | `rgdarr_graft.py` | `rgdg_analysis.py` (SASA and Rg) → `figures/rgdg/` |
| RGD placement, collapsed vs extended (g6-rgd*) | `tea_g6.py` | `g6_rgd_states.py` → `figures/g6/` |
| Grafted vs free exposure | — | `build_sasa_workbook.py` → `sasa-grafted-vs-free.xlsx` |
| Two designs head to head | — | `sasa_two_sample.py <run-a> <run-b>` (Welch t-test across the 16 chains) |

`dwell_slides.py lysramp` / `dwell_slides.py rgd` redraws the dwell figures as
single-panel slide versions.

### Shared building blocks

These modules have no command of their own. The notebook and the scripts above
import them:

- `sasa.py`: solvent-accessible surface area of a motif for coarse-grained
  beads. It uses CALVADOS bead sizes, because standard SASA tools assume atoms.
- `hydrodynamic.py`: hydrodynamic radius (Rh, what dynamic light scattering
  measures) next to Rg.
- `network.py`: network topology of a crosslinked coating, built from
  `crosslink_events.csv`.
- `analysis_workbook.py`: one spreadsheet (contacts, bonds, equilibration)
  from a batch of finished runs' `analysis/` folders.
- `analysis_io.py`, `paths.py`: where analysis output and trajectory data live.
- `elibpy.py`: the ELP sequence shorthand (`V20`, `K1 V4`, ...).

### Shell autocomplete

Install completion once (detects your shell automatically — bash, zsh, fish, or PowerShell):

```bash
sim --install-completion
```

Restart your shell (or `source ~/.bashrc`), then **Tab** after `sim prepare ` or
`sim run ` to autocomplete simulation names.

To see the raw completion script without installing it:

```bash
sim --show-completion
```

### Running without the venv active

Prefix any command with `uv run`:

```bash
uv run sim prepare i70-surfaceattached
uv run sim run     i70-surfaceattached
```
