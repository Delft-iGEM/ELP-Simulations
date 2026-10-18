# Writing and running simulations

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

## Spacing the chains

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

## Space around the lattice

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

## Box height

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

## Run metadata

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

### How long a run takes

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
