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
pymol_i,pymol_j` — `frame` is the saved-trajectory frame the event rounds to,
`chain_*` are 0-based chain indices, `resid_*` are 1-based within their chain,
`x,y,z` is the midpoint of the pair, and `pymol_*` are 1-based for PyMOL's
`index` selector) and `crosslink_summary.txt` (conversion, the
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
- CALVADOS2 has no temperature-dependent hydrophobicity, so it cannot reproduce
  the LCST transition that drives coacervation in real ELP crosslinking. Use
  HPS-T if the transition matters.
- The reaction distance dominates every result. Sweep it (0.5, 0.8, 1.2 nm) and
  report the sensitivity rather than one number.
- Contacts in the first frames are artefacts of the initial placement, not
  encounters the dynamics produced. Set `--crosslink-start-step` past
  equilibration (the notebook's Rg/RMSD check says where that is) and state what
  you used.
- Valence changes the answer: Kawamoto et al. (Macromolecules 2015, 48, 8980)
  found A2+B4 networks far more prone to cyclic defects than A2+B3.

### Run metadata

Every prepared run gets a `runtime/metadata.csv` next to its `.dcd` — the
settings that define it (sequence, chain length, nmol, spacing and the fraction
of the chain length it corresponds to, the margin around the lattice, surface
concentration, box, any crosslinking settings and how many bonds formed, total
steps,
steps per frame, simulated time, temperature, ionic strength, pH, wall depth,
cutoffs, force field, platform, package versions) plus how far the run actually
got. `prepare.py` writes it, `run.py` refreshes it when the run finishes, and
`sim metadata <sim>` (or `--all`) regenerates it — including for simulations
that predate the file, since everything is recovered from `config.yaml`,
`components.yaml`, `molecules.fasta`, `lattice.yaml`, `crosslink.yaml` and the
DCD header. (Only the lattice has to be stated rather than derived: once there is a margin around
it, the box is one number holding two unknowns. Runs from before margins existed
have no `lattice.yaml` and are read back as the seamless lattices they are.)

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
