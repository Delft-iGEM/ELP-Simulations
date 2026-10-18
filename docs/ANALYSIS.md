# Analysis

Analysing runs in `analyze.ipynb`, and the scripts behind each result series.

## Is it equilibrated?

The first cell of the Analyze section in `analyze.ipynb` plots the radius of
gyration and the RMSD of the grafted chains, each per frame and cumulatively,
works out how much burn-in to discard and how many independent samples the rest
is worth, and from that estimates how many steps the next run of a comparable
sequence actually needs. Both observables are there because they fail
differently: a chain reaches its final size long before it has forgotten the
conformation it started in.

## What the run actually crosslinked

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

## Analyzing a whole CSV at once

The Analyze section is pointed at its target by `ANALYZE`, at the top of the
section. It takes one simulation name, a runs CSV (`batches/block-runs.csv`,
`used-runs.csv`, ...), a pattern like `triblock-64-*`, or a list mixing them:

```python
ANALYZE = "batches/block-runs.csv"
sim_names = resolve_targets(ANALYZE)
```

Names out of a CSV are slugified the same way `create_simulation` slugifies them
when it makes the folder, so `My_ELP 2` finds `simulations/my-elp-2/`. A name
that still doesn't match — a row renamed after it was planned, one that was
never submitted, one whose job hasn't written a trajectory yet — is printed with
the closest folder that *does* exist and dropped:

```
19 of 20 run(s) from batches/block-runs.csv can be analyzed
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

## Where one residue sits in z (`sim distribution`)

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
| Lysine-arrangement runs (`lysarr*`) | `batches/lysine-arrangement-sweep.csv` through the notebook | `build_lysine_workbook.py` → `lysine-arrangement-results.xlsx` |

### RGD exposure (functional motif accessibility)

| Series | Generate runs | Analyse / plot |
| --- | --- | --- |
| 24 RGD arrangements (`rgdarr*`, preattached) | `batches/rgd-arrangement-sweep.csv` through the notebook | the notebook's Analyze section |
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
