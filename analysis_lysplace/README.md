# RGD+ lysine placement — which arrangement gives the most inter-chain crosslinks

Sub-submodule of the lysine work. Written 2026-10-16. Everything here is derived; the raw
simulation output under `simulations/*/runtime` was only read.

**In one paragraph.** 130 runs of the RGD+ construct (16 chains, 4 K each, crosslinking on, 400 ns)
were compared across four lysine arrangements. Crosslinking saturates: 19–22 of the 22 possible bonds form in
almost every run, so the *number* of crosslinks says little. What the arrangement decides is
**where the bonds go**. Inter-chain bridges per run: clustered 3, focused 7, random 14, uniform 19
(medians). Every class differs from every other (Holm p < 1e-4). One variable explains it: **the sequence distance
from a lysine to the nearest other lysine on its chain**. Close neighbours react with each other
as loops; far-apart lysines bridge to other chains. Each e-fold increase in that distance raises the
odds of a bridge 3.4-fold. Position along the chain does not matter. A per-lysine model built on
this predicts held-out runs with R² 0.88 and points to the optimum: spread the four lysines as
evenly as the construct allows (minimum spacing ≈ 100 residues, about 20 of 22 bridges). Any design with
minimum spacing ≥ 60 residues comes within ~1 bridge of that. Placement also sets **where the network
sits**: one-block designs crosslink in a thin layer at the wall, while spread designs bridge mid-brush.
It barely changes brush height (h95 within 1.6 nm) and does not change RGD exposure (fraction above h50).

---

## 1. Data — `run_inventory.csv`

| sheet | runs | classes |
|---|---|---|
| `Chr_RGD+_Lysine_Sweep.csv` | 40 | 10 uniform, 10 focused, 10 clustered, 10 random |
| `Chr_RGD+_Lysine_sweep_20.csv` | 20 | 7 uniform, 7 focused, 6 random |
| `Chr_RGD+_Lysine_sweep_30.csv` | 30 | 7 clustered, 8 uniform, 8 focused, 7 random |
| `Chr_RGD+_lysRandom.csv` | 20 | random |
| `Chr_RGD+_lysine_scramble_sweep_20.csv` | 20 | uniform (no class suffix; same one-K-per-block rule, only labelled differently) |

Totals per class: clustered 17, focused 25, random 43, uniform 45 (25 with the `_uniform` suffix + 20
from the "scramble" sheet). *Until 2026-10-17 the scramble sheet was analysed as a fifth class; see the
dated note in `settings.yaml`.*

All 130 runs passed QC: completed, 5714 frames, a single job, 19 of 64 lysines pinned at t = 0
and none during the run, valence 1, box checks pass. Every run name matched the K positions in its simulated
`molecules.fasta`. Constant across runs: 357 residues (LQ + 4 × [RGD domain + ELP block of 14, 11, 15,
11 VPGXG + VP], M at repeat 6 of blocks 2 and 4), C = 0.04 chains/nm², 0.19 M, 293.15 K,
crosslinks at 0.8 nm, p = 1, from 10 ns. **No design was simulated twice**, so noise is judged only
from the spread within each class.

| class | structure | closest K–K spacing on a chain (residues) |
|---|---|---|
| clustered | 4 K in one block, adjacent repeats | 5 |
| focused | 4 K in one block, spread out | 5–15 |
| random | 19 two-block, 19 three-block, 3 one-per-block, 2 one-block | 5–80 |
| uniform | one K per block | 30–105 |

Class and spacing are almost the same variable. Only the random class covers the whole range.

## 2. Frozen definitions — `settings.yaml`

Written before any network, height, timing or per-lysine quantity was computed. **Already seen
at that point:** the per-class mean bridge counts, from a data check of the event logs. The class
ranking on the primary response is therefore not a blind result; everything else is. Dated
notes at the end of the file record what changed in interpretation afterwards.

* **Primary:** `n_bridges`, the number of inter-chain K–K crosslinks at the end of the run (the stated goal).
* **Network:** loops; crosslinks formed; bridge fraction; distinct chain pairs bridged; secondary loops (extra bridges between
  already-joined chains); largest connected cluster of chains (surface not a node); cycle rank;
  unreacted lysines; network at 100 ns. Timing (t50 / t90 of bridges) is descriptive only.
* **Heights:** production frames 1428–5713 (as in SubM4), heights above the wall in nm. Bridge
  and loop height = midpoint of the two bonded beads, every production frame after the bond formed.
  Brush h50 / h95 are mass-weighted. RGD = the R, G, D beads of the four GRGDS motifs.
* **Optimality rule:** (1) most bridges; (2) a connected network in ≥ 90 % of runs; (3) among ties,
  more distinct pairs bridged, then the smaller spread.

## 3. Results — class comparison (`stats/summary.txt`, `figures/class_comparison`)

Medians per class (n = 17 / 25 / 43 / 45):

| | clustered | focused | random | uniform | Kruskal–Wallis p |
|---|---|---|---|---|---|
| **bridges** | **3** | **7** | **14** | **19** | 2e−21 |
| loops | 16 | 14 | 8 | 3 | |
| crosslinks formed (of 22) | 19 | 21 | 22 | 22 | |
| distinct chain pairs bridged | 3 | 7 | 13 | 18 | |
| largest connected cluster / 16 | 0.13 | 0.19 | 0.69 | 0.94 | |
| network connected (all 16 chains), % of runs | 0 | 0 | 14 | 47 | |
| cycle rank | 0 | 0 | 2 | 5 | |
| bridge height (nm) | 1.5 | 2.2 | 6.1 | 7.2 | |
| bridges above h50 (fraction) | 0.00 | 0.00 | 0.55 | 0.70 | |
| brush h95 (nm) | 14.2 | 13.3 | 13.5 | 12.7 | 2e−4 |
| RGD height (nm) | 6.7 | 6.4 | 6.4 | 6.2 | 0.02 |
| Rg (nm) | 5.75 | 5.56 | 5.44 | 5.25 | |

* **Bridges:** every pair of classes differs (Holm p ≤ 1.3e−5; Cliff's δ −0.75 to −1.0; random vs
  uniform δ −0.75 [−0.88, −0.60]).
* **Crosslinks formed** saturate in random and uniform. In clustered designs 7 lysines
  stay unreacted (median), because their partners on the same chain are used up or pinned.
* **Where the network sits** (`figures/z_profiles`): in clustered and focused designs every crosslink,
  bridge or loop, lies within ~3 nm of the wall. The pinned lysine pulls its whole K block down,
  so all reactions happen in the surface layer (0 % of bridges above h50). In spread designs,
  bridges peak at 6–7 nm, above the mass peak, with ~70 % above h50. The brush itself barely
  changes (h50 5.4–5.7 nm, p = 0.06; h95 within 1.6 nm, but uniform is lower than clustered and
  random, Holm p = 0.001 and 0.003).
* **RGD exposure** is not affected by placement: the fraction of RGD above h50 does not differ
  (p = 0.15, all pairs Holm p > 0.5). The absolute RGD height does differ slightly (p = 0.02), but only
  clustered vs uniform survives Holm (6.75 vs 6.19 nm, p = 0.016), i.e. RGD sits ~0.5 nm lower,
  in step with the slightly lower brush (h95) of uniform designs. Choosing a placement for
  crosslinking costs essentially nothing in RGD presentation.
* **Timing** (`figures/bridge_timecourse`): half of the bridges form within ~5 ns of crosslinking
  starting in every class except clustered (17 ns). The last bonds form up to ~390 ns later, so the
  network is still finishing at the end of some runs; event times are not kinetics.
* **Descriptors** (`stats/descriptor_glm.csv`, quasi-binomial on bridges / 22): log minimum spacing
  (OR 1.9 per e-fold, p < 1e-6) and the number of blocks carrying K (OR 1.5 per block) explain the
  runs with pseudo-R² 0.76 (CV R² 0.80); the K centroid does not matter. Within random, only spacing
  matters (number of blocks p = 0.68). The class label still adds information beyond these three
  run-level descriptors (F(3, 123) = 15.2, p = 2e−8). The per-lysine model (§4) resolves this: what counts is each lysine's
  own nearest neighbour, not the chain's minimum.

### The optimality rule, applied (`stats/optimal.txt`)

Step 2, connectivity in ≥ 90 % of runs, **is met by no class.** Even uniform joins all 16 chains in
only 47 % of runs (21 of 45). With 4 lysines per chain, 19 of 64 pinned, a single fully connected network is
not reliable for any arrangement at this density. The rule is left as frozen and reported as
failing; this is a finding in its own right. On the primary response and every tie-break,
**uniform (one lysine per block) comes first** (median 19 bridges, IQR 3). Adding the 3 random
designs that also happen to have one K per block, the 48 one-per-block runs give a median of 19
bridges, and 48 % are connected.

## 4. How to optimise placement — the per-lysine model (`stats/lysine_model.txt`)

Unit: each of the 5850 free (not pinned) lysines. Logistic GEE within run, with predictors fixed in advance.

| P(bridge) | odds ratio [95 % CI] | p |
|---|---|---|
| log distance to the nearest other K on the chain | **3.43 per e-fold** [2.98, 3.95] | < 1e−16 |
| other K within 30 residues | 0.94 [0.81, 1.09] | 0.42 |
| nearest K is pinned | 1.51 [1.20, 1.90] | 5e−4 |
| position along the chain (linear, quadratic) | n.s. | 0.23, 0.41 |

P(loop) mirrors this (OR 0.35 per e-fold; 0.41 if the nearest K is pinned). `figures/lysine_bridge_probability`:
a lysine 5 residues from its neighbour bridges ~13 % of the time; at 20 residues ~50 %; at 60 ~80 %;
beyond 100 ~90 %. A pinned neighbour helps because it cannot be a loop partner.

* **Check:** summing P(bridge) / 2 per run predicts held-out runs (10-fold CV) with R² 0.88, RMSE 2.2
  bridges. Class means are reproduced within 1 bridge (`figures/lysine_model_check`).
* **Design search:** all 211 876 placements of 4 K on the 49 non-M VPGXG positions, pinning averaged.
  Expected bridges rise steadily with the minimum spacing and level off (`figures/design_rule`):

  | minimum K–K spacing (residues) | 5 | 20 | 40 | 60 | 80 | 100 |
  |---|---|---|---|---|---|---|
  | best expected bridges | 12.2 | 16.1 | 18.4 | 19.3 | 19.8 | 20.0 |

  The top designs (e.g. `1R5G10D11`, `1R3G8D11`, `2R4G9D10`; 20.06 [19.6, 20.5] expected bridges)
  all space the lysines ~100–105 residues apart. That is the most even spacing the 14-11-15-11 block layout
  allows, and inside the data range (one simulated design reaches 105). Their predicted advantage over
  a good uniform design (~19.5) is about half a bridge, smaller than the run-to-run spread (SD 2.5 within uniform), so they are a
  plateau, not a single winner.
* **Design rule:** maximise the smallest sequence distance between any two lysines on a chain;
  aim for ≥ 60 residues (within ~1 bridge of the optimum) and ideally even spacing (~100 residues for
  4 K on 357 residues). Never put two lysines within ~20 residues of each other: below that,
  loops win. Where along the chain the lysines sit, and in which block, does not matter beyond that.

## 5. Caveats

* No replicates: a single design's bridge count carries ~±2 bridges of run-to-run noise. Confirm
  the top designs with 3 replicates each before claiming a specific winner.
* Crosslinking here is first-encounter and irreversible (p = 1). Lower reaction probabilities may favour
  bridges more, since close neighbours would no longer win every race. The ranking should hold;
  the size of the gaps may not.
* "Connected" ignores the surface: chains linked only through their pins count as separate.
* The per-lysine model treats lysines as independent given the predictors. Valence-1 competition is
  absorbed into the coefficients, and the run-level check shows it works over this design space.
* The heights are of a brush on a wall at one density (0.04 chains/nm²). How they carry over to a bulk gel is open.

## 6. Files

* `run_inventory.csv` — design table: K positions, descriptors, metadata, QC
* `runs_metrics.csv` — one row per run: design, network (end and 100 ns), timing, heights, chain size
* `per_lysine.csv` — 8320 lysines: block, repeat, state (pinned / bridge / loop / unreacted), partner,
  formation time, height, model predictors
* `per_bond.csv` — 2750 crosslinks: kind, chains, residues, formation time, height over production
* `profiles.npz` — production z profiles per run (mass, lysines, bridges, loops, RGD), 2 Å bins
* `timeseries/<run>.npz` — the per-run extraction (heights, chain size, per-bond / per-lysine heights, profiles)
* `stats/` — `summary.txt`, `class_summary.csv`, `class_tests.csv`, `pairwise.csv`, `descriptor_glm.csv`,
  `optimal.txt`, `lysine_gee_bridge.csv`, `lysine_gee_loop.csv`, `lysine_model_cv.csv`,
  `design_ranking.csv`, `design_top.csv`, `design_rule.csv`, `lysine_model.txt`
* `figures/` — `class_comparison`, `z_profiles`, `lysine_bridge_probability`, `lysine_model_check`,
  `design_rule`, `bridge_timecourse` (HTML with hover + PNG). Transparent background, mid-grey text
  (reads on white and on near-black), no figure titles: the title goes in the caption

## 7. Reproduce

From the repository root (packages layered by `analysis/scripts/run.sh`, fixed seed 20261016):

```bash
analysis_lysplace/scripts/run_all.sh
# = .venv/bin/python analysis_lysplace/scripts/00_inventory.py
#   sbatch analysis_lysplace/scripts/extract.sbatch     (or 01_extract.py --workers 8; ~12 s per run)
#   analysis/scripts/run.sh analysis_lysplace/scripts/02_metrics.py
#   analysis/scripts/run.sh analysis_lysplace/scripts/03_compare.py
#   analysis/scripts/run.sh analysis_lysplace/scripts/04_lysine_model.py
#   analysis/scripts/run.sh analysis_lysplace/scripts/05_figures.py
```

New runs: add the sheet to `SHEETS` in `scripts/common.py` (the class comes from the name suffix)
and rerun. Step 1 only extracts runs without an `.npz`.
