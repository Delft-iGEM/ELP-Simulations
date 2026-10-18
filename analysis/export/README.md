# SubM4 export for the write-up

Generated 2026-10-17 14:07 by `analysis/scripts/09_export.py`, repository HEAD `6539b98`. **The analysis scripts are not committed** (`analysis/` is untracked), so `MANIFEST.md` lists a SHA-256 for each script and `scripts/` holds a copy.

## What is here

* The finished SubM4 analysis: `stats/`, `figures/`, `runs_metrics.csv`, `run_inventory.csv`, `per_lysine.csv`, `timeseries_primary.csv.gz`, `reactive_anchors.csv`, `settings.yaml`, `qc/equilibration_per_run.csv`, and the full write-up as `analysis_README.md`.
* `models/`: six response-surface fits, each with coefficients and 95 % CIs, fit statistics, grouped CV, lack of fit, residual diagnostics, Cook's distance, the one-parameter-at-a-time table and LMG shares (`comparison.csv` and `lmg_shares.csv` put them side by side). Responses: `inter` = log inter-chain free contacts (frozen primary), `intra` = log intra-chain free contacts, `interfrac` = inter / (inter + intra) as a quasi-binomial GLM with logit link on the production pair-frame counts. Each has a `linear` (frozen) and a `log` coding. `models/inter_linear` reproduces `stats/` exactly.
* `figures_logcoded/`: the frozen figure set redrawn for the log-coded primary model.
* `video_runs.md` and `video_runs/`: the runs picked for the videos, with paths, the production start frame, and per-frame contact and pin files (checked frame by frame against the metrics).
* `MANIFEST.md`: every file with its contents and columns. `analysis_export.zip`: all of it.

## New: the log-coded model

The full quadratic was refitted with the natural log of all four parameters, each coded to −1…+1 over its log range. On the primary response it raises R² from 0.959 to 0.981 and grouped-CV R² from 0.890 to 0.936, cuts the lack of fit from F = 4.94 to 2.10 (still significant, p = 0.011), and removes the spurious density maximum: the density effect is ×8.8 [7.8, 9.9] from 0.01 to 0.06 chains/nm², against ×5.5 in the frozen model, and it rises monotonically.

**The log coding was chosen after the lack of fit of the frozen model had been seen.** It is not part of the frozen analysis plan. The frozen linear-coded model stays the baseline, and any claim based on the log coding should say so.

The most influential run under the log coding is `sub4-extra2-l58-s4-c0p01-i0p43-05` (Cook's D = 1.37, studentised residual -6.1), the low-density, narrow-spacing corner already flagged in the frozen fit; it passes QC and is kept.

| response | coding | n | R² | adj. R² | CV R² | CV RMSE (link) | lack of fit F (df) | p |
|---|---|---|---|---|---|---|---|---|
| inter | linear | 122 | 0.959 | 0.953 | 0.890 | 0.231 | 4.94 (75, 32) | 1.9e-06 |
| inter | log | 122 | 0.981 | 0.978 | 0.936 | 0.176 | 2.10 (75, 32) | 0.011 |
| intra | linear | 122 | 0.985 | 0.983 | 0.978 | 0.123 | 1.38 (75, 32) | 0.16 |
| intra | log | 122 | 0.988 | 0.987 | 0.983 | 0.106 | 1.01 (75, 32) | 0.51 |
| interfrac | linear | 122 | 0.977 | 0.973 | 0.934 | 0.258 | 2.55 (75, 32) | 0.0021 |
| interfrac | log | 122 | 0.986 | 0.985 | 0.956 | 0.212 | 1.31 (75, 32) | 0.2 |

R² for `interfrac` is the deviance pseudo-R²; its lack of fit is deviance-based and its residual tests use deviance (Shapiro–Wilk) and Pearson (Breusch–Pagan) residuals. The dispersion is estimated, so the autocorrelated frames are not treated as independent trials.

## What could not be done, and choices made

* **Wet-lab constructs (video item 6): none selected.** No Monoblock, Triblock-32, Triblock-64 or Tetrablock run is in the SubM4 dataset. Folders with those names exist from earlier submodules, but none has a readable trajectory here (see `video_runs.md` §6).
* **Git commit of the scripts:** there is none, because `analysis/` is untracked. Hashes are in `MANIFEST.md`.
* **Titus's 20 `sub4-extra` runs** have summaries only, not trajectories, so they are in none of these fits (they lack the free-only and intra-chain counts). They enter only the team-definition-A sensitivity fit in `stats/`.
* **"Bridge fraction"** is taken as the non-reactive inter / (inter + intra) free-contact fraction (`interfrac`). The crosslink bridge fraction of the reactive anchors A10–A19 has 10 design points, too few for a 15-term quadratic; its main-effects analysis is in `stats/reactive_*`.
* **Random forest** importances in `figures_logcoded/importance_comparison` are the frozen ones, not a refit: trees split on order, so a monotone recoding of the inputs gives the same partitions.
* **Video item 4:** no run has both the longest chain and the widest spacing; the choice rule and the two single-criterion runs are in `video_runs.md`.
* **Zip:** every file is under 25 MB, nothing left out.
* **Not yet included:** 30 prepared SubM4 runs that are in no sheet listed in `scripts/common.py` (e.g. the `sub4-main-*` runs of `Chr_sub4_run5.csv` and the `-nr` A10–A19 reruns). Add their sheets to `LOCAL_SHEETS` when they finish and rerun the commands below.

## Regenerate

From the repository root on DelftBlue (`run.sh` layers the pinned analysis packages on the project venv):

```bash
analysis/scripts/run_all.sh            # steps 00-09 in order
```

or step by step:

```bash
.venv/bin/python analysis/scripts/00_inventory.py
.venv/bin/python analysis/scripts/01_extract_timeseries.py --workers 8   # only runs without a .npz; or sbatch analysis/scripts/extract.sbatch
analysis/scripts/run.sh analysis/scripts/02_equilibration.py
analysis/scripts/run.sh analysis/scripts/03_metrics.py
analysis/scripts/run.sh analysis/scripts/04_stats.py
analysis/scripts/run.sh analysis/scripts/05_figures.py
analysis/scripts/run.sh analysis/scripts/06_reactive.py
analysis/scripts/run.sh analysis/scripts/07_models.py
analysis/scripts/run.sh analysis/scripts/05_figures.py --stats analysis/export/models/inter_log --out analysis/export/figures_logcoded --coding log
analysis/scripts/run.sh analysis/scripts/08_video_runs.py
analysis/scripts/run.sh analysis/scripts/09_export.py
```

New runs: add the sheet to `LOCAL_SHEETS` in `scripts/common.py`, then run everything. Nothing in `settings.yaml` changes: the burn-in and every other definition stay frozen. The video selection is re-derived by rule, so it can change when new runs arrive (e.g. a new lowest density).
