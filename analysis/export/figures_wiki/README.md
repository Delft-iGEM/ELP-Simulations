# SubM4 wiki figures

Generated 2026-10-16 by `analysis/scripts/10_wiki_figures.py` (repo HEAD `6725a2d`; `analysis/` is untracked, so the script copies in `scripts/` are the record). Data as exported on 2026-10-16 08:05 (`runs_metrics.csv` built 07:07).

**Nothing was refitted differently and no simulation or model changed.** The script refits the same three log-coded models with the export's own code (`07_models.py`, `model.py`) and asserts that every coefficient equals `models/{inter,intra,interfrac}_log/coefficients.csv`. The log coding was chosen after the frozen model's lack of fit had been seen; any claim built on it should say so. The frozen linear model stays the baseline.

**Runs in each fit:** 104 finished runs (84 design points, 20 pure-error df) for `inter_log`, `intra_log` and `interfrac_log`. Reactive figure: 30 runs, 10 design points (A10–A19). Titus's 20 `sub4-extra` runs are in none of these fits.

## Style

PNG with transparent background at 2x (1200 px wide figure = 2400 px; the older figures keep their original, smaller canvases), plus SVG. Text, axes and grid are mid-grey `#7c8089` (about 4.2:1 on both white and `#1a1a1a`); data use the team palette. White marker outlines were removed. Checked by compositing every PNG over white and over `#1a1a1a`. Poppins. `data/` holds the plotted numbers: one compact JSON per figure (`<name>.json`, all traces) plus the CSVs named below.

## Figures

| file | what | model / data |
|---|---|---|
| `coefficients` | standardised coefficients, 95 % CI | `models/inter_log/coefficients.csv` |
| `predicted_vs_heldout` | 5-fold grouped CV | `models/inter_log/run_diagnostics.csv` |
| `noise_check_anchors` | replicate SD / CV vs mean | `runs_metrics.csv` |
| `response_surface_2x2` | density x spacing at four L, I corners | refit `inter_log` |
| `importance_comparison` | LMG, RF permutation, SHAP | `models/inter_log/parameter_tests.csv`, frozen RF |
| `residual_diagnostics` | residuals, Q-Q, Cook's D | `models/inter_log/run_diagnostics.csv` |
| `qc_timeseries_sample_inter_free_12A` | 12 random runs (seed from `02_equilibration.py`), per frame and 100-frame mean | `timeseries/*.npz` |
| `effect_curves` | one-parameter-at-a-time predictions, 95 % CI of the mean, base dot, fold change vs base, raw runs | `inter_log`; `data/effect_curves_endpoints.csv` |
| `knobs_fold_change` | low and high end vs base, CI whiskers | `inter_log`, matches `effects_oat.csv`; `data/knobs_fold_change.csv` |
| `response_surface_base` | stiffness map + tissue bands (k = 1 drawn) | `inter_log`; `data/tissue_bands.csv` (k = 1 and 2) |
| `bridges_vs_loops` | A: fold change in bridges (inter) and loops (intra); B: bridge share map, 50 % line | `inter_log`, `intra_log`, `interfrac_log` (quasi-binomial logit); `data/bridges_vs_loops_panelA.csv` |
| `reactive_bridge_fraction` | observed bridge fraction vs predicted bridge share, A10–A19 | `stats/reactive_design_points.csv` |

`data/model_comparison.csv` is the two-model table for the stats dropdown.

Choices to know about:
* **Axes in `effect_curves`:** all four x-axes are log, because the model is linear in the log of every parameter. The design itself is space-filling, not geometric, so this is a modelling choice, stated in each axis label.
* **Raw dots in `effect_curves`:** runs whose other three parameters are within 0.5 coded log units of the base (a quarter of each log range). That leaves 10 (density), 20 (length), 14 (spacing) and 20 (ionic strength) runs, enough to show, but they are not exactly at base, so they scatter around the curve; the upper edge of spacing is truncated because the base is already at 0.77 coded.
* **Fold labels** are curve end / curve at base (e.g. density x0.16 at 0.01 and x1.3 at 0.06); the high/low ratios are in `knobs_fold_change` and `model_comparison`.
* **Tissue bands:** G'0 = 3 kPa at the base-case proxy (0.429), E = 3G. **The tissue E ranges (brain 0.1–1, fat 1–4, liver 0.5–2, muscle 8–17 kPa) are placeholders marked `# UNVERIFIED` in the script, with no citation. Check them before the wiki uses them.** With k = 1 the map covers proxy 0.025 to 1.23: fat and muscle fall inside, liver and brain only partly (brain 0.005–0.048 is mostly below the lowest value the map reaches, i.e. softer than any simulated construct). With k = 2 all four fall inside. Everything hinges on G'0 and k, which one rheology measurement fixes.
* **`reactive_bridge_fraction`** uses the `pred_bridge_share` column of `stats/reactive_design_points.csv`, which comes from the frozen linear-coded OLS on the inter/intra ratio (06_reactive.py), not the log-coded GLM of `bridges_vs_loops` panel B. Side check with the log GLM: Spearman 0.87, Pearson 0.98 (in the JSON).
* **Not re-exported:** `qc_timeseries_sample_rg_mean_A` (not in the list).
* The `style.py` change is a `wiki_mode()` switch; the default behaviour of every older script is unchanged. The SHA-256 of `style.py` in `export/MANIFEST.md` is now stale until `09_export.py` is rerun.

## Do the two models agree? (`data/model_comparison.csv`)

| | frozen linear | log, all four |
|---|---|---|
| R², adj. R² | 0.960, 0.954 | 0.984, 0.981 |
| grouped 5-fold CV R² | 0.882 | 0.944 |
| lack of fit F (69, 20), p | 6.61, 1e-5 | 2.51, 0.012 |
| without extra2-05 (68, 20) | 4.59, p = 2e-4 | 1.69, p = 0.0947 |
| density, 0.01 to 0.06 | x5.4 [4.2, 6.8] | x8.3 [7.2, 9.5] |
| chain length, 50 to 150 | x0.46 [0.37, 0.57] | x0.40 [0.35, 0.45] |
| spacing, 5 to 25 % | x0.62 [0.50, 0.78] | x0.63 [0.55, 0.73] |
| ionic strength, 0.05 to 0.5 M | x1.27 [1.04, 1.54] | x1.25 [1.09, 1.43] |

**Yes on sign and ranking, and for three of the four sizes; no on the size of the density effect.** Signs and the order of importance (density >> length ~ spacing > ionic strength) are the same, and the variance shares are the same (81 / 6 / 7 / 6 % for both). Length, spacing and ionic strength agree within overlapping CIs. The density effect is about 1.5 times larger in the log model and the two CIs do not overlap: the frozen model underestimates it because it forces a spurious maximum near the top of the density range. Use x8.3 only with the log-coding caveat; x5.4 is the pre-specified number.

## Checks against the wiki text

| wiki text | export | verdict |
|---|---|---|
| density x8.3 [7.2, 9.5] | 8.27 [7.22, 9.47] | matches |
| length x0.40, spacing x0.63, ionic x1.25 | 0.396, 0.634, 1.247 | match |
| variance shares 81 / 6 / 7 / 6 % | 81.3 / 6.2 / 6.8 / 5.7 % (LMG share of explained R²) | match |
| R² 0.984, adj. 0.981, CV 0.944 | 0.9838, 0.9812, 0.9435 | match |
| lack of fit F 2.5 (p = 0.012) | 2.51 (p = 0.0120) | matches |
| without extra2-05: F 1.7 (p = 0.09) | 1.69 (p = 0.0947) | matches |
| loops: spacing x0.072, length x0.121 | 0.0724, 0.1214 | match |
| bridge share base 86 % | 85.5 % | matches (rounds to 86, or say 85.5) |
| 46 to 87 % with spacing | 45.8 to 87.2 % | match |
| 81 to 93 % with length | 80.8 to 92.9 % | match |
| 48 to 89 % with density | 48.2 to 89.4 % | match |
| reactive Spearman 0.90, Pearson 0.97 | 0.903, 0.972 (against the frozen-model bridge share) | match; with the log GLM they are 0.87 and 0.98, which differ by more than rounding, so state which model the wiki means |

Nothing differs by more than rounding for the stated model choices. Two things the wiki should say: the "without extra2-05" figure is a refit on 103 runs (68, 20 df), and the reactive correlation depends on which bridge-share model is used.

## Regenerate

From the repository root (the analysis steps 00–07 must already have run; they have):

```bash
analysis/scripts/run.sh analysis/scripts/10_wiki_figures.py     # about 2 minutes, writes analysis/export/figures_wiki/
# or everything, end to end:  analysis/scripts/run_all.sh
```

## When the remaining runs arrive

* **main-100 runs** (`Chr_sub4_run5.csv`): change the model, so **regenerate everything**. Add the sheet to `LOCAL_SHEETS` in `scripts/common.py`, rerun `run_all.sh` (steps 00–07 refit, then this script). Every number in the checks table above will move, so recheck it; the lack-of-fit and the "without extra2-05" lines especially.
* **Non-reactive A10–A19 reruns** (`rerun_A10-A19_nonreactive_30.csv`, `-nr` runs): they become ordinary contact runs, so they enter the fits too (regenerate everything), and `reactive_bridge_fraction` can then be compared against measured, not only predicted, non-reactive bridge shares.
* **Titus's data**: his 20 `sub4-extra` runs have summaries only, so they affect only the team-definition-A fit, not these figures; if his trajectories arrive, they join the fits (regenerate everything). New reactive runs change only `reactive_bridge_fraction`.
* **A rheology measurement:** only `response_surface_base`, `data/tissue_bands.csv` (set `G0_KPA` and the exponent k), and the tissue dict need updating.
