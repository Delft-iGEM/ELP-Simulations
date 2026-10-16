# Manifest — SubM4 export

Generated 2026-10-16 08:05 by `analysis/scripts/09_export.py`. Repository HEAD: `6725a2d` (2026-10-15 18:50:54 +0200).

The analysis scripts are **not committed** (`analysis/` is untracked in git), so no commit identifies them; their SHA-256 hashes at export time are listed below, and `scripts/` holds the copies.

| script | SHA-256 |
|---|---|
| `00_inventory.py` | `83fba234b9c2d0b7b1d6f4ade93b10ab749297796d8bbc60b82ca0651b70385d` |
| `01_extract_timeseries.py` | `fee995eb5ee8d163cd7f68015df7c3aa15b644c4fc1fc5cada39a71de5118e79` |
| `02_equilibration.py` | `9164998e51e61e98e317caa17be3a22ba3da3355625ba56939722b93ed3fb425` |
| `03_metrics.py` | `8c24e761ba4aa430e989b48c6c2311f565870107f1bd48376edcc490442b24ef` |
| `04_stats.py` | `fe5c3b845e99c2f0a664e24dd2e09f475300afefb6c1092a0719605b245407be` |
| `05_figures.py` | `0f0e763a29771e3f423dd6df3d60e2d2094ac2fc4464f975dc0d8ccf6c5660bd` |
| `06_reactive.py` | `9b461d825a0d18e479e9d42c1f1a40c5a8faf37c0b2dfb40d3b9ff8e9883bb06` |
| `07_models.py` | `1e38c00e81e178f0a13ac34fe5b4c74960c0ff6a4283b418820481deb31730ed` |
| `08_video_runs.py` | `4590666a77e43acadfbbd5561e811ff207b12b817baaa8a6b67f64b934960afc` |
| `09_export.py` | `15cd5eaa6a5af3762afc0625e880141f1ca6ef55a5ef41771a609745743bbccd` |
| `common.py` | `09931915ed7e6b10335528e9282170d03f45bee08b225c6a5791358056456021` |
| `extract.sbatch` | `ac2d22c1d0181d60530ba1ef71171eb8ee3452ece19f1cdabb10ae019e443385` |
| `model.py` | `8df9474b17e07fdbc492226bf7e480fbf560609990562a48cb50a9a23a5f2235` |
| `run.sh` | `44dbb731428e4ef414b14d1931aa074abac2aaaf1a290233f5de65a8907181b5` |
| `run_all.sh` | `08b8311fb5547f58ac4d3346ecb3978cfba2a2a8bcf2356a600e38c0f75ea35b` |
| `style.py` | `dee0745247fd93078cc376470c53266856bfe7ff017324f115b46c2be2340fcc` |

## Files

Paths are relative to `export/`. Sizes in kB. Columns are given for CSV files.

### `MANIFEST.md` (34 kB)

This file.

### `README.md` (6 kB)

Handoff: what is here, what is new, what could not be done, how to regenerate.

### `analysis_README.md` (27 kB)

The full analysis write-up (analysis/README.md): methods, QC, results, open questions.

### `figures/coefficients.html` (13 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/coefficients.png` (112 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/importance_comparison.html` (11 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/importance_comparison.png` (120 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/noise_check_anchors.html` (12 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/noise_check_anchors.png` (177 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/predicted_vs_heldout.html` (19 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/predicted_vs_heldout.png` (163 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/qc_timeseries_sample_inter_free_12A.html` (2507 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/qc_timeseries_sample_inter_free_12A.png` (436 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/qc_timeseries_sample_rg_mean_A.html` (4298 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/qc_timeseries_sample_rg_mean_A.png` (672 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/residual_diagnostics.html` (35 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/residual_diagnostics.png` (193 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/response_surface_2x2.html` (337 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/response_surface_2x2.png` (570 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/response_surface_base.html` (97 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures/response_surface_base.png` (292 kB)

Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8.

### `figures_logcoded/coefficients.html` (13 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/coefficients.png` (118 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/importance_comparison.html` (11 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/importance_comparison.png` (124 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/noise_check_anchors.html` (12 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/noise_check_anchors.png` (177 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/predicted_vs_heldout.html` (19 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/predicted_vs_heldout.png` (164 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/residual_diagnostics.html` (35 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/residual_diagnostics.png` (194 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/response_surface_2x2.html` (337 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/response_surface_2x2.png` (578 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/response_surface_base.html` (97 kB)

The same figures for the log-coded primary model.

### `figures_logcoded/response_surface_base.png` (303 kB)

The same figures for the log-coded primary model.

### `models/comparison.csv` (2 kB)

All six new fits side by side: R², adjusted R², grouped CV, lack of fit, residual tests, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`, `cv5_R2`, `cv5_RMSE_link`

### `models/inter_linear/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_log`, `se`, `t`, `p`, `ci_low`, `ci_high`, `vif`, `fold_change_per_unit_coded`

### `models/inter_linear/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/inter_linear/importance_rf.csv` (0 kB)

Random-forest importances copied from stats/ (trees are invariant to a monotone recoding).

Columns: `parameter`, `perm_importance_train`, `perm_importance_train_sd`, `perm_importance_heldout`, `mean_abs_shap`

### `models/inter_linear/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `sd_pure_error_log`, `ss_lack_of_fit`, `df_lack_of_fit`, `sd_residual_log`, `F`, `p`, `n_design_points`, `scale`

### `models/inter_linear/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2_ols`, `cv5_grouped_RMSE_log_ols`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/inter_linear/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/inter_linear/run_diagnostics.csv` (27 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_ols`

### `models/inter_log/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_log`, `se`, `t`, `p`, `ci_low`, `ci_high`, `vif`, `fold_change_per_unit_coded`

### `models/inter_log/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/inter_log/importance_rf.csv` (0 kB)

Random-forest importances copied from stats/ (trees are invariant to a monotone recoding).

Columns: `parameter`, `perm_importance_train`, `perm_importance_train_sd`, `perm_importance_heldout`, `mean_abs_shap`

### `models/inter_log/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `sd_pure_error_log`, `ss_lack_of_fit`, `df_lack_of_fit`, `sd_residual_log`, `F`, `p`, `n_design_points`, `scale`

### `models/inter_log/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2_ols`, `cv5_grouped_RMSE_log_ols`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/inter_log/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/inter_log/run_diagnostics.csv` (27 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_ols`

### `models/interfrac_linear/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_logit`, `se`, `z`, `p`, `ci_low`, `ci_high`, `vif`, `odds_ratio_per_unit_coded`

### `models/interfrac_linear/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/interfrac_linear/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `dispersion_pure_error`, `ss_lack_of_fit`, `df_lack_of_fit`, `dispersion_model`, `F`, `p`, `n_design_points`, `scale`

### `models/interfrac_linear/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2`, `cv5_grouped_RMSE_logit`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/interfrac_linear/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/interfrac_linear/run_diagnostics.csv` (28 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `n_hits`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_logit`

### `models/interfrac_log/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_logit`, `se`, `z`, `p`, `ci_low`, `ci_high`, `vif`, `odds_ratio_per_unit_coded`

### `models/interfrac_log/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/interfrac_log/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `dispersion_pure_error`, `ss_lack_of_fit`, `df_lack_of_fit`, `dispersion_model`, `F`, `p`, `n_design_points`, `scale`

### `models/interfrac_log/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2`, `cv5_grouped_RMSE_logit`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/interfrac_log/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/interfrac_log/run_diagnostics.csv` (28 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `n_hits`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_logit`

### `models/intra_linear/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_log`, `se`, `t`, `p`, `ci_low`, `ci_high`, `vif`, `fold_change_per_unit_coded`

### `models/intra_linear/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/intra_linear/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `sd_pure_error_log`, `ss_lack_of_fit`, `df_lack_of_fit`, `sd_residual_log`, `F`, `p`, `n_design_points`, `scale`

### `models/intra_linear/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2_ols`, `cv5_grouped_RMSE_log_ols`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/intra_linear/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/intra_linear/run_diagnostics.csv` (27 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_ols`

### `models/intra_log/coefficients.csv` (2 kB)

Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit.

Columns: `term`, `coef_log`, `se`, `t`, `p`, `ci_low`, `ci_high`, `vif`, `fold_change_per_unit_coded`

### `models/intra_log/effects_oat.csv` (1 kB)

One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio.

Columns: `parameter`, `base_value`, `low_value`, `high_value`, `base_pred`, `base_ci_low`, `base_ci_high`, `low_pred`, `low_ci_low`, `low_ci_high`, `high_pred`, `high_ci_low`, `high_ci_high`, `ratio_kind`, `ratio_high_over_low`, `ratio_ci_low`, `ratio_ci_high`

### `models/intra_log/lack_of_fit.csv` (0 kB)

Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM).

Columns: `ss_pure_error`, `df_pure_error`, `sd_pure_error_log`, `ss_lack_of_fit`, `df_lack_of_fit`, `sd_residual_log`, `F`, `p`, `n_design_points`, `scale`

### `models/intra_log/model_fit.csv` (1 kB)

R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance.

Columns: `response`, `coding`, `coding_note`, `model`, `n_runs`, `n_design_points`, `n_params`, `R2`, `adj_R2`, `R2_kind`, `cv5_grouped_R2_ols`, `cv5_grouped_RMSE_log_ols`, `lof_F`, `lof_df_num`, `lof_df_den`, `lof_p`, `dispersion`, `shapiro_W`, `shapiro_p`, `shapiro_on`, `breusch_pagan_p`, `breusch_pagan_on`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `models/intra_log/parameter_tests.csv` (1 kB)

Per-parameter F-test (all terms containing the parameter dropped) and LMG shares.

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `models/intra_log/run_diagnostics.csv` (27 kB)

Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_ols`

### `models/lmg_shares.csv` (3 kB)

Per-parameter F-tests and LMG shares for all six fits.

Columns: `response`, `coding`, `parameter`, `F`, `p`, `lmg_share_R2`, `lmg_share_of_explained`

### `per_lysine.csv` (429 kB)

Per lysine per run: chain, residue (0-based in chain), K index, pinned, mean free inter-chain partners at 12 Å.

Columns: `slug`, `chain`, `residue_in_chain`, `pinned`, `mean_inter_partners_12A`, `k_index`

### `qc/equilibration_per_run.csv` (35 kB)

Chodera burn-in per run and observable (frame index).

Columns: `slug`, `observable`, `t0_frame`, `t0_ns`, `g_after_t0`, `n_eff`, `g_full`

### `reactive_anchors.csv` (10 kB)

Titus's reactive anchors A10-A19: team-A contacts and crosslinks formed / loops / bridges. Never in the contact fits.

Columns: `Name`, `slug`, `batch`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_repeats`, `spacing_pct`, `ionic_strength`, `lead_repeats`, `trail_repeats`, `Sequence`, `concentration_planned`, `concentration_rel_change`, `meta_n_residues`, `meta_box_x_nm`, `meta_box_z_nm`, `meta_box_headroom_nm`, `meta_walltime_used_fraction`, `crosslink_distance`, `crosslink_prob`, `inter_teamA_12A`, `intra_teamA_12A`, `inter_teamA_pairs_ever`, `xl_formed`, `xl_loops`, `xl_bridges`, `xl_max`, `xl_conversion`, `xl_bridge_frac`

### `run_inventory.csv` (70 kB)

Every run-sheet row matched to its output, with status (finished / summary_only / excluded / missing) and the reason.

Columns: `Name`, `slug`, `sheet`, `batch`, `anchor_id`, `replicate`, `chain_length`, `n_lys`, `spacing_repeats`, `spacing_pct`, `lead_repeats`, `trail_repeats`, `Sequence`, `concentration_planned`, `ionic_strength`, `steps`, `mode`, `preattached_fraction`, `walltime`, `crosslink_distance`, `crosslink_prob`, `meta_status`, `meta_n_residues`, `meta_n_frames_saved`, `meta_run_legs`, `meta_walltime_used_fraction`, `meta_box_x_nm`, `meta_box_z_nm`, `meta_z_max_nm`, `meta_box_headroom_nm`, `meta_box_tall_enough`, `meta_spacing_nm`, `meta_concentration_chains_per_nm2`, `meta_lattice_matches_box`, `meta_margin_nm`, `meta_ionic_strength_M`, `meta_lysines_total`, `meta_surface_bonds`, `meta_surface_bonds_run`, `meta_surface_chains_with_1`, `meta_surface_chains_with_2`, `meta_surface_chains_with_3plus`, `meta_crosslinking`, `meta_surface_seed`, `meta_total_time_ns`, `concentration`, `concentration_rel_change`, `status`, `reason`

### `runs_metrics.csv` (103 kB)

One row per run: design parameters, all per-run metrics with block SEs, QC flags. `data_source` = trajectory (104) or teammate_summary (20, team-definition-A columns only).

Columns: `Name`, `slug`, `data_source`, `batch`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_repeats`, `spacing_pct`, `ionic_strength`, `lead_repeats`, `trail_repeats`, `Sequence`, `concentration_planned`, `concentration_rel_change`, `meta_n_residues`, `meta_box_x_nm`, `meta_box_z_nm`, `meta_box_headroom_nm`, `meta_walltime_used_fraction`, `production_frames`, `discard_frames`, `inter_free_12A`, `inter_free_12A_se`, `inter_free_12A_per_chain`, `inter_free_12A_frac_zero_frames`, `inter_free_12A_drift`, `inter_free_12A_drift_z`, `inter_free_12A_from200ns`, `inter_free_12A_from200ns_se`, `inter_free_8A`, `inter_free_8A_se`, `inter_free_16A`, `inter_free_16A_se`, `inter_teamA_12A`, `inter_teamA_12A_se`, `intra_free_12A`, `intra_free_12A_se`, `bridge_loop_ratio_12A`, `largest_cluster_frac`, `largest_cluster_frac_se`, `inter_free_pairs_ever`, `inter_free_pairs_total`, `h95_nm`, `h95_nm_se`, `h95_nm_sd_time`, `h50_nm`, `zmean_nm`, `h95_drift_z`, `rg_nm`, `rg_nm_se`, `rg_nm_sd`, `rg_nm_sd_between_chains`, `ree_nm`, `ree_nm_se`, `ree_nm_sd`, `ree_nm_sd_between_chains`, `rg_drift_z`, `design_point`, `t0_detected_h95_A`, `t0_detected_inter_free_12A`, `t0_detected_rg_mean_A`, `qc_drift_primary`, `qc_drift_rg`, `qc_drift_h95`, `qc_burnin_exceeds_discard`, `qc_plan_box_changed`, `qc_box_headroom_lt_4nm`, `qc_any`, `excluded`, `exclusion_reason`, `intra_teamA_12A`, `inter_teamA_pairs_ever`, `inter_teamA_pairs_total`, `checks_not_possible`

### `scripts/00_inventory.py` (7 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/01_extract_timeseries.py` (9 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/02_equilibration.py` (3 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/03_metrics.py` (10 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/04_stats.py` (15 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/05_figures.py` (14 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/06_reactive.py` (5 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/07_models.py` (13 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/08_video_runs.py` (14 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/09_export.py` (18 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/common.py` (7 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/extract.sbatch` (1 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/model.py` (4 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/run.sh` (0 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/run_all.sh` (1 kB)

Copy of the analysis scripts that produced everything here.

### `scripts/style.py` (4 kB)

Copy of the analysis scripts that produced everything here.

### `settings.yaml` (2 kB)

Frozen definitions (response, pair definition, cutoff, burn-in, QC thresholds, coding ranges, seed).

### `stats/coefficients.csv` (2 kB)

Frozen model: OLS coefficients on log(inter_free_12A), SE, t, p, 95 % CI, VIF.

Columns: `term`, `coef_log`, `se`, `t`, `p`, `ci_low`, `ci_high`, `vif`, `fold_change_per_unit_coded`

### `stats/glm_gamma.csv` (2 kB)

Gamma GLM (log link) coefficients next to OLS on log y.

Columns: `term`, `coef_log`, `se`, `p`, `ci_low`, `ci_high`, `ols_coef_log`

### `stats/importance_rf.csv` (0 kB)

Random-forest permutation importance and mean |SHAP| (frozen analysis).

Columns: `parameter`, `perm_importance_train`, `perm_importance_train_sd`, `perm_importance_heldout`, `mean_abs_shap`

### `stats/lack_of_fit.csv` (0 kB)

Frozen model: lack-of-fit F-test against the anchor pure error.

Columns: `ss_pure_error`, `df_pure_error`, `sd_pure_error_log`, `ss_lack_of_fit`, `df_lack_of_fit`, `sd_residual_log`, `F`, `p`, `n_design_points`

### `stats/model_fit.csv` (0 kB)

Frozen model: R², grouped CV (OLS and random forest), residual tests, Cook's distance.

Columns: `n_runs`, `n_params`, `R2`, `adj_R2`, `cv5_grouped_R2_ols`, `cv5_grouped_RMSE_log_ols`, `cv5_grouped_R2_rf`, `cv5_grouped_RMSE_log_rf`, `shapiro_W`, `shapiro_p`, `breusch_pagan_p`, `max_vif`, `n_cooks_gt_4_over_n`, `max_cooks`

### `stats/parameter_tests.csv` (1 kB)

Frozen model: per-parameter F-tests and LMG shares of R².

Columns: `parameter`, `terms_dropped`, `F`, `df_num`, `df_den`, `p`, `unique_dR2`, `lmg_share_R2`, `lmg_share_of_explained`

### `stats/reactive_anchors.csv` (11 kB)

Reactive anchors vs the non-reactive team-A model at the same design points.

Columns: `Name`, `slug`, `batch`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_repeats`, `spacing_pct`, `ionic_strength`, `lead_repeats`, `trail_repeats`, `Sequence`, `concentration_planned`, `concentration_rel_change`, `meta_n_residues`, `meta_box_x_nm`, `meta_box_z_nm`, `meta_box_headroom_nm`, `meta_walltime_used_fraction`, `crosslink_distance`, `crosslink_prob`, `inter_teamA_12A`, `intra_teamA_12A`, `inter_teamA_pairs_ever`, `xl_formed`, `xl_loops`, `xl_bridges`, `xl_max`, `xl_conversion`, `xl_bridge_frac`, `pred_nonreactive`, `ratio_obs_over_pred`

### `stats/reactive_design_points.csv` (2 kB)

Reactive anchors per design point, with non-reactive model predictions.

Columns: `anchor_id`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `bridge_frac_mean`, `bridge_frac_sd`, `formed_mean`, `pred_inter_free_12A`, `pred_bridge_share`, `glm_bridge_frac`

### `stats/reactive_glm.csv` (1 kB)

Reactive anchors: quasi-binomial GLM of the crosslink bridge fraction (main effects).

Columns: `term`, `coef_logit`, `se`, `p`, `ci_low`, `ci_high`, `odds_ratio_per_unit_coded`, `vif`

### `stats/reactive_summary.txt` (3 kB)

Reactive-anchor analysis in plain text.

### `stats/run_diagnostics.csv` (29 kB)

Frozen model, per run: fitted, residual, leverage, Cook's D, held-out predictions.

Columns: `slug`, `design_point`, `anchor_id`, `replicate`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `logy`, `fitted`, `resid`, `resid_studentized`, `leverage`, `cooks_d`, `qc_cooks_gt_4_over_n`, `cv_pred_ols`, `cv_pred_rf`

### `stats/secondary.csv` (1 kB)

The frozen model on secondary metrics (R² and LMG shares).

Columns: `metric`, `scale`, `R2`, `lmg_C`, `lmg_L`, `lmg_S`, `lmg_ion`

### `stats/sensitivity.csv` (2 kB)

The frozen model on the sensitivity responses (8 / 16 Å, team definition A with and without Titus's runs, 200 ns burn-in).

Columns: `response`, `runs`, `n`, `R2`, `adj_R2`, `lmg_C`, `p_C`, `lmg_L`, `p_L`, `lmg_S`, `p_S`, `lmg_ion`, `p_ion`, `coef_C`, `coef_L`, `coef_S`, `coef_ion`

### `stats/shap_values.npy` (3 kB)

SHAP values of the random forest (NumPy array, runs x parameters).

### `stats/summary.txt` (10 kB)

Frozen first-pass statistics in plain text (n = 104).

### `stats/teammate_check.csv` (4 kB)

Titus's 20 runs predicted by the team-A model fitted on our runs alone (held out).

Columns: `slug`, `concentration`, `chain_length`, `spacing_pct`, `ionic_strength`, `y`, `pred_log`, `pi_low_log`, `pi_high_log`, `resid_log`, `outside_95pi`, `pred`

### `timeseries_primary.csv.gz` (1679 kB)

Per-frame inter / intra free contacts at 12 Å for every trajectory run, all frames including burn-in.

Columns: `slug`, `frame`, `inter_free_12A`, `intra_free_12A`

### `video_runs/sub4-anchor-A02-L64-S7-C0p01-I0p19-rep1_contacts.csv` (51 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-anchor-A02-L64-S7-C0p01-I0p19-rep1_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-anchor-A02-L64-S7-C0p01-I0p19-rep1_surface_state.json` (5 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-anchor-A02-L64-S7-C0p01-I0p19-rep1_top.pdb` (406 kB)

Copy of the run's topology (Å).

### `video_runs/sub4-anchor-BASE-CASE-L72-uniform-C0p04-I0p20-rep1_contacts.csv` (85 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-anchor-BASE-CASE-L72-uniform-C0p04-I0p20-rep1_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-anchor-BASE-CASE-L72-uniform-C0p04-I0p20-rep1_surface_state.json` (4 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-anchor-BASE-CASE-L72-uniform-C0p04-I0p20-rep1_top.pdb` (456 kB)

Copy of the run's topology (Å).

### `video_runs/sub4-extra2-L58-S4-C0p01-I0p43-05_contacts.csv` (122 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-extra2-L58-S4-C0p01-I0p43-05_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-extra2-L58-S4-C0p01-I0p43-05_surface_state.json` (5 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-extra2-L58-S4-C0p01-I0p43-05_top.pdb` (368 kB)

Copy of the run's topology (Å).

### `video_runs/sub4-extra3-L149-S37-C0p02-I0p24-09_contacts.csv` (27 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-extra3-L149-S37-C0p02-I0p24-09_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-extra3-L149-S37-C0p02-I0p24-09_surface_state.json` (4 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-extra3-L149-S37-C0p02-I0p24-09_top.pdb` (943 kB)

Copy of the run's topology (Å).

### `video_runs/sub4-extra4-L138-S22-C0p01-I0p08-16_contacts.csv` (21 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-extra4-L138-S22-C0p01-I0p08-16_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-extra4-L138-S22-C0p01-I0p08-16_surface_state.json` (5 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-extra4-L138-S22-C0p01-I0p08-16_top.pdb` (874 kB)

Copy of the run's topology (Å).

### `video_runs/sub4-extra4-L76-S4-C0p06-I0p35-28_contacts.csv` (355 kB)

Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md.

Columns: `frame`, `type`, `chain_A`, `resid_A`, `chain_B`, `resid_B`, `index_A`, `index_B`

### `video_runs/sub4-extra4-L76-S4-C0p06-I0p35-28_pins.csv` (0 kB)

Surface-pinned lysines: chain, resid, K index, atom index.

Columns: `chain`, `resid`, `k_index`, `index`

### `video_runs/sub4-extra4-L76-S4-C0p06-I0p35-28_surface_state.json` (4 kB)

Copy of the run's surface state (pinned lysines, tether points).

### `video_runs/sub4-extra4-L76-S4-C0p06-I0p35-28_top.pdb` (482 kB)

Copy of the run's topology (Å).

### `video_runs.md` (11 kB)

The runs chosen for the videos: reasons, absolute paths, production start frame, frame count.

### `analysis_export.zip`

Everything above in one archive (files over 25 MB left out; see README.md).
