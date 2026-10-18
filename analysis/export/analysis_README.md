# SubM4 control sweep — processed dataset and first-pass analysis

Handoff for the statistics/interpretation round. Written 2026-10-15. Everything here is
derived; the raw simulation output under `simulations/*/runtime` (→ `/scratch`) was only read.

**Status in one paragraph.** 104 of the 260 planned runs exist and were analysed (all
complete, all pass QC). The primary response (inter-chain contacts between free lysines,
12 Å) is fitted well by the quadratic model (R² 0.96, grouped-CV R² 0.88), and the
conclusions do not depend on cutoff, pair definition or burn-in. Grafting density dominates
(≈ 80 % of explained variance), followed by lysine spacing and chain length, which both reduce
the proxy; ionic strength has a small positive effect. **The lack-of-fit test is significant**,
so the quadratic is not the final model: it misses the low-density corner and puts a spurious
maximum near the top of the density range. Treat every number below as preliminary.

**Update 2026-10-17 — non-reactive A10–A15 reruns.** The `-nr` reruns of Titus's anchors have
finished for **A10–A15 (18 runs)**; they pass the same QC and enter the primary free-only
analysis, so it now has **122 runs, 90 design points, 16 anchors** (pure-error df 20 → 32).
A16–A19 (12 runs) are still running or queued and stay `incomplete` in the inventory; rerun
`analysis/scripts/run_all.sh` once they finish. The new §7d has the numbers. Everything in §5–§7
below was written for 104 runs; the figures, `stats/` and `export/` were regenerated for 122,
so the headline numbers moved slightly (R² 0.960 → 0.959, grouped-CV R² 0.882 → 0.890; lack of
fit still significant, F 6.6 → 4.9). Ranking and signs of the parameters are unchanged.

**Update 2026-10-16 — Titus's runs.** His folders arrived with the merge of `origin/main`
(`6725a2d`). They hold his notebook's per-run summaries, not trajectories, so his 20
non-reactive `sub4-extra` runs enter only the **team-definition-A** fit (124 runs; §7a). The
primary free-only analysis stays at 104 runs, and every number in §5–§7 is unchanged.
Fitted on our runs alone, the model predicts his runs held out with no bias (held-out R² 0.83).
His two misses are exactly the two lack-of-fit spots below. A10–A19 (crosslinking on) stay out
of the contact fits but are analysed on their own (§7b). There, the crosslinked bridge fraction
follows the non-reactive *bridge share* (Spearman 0.90), not the primary contact proxy (0.32).

---

## 1. What was run

| | |
|---|---|
| Model | CALVADOS2, one bead per residue, implicit solvent, OpenMM 8.2 (CUDA, DelftBlue A100) |
| System | 16 chains (VPGIG)ₙ with 4 VPGKG, `preattached` mode, 293.15 K, pH 7.5 (not used by the model) |
| Length | 40 000 000 steps × 0.01 ps = 400 ns, 5714 frames saved every 70 ps |
| Box | lattice 4 × 4, spacing = 1/√C, no margin; height 50–90 nm (0.3 × contour length, ≥ 50 nm) |
| Crosslinking | off (no K–K bonds form); the proxy counts *contacts*, i.e. crosslinking opportunity |

Run sheets submitted from this checkout: `Chr_sub4_RUN_1.csv` (anchors BASE + A01–A09 =
`submodule4_anchor_replicates_part1_30`), `CHR_sub4_RUN_2.csv` (extra2 + extra3),
`CHR_sub4_RUN_3.csv` (extra4), `CHR_sub4_RUN_4.csv` (extra5 01–14). Titus's sheets (commit
`0d88132`, merged in `6725a2d`): `submodule4_extra_20.csv`,
`submodule4_anchor_replicates_part2_30.csv`.

## 2. Inventory — `run_inventory.csv` (one row per sheet row, 154 rows)

| batch | finished & analysed | summary only (team def. A) | excluded |
|---|---|---|---|
| anchor BASE, A01–A09 (×3) | 30 | | |
| anchor A10–A19 (×3, Titus), reactive | | | **30** |
| anchor A10–A15 `-nr` rerun (×3) | **18** | | |
| anchor A16–A19 `-nr` rerun (×3) | not finished (12) | | |
| extra (Titus) | | 20 | |
| extra2 / extra3 / extra4 / extra5 | 20 / 10 / 30 / 14 | | |
| **total** | **122** | 20 | 30 |

Not in any sheet available here: the 100 `sub4-main-*` rows of `submodule4_full_design_160.csv`
and `sub4-extra5` rows 15–20. The superseded `batch1/batch2/milliQ/PBS` sheets have no runs on
disk.

**Excluded — anchors A10–A19:** run with reactive K–K crosslinking on (0.8 nm, prob 1.0, 2:30
walltime), against the protocol. Formed bonds hold lysine pairs together, so the contact count
is inflated (Titus's summary, same 12 Å cutoff and team pair definition: ~8 inter-chain
contacts/frame at C = 0.014, against ≤ 0.18 in our non-reactive runs below C = 0.018). The non-reactive rerun sheet
`rerun_A10-A19_nonreactive_30.csv` has been submitted (A10–A15 done, A16–A19 pending; see §7d) (crosslink columns blank, names `...-nr-repN`, walltime
1:45 for L ≥ 130 because L150 used 94 % of 1:30).

**Summary only — `sub4-extra` (20 runs):** non-reactive, same template as ours (`prepare.py`
differs only in run parameters and seeds). The trajectories are on `/scratch/tnartey` (mode 700).
What we have is `simulations/<name>/analysis/crosslink_contacts_K.txt` and `load_summary.txt`,
written by `analyze.ipynb` (cell 20) at 12 Å from frame 1428. That cell drops only pairs
with both lysines pinned (team definition A). Recomputed on our own runs, it reproduces
`inter_teamA_12A` exactly (checked on extra2-05 and a05-rep2). The free-only primary response cannot be
recovered from it. Checks possible from the summaries: 5714 frames, cutoff and window, residues
per chain, concentration within 0.02 % of planned (from the box), observed headroom ≥ 7.4 nm.
Not possible: equilibration, drift, block SE.

With the 122 runs: 90 distinct design points, 16 of them in triplicate (pure-error df = 32).
The team-A fit adds 20 more points (110 distinct; pure error unchanged).

## 3. Frozen definitions — `settings.yaml`

Frozen before any response was summarised by parameter; confirmed by Chris on 2026-10-15.

| item | value | source / reason |
|---|---|---|
| Primary response | mean inter-chain K–K contacts per frame (`inter_free_12A`), also ÷16 per chain | stiffness proxy |
| Pair definition | **free lysines only**: neither lysine surface-pinned | a pinned lysine's amine is used by the surface bond |
| Cutoff | **12 Å** CA–CA, minimum image in x, y, z (`box=`) | team convention (lysine-arrangement, RGD submodules) |
| Sensitivity | 8 Å, 16 Å; team definition A (drop only both-pinned pairs); 200 ns burn-in | listed before fitting |
| Intra-chain exclusion | pairs < 4 residues apart (never binds; min K–K gap is 20 residues) | notebook convention |
| Burn-in | **1428 frames = 100 ns**; production 300 ns (4286 frames) | §5 |
| Standard error | block average, 10 blocks of 30 ns | ≫ correlation times |
| Drift flag | \|2nd − 1st half of production\| / combined block SE > 3 | protocol |
| Influence flag | Cook's D > 4/n | protocol |
| Random seed | 20261015 (RF, permutation, CV) | |

Pipeline facts read from the code, not assumed:

* **Seeds are independent.** `random_number_seed: null` → OpenMM picks a fresh integrator seed per
  run; each run also has its own `surface_seed` (all 104 unique), so replicates pin different
  lysines and start from different structures (`top.pdb` hashes differ).
* **`preattached`, fraction 0.3:** at t = 0, round(0.3 × 64) = 19 lysines (≥ 1 per chain; mostly
  1, a few 2) are tethered (k = 2000 kJ mol⁻¹ nm⁻²) to fixed points 0.6 nm above the wall
  (z = 2.5 nm). `dynamic` is off, so no lysine binds during the run (verified in
  `surface_state.json` for every run).
* **Topology:** one `CA` bead per residue; lysine = resname `LYS`; chains = segments A–P.
  DCD/PDB in **Å**; CALVADOS settings in nm. **z is the surface normal** (wall at 1.9 nm).
  Chains are written whole (max bond 4.6 Å), so Rg and end-to-end distance need no unwrapping.
* **Box planning:** for every run the lattice tiles the box with no margin and the
  concentration is the planned one to within 0.02 % (spacing rounded to 3 decimals). Box
  headroom ≥ 8.8 nm (> 4 nm cutoff) everywhere.

## 4. Per-run outputs

* `runs_metrics.csv` — one row per run: name, `anchor_id`, `replicate`, `design_point` (CV
  group), the four parameters (concentration as simulated, chain length, spacing in repeats and
  %, ionic strength), lead/trail background repeats, then the metrics with block SEs:
  * `inter_free_12A` (+`_se`, `_per_chain`, `_frac_zero_frames`, `_drift`, `_drift_z`, `_from200ns`)
  * sensitivity: `inter_free_8A`, `inter_free_16A`, `inter_teamA_12A`
  * loops: `intra_free_12A`; `bridge_loop_ratio_12A` = inter / intra
  * `largest_cluster_frac` — largest connected group of chains in the contact graph, ÷ 16
  * `inter_free_pairs_ever` / `_total` — distinct free inter-chain pairs ever within 12 Å
  * brush height: `h95_nm` (95 % of chain mass below it, measured from the wall), `h95_nm_sd_time`,
    `h50_nm`, `zmean_nm`
  * `rg_nm`, `ree_nm` — mean, SE, SD over chains × frames, SD between chains
  * QC columns `qc_*`, `t0_detected_*`, and `excluded` / `exclusion_reason` (all False / empty)
  * `data_source`: `trajectory` (our 104) or `teammate_summary` (Titus's 20). Teammate rows
    have only the design columns, `inter_teamA_12A`, `intra_teamA_12A`,
    `inter_teamA_pairs_ever/_total`, the box QC and `checks_not_possible`. Every other metric is
    empty, so `model.load()` keeps them out of all fits except the team-A one
* `reactive_anchors.csv` — Titus's A10–A19 (crosslinking on): team-A contacts per run plus
  crosslinks formed / loops / bridges (`xl_*`) from his `figures/sub4/crosslinks-submodule4-anchor.csv`.
  Never in the contact fits; analysed separately in §7b
* `per_lysine.csv` — per lysine per run: chain, residue, K1–K4 index, pinned, mean number of
  free inter-chain partners at 12 Å (NaN for pinned lysines). For the structure figures.
* `timeseries_primary.csv.gz` — per-frame inter/intra free contacts at 12 Å, every run, all
  frames (including burn-in), for re-checking equilibration.
* `timeseries/<run>.npz` — full per-frame extraction (all cutoffs, both pair definitions,
  cluster sizes, per-lysine degrees, 100-frame block pair hits, Rg/Ree per chain, height
  quantiles, block density profiles). 83 MB; these regenerate everything above without
  touching the trajectories.

## 5. Quality control

* **Completion:** 104/104 completed, 5714/5714 frames, one SLURM leg each (none hit the walltime;
  maximum use 94 %, L150).
* **Equilibration** (`qc/equilibration_per_run.csv`, Chodera max-N_eff burn-in):

  | observable | median | 90th pct | 95th pct | max |
  |---|---|---|---|---|
  | primary response | 0 ns | 8.5 ns | 22 ns | 57 ns |
  | 95 % height | 1.4 ns | 51 ns | 96 ns | 148 ns |
  | chain-mean Rg | 7 ns | 94 ns | 139 ns | 199 ns |

  Rg visibly relaxes in 20–40 ns (`figures/qc_timeseries_sample_*`); the later detections sit on
  slow fluctuations: mean Rg over 100–250 ns vs 250–400 ns differs by −3.8 % to +3.7 % (mean
  0.0 %), with no direction. 100 ns was therefore kept; the 14 runs with a later Rg/height
  detection are flagged (`qc_burnin_exceeds_discard`), not excluded. Primary response with a
  200 ns burn-in correlates r = 0.9992 with the 100 ns one.
* **Drift:** 4 runs flagged on the primary response (|z| 3.1–4.6) — a7-rep3, extra3-07,
  extra4-07, extra4-28. The shifts are 3–8 % of the mean; they are flagged only because
  within-run SEs are tiny (median 3 %). 3 runs flagged on height, 1 on Rg. All kept.
* **Influence** (`stats/run_diagnostics.csv`): 9 runs have Cook's D > 4/n. One stands out:
  `sub4-extra2-l58-s4-c0p01-i0p43-05` (D = 0.95, studentised residual −6.2: observed 0.040, fitted
  0.076). It is at the low-density, narrow-spacing corner, has no QC problem, and is kept.
  The A02 replicates (C = 0.010) are the next most influential.
* Nothing was excluded on QC grounds.

## 6. Noise check (step 4) — `figures/noise_check_anchors`

Replicate CV across the 16 anchors (10 original + A10–A15 reruns): 0.9–21 % (pooled SD of
log y = 0.078, i.e. ~8 %). A10 (C = 0.014, the lowest-density rerun) is the noisy one at 21 %; the
other five reruns are at 2–6 %. The replicate SD rises with the mean (Spearman +0.35, n.s.) while the CV does not (−0.34, n.s.): the noise is
multiplicative, so **the log transform is appropriate**. No run has a zero mean (range
0.040–1.36), so log y needs no offset.

Within-run block SE understates run-to-run noise: replicate SD is 0.3–4.4× the block SE (A05
highest). Different replicates pin different random lysines, which is real design-relevant
variability, so use the replicate pure error, never the per-run SE, as the noise floor.

## 7. First-pass statistics (step 5) — `stats/`

Model: OLS on log(`inter_free_12A`), four parameters coded to [−1, +1] over the protocol
ranges, full quadratic (14 terms + intercept), n = 122 (104 before the A10–A15 reruns).

* **Fit:** R² 0.960, adjusted R² 0.954; 5-fold grouped CV (replicates kept together) R² 0.882,
  RMSE 0.247 in log units (~28 %). Random forest: CV R² 0.871. The quadratic generalises as well
  as the forest. Maximum VIF 1.63 (the LHS design is close to orthogonal).
* **Per-parameter F-tests** (dropping all 5 terms containing the parameter) and LMG share of R²:

  | parameter | F (5, 89) | p | share of explained variance | random forest held-out permutation |
  |---|---|---|---|---|
  | grafting density C | 314 | 7e−55 | 82 % | dominant |
  | lysine spacing S | 29 | 2e−17 | 7 % | 2nd |
  | chain length L | 17 | 1e−11 | 6 % | 3rd |
  | ionic strength I | 4.8 | 6e−4 | 6 % | small |

  The RF and SHAP give the same ranking (C ≫ S > L > I). The RF puts less weight on I than the
  LMG share does.
* **Effect sizes**, one at a time from the base case (model predictions, 95 % CI; base case
  predicted 0.44 [0.40, 0.49], observed BASE anchor 0.41):

  | change, others at base | proxy low end | proxy high end | ratio |
  |---|---|---|---|
  | C 0.01 → 0.06 | 0.080 [0.068, 0.094] | 0.43 [0.36, 0.51] | ×5.4 |
  | L 50 → 150 | 0.56 [0.47, 0.66] | 0.26 [0.22, 0.30] | ×0.46 |
  | S 5 → 25 % | 0.61 [0.52, 0.71] | 0.38 [0.33, 0.45] | ×0.62 |
  | I 0.05 → 0.50 M | 0.37 [0.32, 0.43] | 0.46 [0.41, 0.52] | ×1.27 |

* **Interactions with CIs excluding 0:** C×S (−0.47), C×L (−0.22), L×S (−0.24): spacing and
  length cost more at high density. Nothing involving ionic strength interacts (C×I, L×I, S×I
  CIs all include 0); I² is negative (−0.13).
* **Lack of fit: significant.** F(69, 20) = 6.6, p = 1e−5; residual SD 0.156 vs pure-error SD
  0.067 (log units). The quadratic leaves systematic structure:
  1. *Low-density corner.* Residuals curve at low fitted values; A02 (C = 0.010) is
     over-predicted ~2.4× when held out.
  2. *Spurious optimum.* C² = −0.74 puts a maximum at C ≈ 0.052. The binned raw means
     (0.105, 0.25, 0.39, 0.52, 0.62, 0.65 across C bins) and the forest's partial dependence rise
     monotonically and flatten. Nothing turns over. The bend-back at high C and wide spacing in the
     surface figures is an artefact.
  3. *Exploratory only, not part of the frozen analysis:* coding C on a log scale raises R² to
     0.979 and halves the lack of fit (F 6.6 → 3.4, residual SD 0.114). It is still significant.
* **Residuals:** Shapiro–Wilk W = 0.946, p = 4e−4 (heavy tails, driven by the extra2-05 outlier);
  Breusch–Pagan p < 1e−4 (heteroscedastic even on the log scale).
* **Gamma GLM, log link** (`glm_gamma.csv`), used instead of the negative binomial in the brief:
  the per-run response is a time average of autocorrelated frames, not a count of independent
  events. All coefficients agree with OLS-on-log to within ~0.04 and the same terms are
  significant.
* **Sensitivity** (`sensitivity.csv`): at 8 Å / 16 Å, team definition A, and 200 ns burn-in,
  R² is 0.94–0.96, all four parameters have p < 0.02, LMG shares change by ≤ 2 points, and every
  linear coefficient keeps its sign. The conclusions do not depend on these choices.
* **Secondary metrics** (`secondary.csv`, same model, LMG shares): Rg, end-to-end distance and brush
  height are governed by chain length (76–91 %), as a sanity check expects. Loops
  (`intra_free_12A`) are governed by spacing (54 %) and length (40 %), barely by density. The
  bridge-to-loop ratio splits density 51 %, spacing 25 %, length 20 %. Largest cluster
  fraction: density 73 %.

**Plain-language first pass (preliminary):** within the sampled ranges, the crosslinking
opportunity between chains is set mainly by how densely the chains are grafted. It rises about
5-fold from 0.01 to 0.06 chains/nm² and levels off at the top of that range. The other levers
give finer tuning. Spacing the four lysines more widely and making the chain longer both lower
inter-chain contacts (each by about a factor of 2 across its range), and both effects grow with
density. Higher salt raises contacts slightly (~27 % from 0.05 to 0.5 M), consistent with
screening the K⁺–K⁺ repulsion, and it does not interact with anything. Lysine spacing mainly
shifts lysines from bridges to loops, which matters because loops carry no load.

### 7a. Titus's runs (added 2026-10-16) — `stats/teammate_check.csv`, `stats/sensitivity.csv`

* **Held-out agreement.** The team-A model fitted on our 104 runs predicts his 20 runs with
  R² 0.83, RMSE 0.41 log units (our own grouped CV: 0.25), mean residual −0.07 (p = 0.44, no
  offset between the two accounts). 18/20 are inside the 95 % prediction interval.
* **The two misses are the known lack-of-fit spots:**
  `extra-04` (C 0.011, spacing 4.9 %, just below the 5 % range): observed 0.009, predicted 0.048.
  This is the low-density, narrow-spacing corner again (as extra2-05, A02), now more extreme.
  It is the lowest response of all 124 runs.
  `extra-15` (C 0.060, spacing 24 %): observed 0.63, predicted 0.37. This is where the quadratic
  bends back, and the run shows no turnover.
* **Team-A fit on 124 runs:** R² 0.938 (104 runs: 0.962). Ranking and signs unchanged
  (share of explained variance C 82 %, L 7 %, S 7 %, I 4 %; all p ≤ 0.001). Lack of fit gets worse (F 7.1 → 12.1,
  residual SD 0.158 → 0.209), because the new points probe the weak corners.
* *Exploratory, not frozen:* with log-coded C, held-out R² rises to 0.88 and extra-15 is fitted
  (residual +0.55 → +0.16), but extra-04 is not (−1.63 → −1.48). On 124 runs, log C gives
  R² 0.963 and LOF F 7.2. The low-density, narrow-spacing drop is steeper than either form allows.
* **Reactive anchors A10–A19:** 7.5–20 team-A contacts per frame, 26–95× the non-reactive
  model at the same design points. Crosslinks lock pairs together, so these runs measure
  something else. They are not usable as replicates (§7b analyses them on their own).

### 7b. Reactive anchors A10–A19 on their own (added 2026-10-16) — `scripts/06_reactive.py`, `stats/reactive_*`

30 runs, 10 design points × 3, crosslinking on (0.8 nm, p = 1, valence 1). The data come from
Titus's table. His crosslink event files are on his scratch, so the counts could not be re-derived.

* **Conversion is saturated:** 20–22 of the 22 possible bonds form in every run (45 free lysines,
  19 pinned; checked in all 30 summaries). The number of crosslinks carries no design information.
  What varies is **where** they form: bridge fraction (inter-chain / all) 0.30–1.00.
* **Quasi-binomial GLM on the bridge fraction, main effects only** (10 points leave no room for
  curvature or interactions). Odds ratios per coded unit, with 95 % CIs:
  spacing ×2.6 [1.7, 4.1], p < 1e−4; density ×1.9 [1.1, 3.3], p = 0.02; length ×2.2 [0.9, 5.4],
  p = 0.09; salt ×1.5 [0.6, 3.8], p = 0.44. Dispersion 1.35. Lack of fit against the
  replicates is not significant (F(5, 20) = 1.09, p = 0.40).
* **Spacing has the opposite sign here.** It lowers non-reactive inter-chain contacts (§7) but
  raises the bridge fraction. Wider spacing removes loops faster than it removes bridges, as the
  secondary metrics already showed (loops: spacing 54 % of R²).
* **Proxy check (§9.6).** At the 10 design points, the non-reactive model's predictions (our 104 runs)
  compare with the reactive bridge fraction as follows:

  | non-reactive prediction | Spearman | Pearson |
  |---|---|---|
  | inter-chain free contacts (primary response) | 0.32 (p = 0.37) | 0.45 |
  | bridge share = inter / (inter + intra) free contacts | **0.90** (p = 3e−4) | **0.97** |

  The bridge share matches even in absolute value (A10 0.36 predicted vs 0.39 observed; A12
  0.78 vs 0.79). If crosslinking saturates in the experiment too, the loop-versus-bridge
  competition sets the network, and the raw inter-chain contact count is the weaker proxy.
  This argues for promoting the bridge share to a co-primary response. That is a change to the
  frozen plan, for Chris and the team to decide. Caveats: 10 points, a teammate's table not
  re-derived here, model predictions rather than paired non-reactive runs, and bridge fraction
  is not stiffness either. The non-reactive A10–A19 reruns would give a paired comparison.

### 7d. Non-reactive A10–A15 reruns (added 2026-10-17) — `-nr` runs, `stats/`, `export/models/`

18 runs, 6 anchors × 3, same design points as Titus's reactive A10–A15 but crosslinking off, so
they are proper replicates for the primary analysis.

* **Out-of-sample check.** The quadratic fitted on the earlier 104 runs predicts the 6 new anchor
  means with R² 0.90, RMSE 0.16 log units (below the 0.25 grouped-CV RMSE), mean residual +0.08;
  17/18 runs lie inside the 95 % prediction interval. Per anchor, observed / predicted:
  A10 1.33, A11 1.08, A12 0.98, A13 1.17, A14 1.05, A15 0.95. A10, at C = 0.014 in the
  low-density corner, is again the worst, and is under-predicted, the opposite direction to
  the over-prediction at the narrow-spacing corner (`extra2-05`, `extra-04`).
* **Refit on 122 runs:** R² 0.959, adj. 0.953, grouped-CV R² 0.890 (RMSE 0.23), RF CV R² 0.862.
  Variance shares: C 83 %, L 6 %, S 5 %, I 6 % (as before). Lack of fit F(75, 32) = 4.9,
  p = 2e−6 (was F 6.6 with df_pe = 20); sd of pure error 0.078 vs residual sd 0.151. The
  quadratic is still not adequate. Log-coded C: R² 0.981, CV 0.936, LOF F 2.1, p = 0.011.
* **Paired reactive vs non-reactive, now measured instead of predicted.** Team-A inter-chain
  contacts per frame, reactive / non-reactive, are ×38–91 across A10–A15 (reactive 7.5–20 per
  frame, non-reactive 0.12–0.49). This confirms the exclusion
  from the contact fits.
* **Proxy check on the 6 paired points.** The reactive bridge fraction against the *measured*
  non-reactive quantity: contacts Spearman 0.83 (p = 0.04), bridge share inter / (inter + intra)
  Spearman 0.83, Pearson **0.97** (p = 0.001). With six points that is suggestive only. The
  bridge share matches in absolute value (A12 0.79 vs 0.79, A15 0.92 vs 0.92, A10 0.45 vs 0.39),
  while the contact count is rank-consistent over these six points. §7b found weak rank agreement for
  the contact count over all 10 anchors, using model predictions; whether that persists with
  measured values has to wait for A16–A19.
* **QC:** 18/18 completed, 5714 frames, one leg. No drift flag at |z| > 3 (max |z| 2.6, A12-rep1);
  4 runs have a late Rg/height burn-in detection (flagged, kept). One run, A10-rep2, has Cook's D
  above 4/n (studentised residual +2.6, observed 0.140 vs 0.106 / 0.094 for its siblings).

### 7c. Log-coded refit and two extra responses (added 2026-10-16) — `export/models/`, `scripts/07_models.py`

**Not frozen: the log coding was chosen after the frozen model's lack of fit was seen.** *The table
below is the 104-run version; the 122-run values are in `export/models/comparison.csv` (inter:
linear R² 0.959, CV 0.890, LOF F 4.9; log R² 0.981, CV 0.936, LOF F 2.1, p = 0.011).* The frozen
linear-coded model stays the baseline. Same full quadratic, n = 104, with the natural log of each
parameter coded to −1…+1 over its log range:

| response | coding | R² | grouped-CV R² | lack of fit F, p |
|---|---|---|---|---|
| inter | linear | 0.960 | 0.882 | 6.61, 1e-05 |
| inter | log | 0.984 | 0.944 | 2.51, 0.012 |
| intra | linear | 0.986 | 0.976 | 2.14, 0.029 |
| intra | log | 0.989 | 0.982 | 1.61, 0.12 |
| interfrac | linear | 0.978 | 0.939 | 2.62, 0.0092 |
| interfrac | log | 0.988 | 0.962 | 1.23, 0.31 |

* Primary response, log coding: no spurious maximum. Density gives ×8.3 [7.2, 9.5] from 0.01 to
  0.06, against ×5.4 frozen, and rises monotonically. Length ×0.40, spacing ×0.63, salt ×1.25. The
  lack of fit is smaller but still significant (p = 0.012). Without `extra2-05` (Cook's D 1.47 under
  this coding, leverage 0.45) it is not (F 1.69, p = 0.09). The run is QC-clean and is kept.
* `interfrac` (inter / (inter + intra), quasi-binomial GLM, logit link, R² = deviance pseudo-R²): density
  52 %, spacing 27 %, length 19 %, salt ~2 % (n.s.). Longer chains *raise* the fraction (odds ×3.2
  across the range), although they lower the inter-chain count, because they remove loops faster.
* `intra`: spacing 54 %, length 41 %; density and salt small.

## 8. Figures (draft; Plotly HTML with hover + PNG) — `figures/`

Palette check: the five team colours were run through the dataviz validator on white. Yellow
`#FFAF00` fails for thin marks (contrast 1.8:1), so it appears only inside the continuous
surface scale. Green and orange are close for deuteranopes (ΔE 7.2), so the figures never use
them as a pair without a second cue.

| file | proposed caption |
|---|---|
| `coefficients` | Standardised coefficients of the quadratic model on log(inter-chain free K–K contacts per frame), with 95 % CIs; parameters coded to [−1, +1] over their sampled ranges. Blue: CI excludes 0. C grafting density, L chain length, S lysine spacing, I ionic strength. n = 104 runs. |
| `predicted_vs_heldout` | Observed vs held-out predicted proxy, 5-fold cross-validation with the three replicates of each anchor kept in one fold; log axes; dotted line y = x. Orange squares: anchor replicates. |
| `response_surface_base` | Predicted proxy vs grafting density and lysine spacing, with L = 72 repeats and I = 0.20 M (base case). Arrows: direction of steepest increase in coded units. Preliminary: the fall-off at the top right is a quadratic artefact not seen in the raw data (README §7). |
| `response_surface_2x2` | As above, with chain length and ionic strength at 10 % / 90 % of their ranges (L 60 / 140, I 0.095 / 0.455 M). Shared colour scale. Corners of the 4-D space are thinly sampled. |
| `residual_diagnostics` | Studentised residuals vs fitted, normal Q–Q, and Cook's distance (dashed: 4/n; orange: above it). |
| `noise_check_anchors` | Replicate SD and CV vs the anchor mean, 10 anchors × 3 replicates. |
| `importance_comparison` | Parameter importance normalised to sum to 1: LMG share of R² (quadratic), held-out permutation importance and mean \|SHAP\| (random forest). |
| `qc_timeseries_sample_*` | Primary response and chain-mean Rg vs time for 12 randomly chosen runs (seed 20261015). Light: per frame; dark: 100-frame running mean; dashed: per-run Chodera burn-in (the analysis discards a fixed 100 ns). |

## 9. Unresolved — for Chris / the next round

1. **Data:** where are the 100 `sub4-main-*` runs (and extra5 15–20)? Were they ever submitted?
   Without them the design has 84 points, not 220.
2. **A16–A19 reruns:** A10–A15 are in (§7d); A16–A19 (12 runs) are not finished. When they are,
   rerun `analysis/scripts/run_all.sh` and update §7d (pure error then reaches 20 anchors / 40 df
   as planned, and the proxy check covers all 10 paired points). Add the rerun sheet to
   `used-runs.csv` via the notebook as usual.
3. **Titus's 20 `sub4-extra` runs:** now in the team-A fit from his summaries (§7a). To put them in
   the primary (free-only) analysis and QC them, he still needs to share the run folders or run
   step 2a and send the `.npz` files.
4. **Model form:** the lack of fit is real, and Titus's runs confirm both failure spots out of
   sample (§7a). Log-coded C fixes the high-density turnover but not the low-density,
   narrow-spacing corner. Candidates: log C plus a C×S term on the log scale, or a saturating form
   in C. Re-test after (1)–(3).
5. **Low-density corner:** at C ≈ 0.01 the response is ~0.04–0.1 contacts/frame and the model is
   least reliable there. If target tissues land there, it needs more runs.
6. **Proxy validity:** the proxy is the instantaneous contact count of a non-reactive brush.
   Reactive runs (like Titus's anchors) lock in the first encounters instead. In those runs the
   bonds saturate, and the bridge fraction tracks the non-reactive bridge share, not the contact
   count (§7b). Decide whether the bridge share becomes co-primary. Rheology settles it later.
7. **Block position:** the 4-lysine block is centred at 38–50 % of the chain (never past the
   middle), and the lead/trail split is tied to spacing (r = −0.91). Spacing effects therefore
   include a small shift of the block toward the N-terminus. The lysine-arrangement submodule
   found placement matters a lot, so state this when generalising.

**Known issue (found 2026-10-16):** `01_extract_timeseries.py` looked residue masses up by three-letter
code, and placeholder rows in `residues_CALVADOS2.csv` (Z = VAL, X = LYS, MW −2) overwrote valine and
lysine. The mass-weighted heights (`h50_nm`, `h95_nm`, `zmean_nm`) and `rg_nm` in the current
`.npz` files carry this error. Recomputed on three runs it moves them by ≤ 0.05 nm (≤ 0.4 %), so no
conclusion changes. Contacts, Ree and everything in the primary analysis are unaffected. The script is
fixed; regenerate with `01_extract_timeseries.py --overwrite` at the next full rerun, so old and new
runs are consistent.

## 10. Reproduce

From the repository root on DelftBlue (the analysis packages are layered on the project venv
by `uv run --with`; `pyproject.toml` and `uv.lock` are untouched):

```bash
python analysis/scripts/00_inventory.py                    # run_inventory.csv
sbatch analysis/scripts/extract.sbatch                     # step 2a -> timeseries/*.npz (~15 min, 16 cores)
#   or locally: python analysis/scripts/01_extract_timeseries.py --workers 8
analysis/scripts/run.sh analysis/scripts/02_equilibration.py   # qc/equilibration_per_run.csv, QC plots
analysis/scripts/run.sh analysis/scripts/03_metrics.py         # runs_metrics.csv, per_lysine.csv, timeseries_primary.csv.gz
analysis/scripts/run.sh analysis/scripts/04_stats.py           # stats/*  (~5 min, SHAP + 1000-tree forests)
analysis/scripts/run.sh analysis/scripts/05_figures.py         # figures/*
analysis/scripts/run.sh analysis/scripts/06_reactive.py        # stats/reactive_* (Titus's A10-A19)
analysis/scripts/run.sh analysis/scripts/07_models.py          # export/models/* (frozen + log-coded, intra, interfrac)
analysis/scripts/run.sh analysis/scripts/05_figures.py --stats analysis/export/models/inter_log \
    --out analysis/export/figures_logcoded --coding log           # export/figures_logcoded/*
analysis/scripts/run.sh analysis/scripts/08_video_runs.py      # export/video_runs.md, export/video_runs/*
analysis/scripts/run.sh analysis/scripts/09_export.py          # export/ copies, MANIFEST.md, README.md, zip
# or all of it:  analysis/scripts/run_all.sh
```

`python` = `.venv/bin/python`. `run.sh` pins scikit-learn 1.5.2, plotly 5.24.1, kaleido 0.2.1 and shap 0.46.0.
For PNGs in Poppins, the font has to be in `~/.fonts` (kaleido does not read `~/.local/share/fonts`).
When new runs arrive, add their sheet to `LOCAL_SHEETS` in `scripts/common.py` (or
`TEAMMATE_SHEETS` if only the notebook summaries are available) and rerun from the top. Step 2a
only extracts runs that have no `.npz` yet. If Titus's trajectories become readable, the
`summary_only` rows switch to `finished` and go through the full pipeline automatically.

Scripts: `common.py` (paths, sheet parsing, coding), `model.py` (the shared model definition),
`style.py` (palette, font, layout).
