"""Step 5 - first-pass statistics (preliminary). Writes analysis/stats/.

  coefficients.csv        OLS (log response) coefficients, SE, t, p, 95 % CI, VIF
  parameter_tests.csv     per-parameter F-test (drop every term containing it), unique dR2, LMG share
  lack_of_fit.csv         lack-of-fit F-test against the anchor replicates' pure error
  model_fit.csv           R2, adj R2, grouped-CV R2/RMSE for OLS and random forest, diagnostics
  run_diagnostics.csv     per run: fitted, residual, leverage, Cook's D, CV prediction
  importance_rf.csv       random-forest permutation importance and mean |SHAP|
  glm_gamma.csv           Gamma GLM (log link) coefficients, the variance-grows-with-mean check
  sensitivity.csv         the same model on the frozen sensitivity responses (team definition A
                          twice: own 104 runs, and with the teammate's 20 summary-only runs)
  teammate_check.csv      teammate runs vs the team-A model fitted on our runs alone (held out)
  reactive_anchors.csv    teammate's reactive anchors A10-A19 vs the non-reactive team-A model
  secondary.csv           the same model on the secondary metrics (R2 and LMG shares only)
  summary.txt             everything above in plain text

    analysis/scripts/run.sh analysis/scripts/04_stats.py
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import yaml
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor

sys.path.insert(0, str(Path(__file__).parent))
import model as M  # noqa: E402
from common import ANALYSIS  # noqa: E402

SEED = yaml.safe_load((ANALYSIS / "settings.yaml").read_text())["stats"]["random_seed"]
OUT = ANALYSIS / "stats"
OUT.mkdir(exist_ok=True)
X_COLS = [M.SHORT[p] for p in M.PARAMS]


def lmg_shares(df, y="logy"):
    """Shapley (LMG) decomposition of R2 over the four parameter groups; sums to R2."""

    def r2(subset):
        if not subset:
            return 0.0
        # a term enters when ALL of its variables are in the subset
        allowed = {M.SHORT[p] for p in subset}
        terms = [t for t in M.TERMS if M.term_vars(t) <= allowed]
        return M.fit(df, y, terms).rsquared

    cache = {}
    share = {}
    for p in M.PARAMS:
        others = [q for q in M.PARAMS if q != p]
        tot = 0.0
        for k in range(len(others) + 1):
            for sub in itertools.combinations(others, k):
                w = 1 / (len(M.PARAMS) * len(list(itertools.combinations(others, k))))
                a = cache.setdefault(frozenset(sub), r2(list(sub)))
                b = cache.setdefault(frozenset(sub + (p,)), r2(list(sub) + [p]))
                tot += w * (b - a)
        share[p] = tot
    return share


def grouped_cv(df, make_pred, n_splits=5):
    from sklearn.model_selection import GroupKFold
    pred = np.full(len(df), np.nan)
    for tr, te in GroupKFold(n_splits).split(df, groups=df["design_point"]):
        pred[te] = make_pred(df.iloc[tr], df.iloc[te])
    resid = df["logy"].to_numpy() - pred
    r2 = 1 - (resid ** 2).sum() / ((df["logy"] - df["logy"].mean()) ** 2).sum()
    return pred, r2, float(np.sqrt((resid ** 2).mean()))


def lack_of_fit(df, res):
    rep = df[df.anchor_id.notna()]
    ss_pe = float(((rep["logy"] - rep.groupby("anchor_id")["logy"].transform("mean")) ** 2).sum())
    df_pe = int((rep.groupby("anchor_id").size() - 1).sum())
    ss_res, df_res = float(res.ssr), int(res.df_resid)
    ss_lof, df_lof = ss_res - ss_pe, df_res - df_pe
    F = (ss_lof / df_lof) / (ss_pe / df_pe)
    return {"ss_pure_error": ss_pe, "df_pure_error": df_pe, "sd_pure_error_log": np.sqrt(ss_pe / df_pe),
            "ss_lack_of_fit": ss_lof, "df_lack_of_fit": df_lof, "sd_residual_log": np.sqrt(ss_res / df_res),
            "F": F, "p": float(stats.f.sf(F, df_lof, df_pe)),
            "n_design_points": int(df["design_point"].nunique())}


def param_tests(df, res, y="logy"):
    rows = []
    shares = lmg_shares(df, y)
    for p in M.PARAMS:
        red = M.fit(df, y, M.terms_without(p))
        F, pval, ddf = res.compare_f_test(red)
        rows.append({"parameter": p, "terms_dropped": len(M.TERMS) - len(M.terms_without(p)),
                     "F": F, "df_num": ddf, "df_den": res.df_resid, "p": pval,
                     "unique_dR2": res.rsquared - red.rsquared, "lmg_share_R2": shares[p],
                     "lmg_share_of_explained": shares[p] / res.rsquared})
    return pd.DataFrame(rows)


def main():
    df = M.load()
    res = M.fit(df)
    txt = [f"SubM4 first-pass statistics - n = {len(df)} runs, {df.design_point.nunique()} design points, "
           f"response = log({M.RESPONSE})", ""]

    # ---- coefficients, VIF ------------------------------------------------------------
    ci = res.conf_int()
    X = res.model.exog
    vif = {name: (variance_inflation_factor(X, i) if name != "Intercept" else np.nan)
           for i, name in enumerate(res.model.exog_names)}
    coef = pd.DataFrame({"term": res.params.index, "coef_log": res.params.values, "se": res.bse.values,
                         "t": res.tvalues.values, "p": res.pvalues.values, "ci_low": ci[0].values,
                         "ci_high": ci[1].values, "vif": [vif[n] for n in res.params.index]})
    coef["fold_change_per_unit_coded"] = np.exp(coef.coef_log)
    coef.to_csv(OUT / "coefficients.csv", index=False)
    txt += [str(res.summary()), ""]

    # ---- per-parameter tests and variance shares ---------------------------------------
    pt = param_tests(df, res)
    pt.to_csv(OUT / "parameter_tests.csv", index=False)
    txt += ["Per-parameter F-tests (drop every term containing the parameter) and R2 shares:",
            pt.round(4).to_string(index=False), ""]

    # ---- lack of fit ------------------------------------------------------------------
    lof = lack_of_fit(df, res)
    pd.DataFrame([lof]).to_csv(OUT / "lack_of_fit.csv", index=False)
    txt += ["Lack-of-fit F-test (pure error = anchor replicates):",
            pd.Series(lof).round(5).to_string(), ""]

    # ---- diagnostics ------------------------------------------------------------------
    infl = res.get_influence()
    cooks = infl.cooks_distance[0]
    sw = stats.shapiro(res.resid)
    bp = sm.stats.het_breuschpagan(res.resid, res.model.exog)

    # ---- grouped CV: OLS and random forest --------------------------------------------
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.inspection import permutation_importance

    ols_pred, ols_cv_r2, ols_cv_rmse = grouped_cv(df, lambda tr, te: M.fit(tr).predict(te).to_numpy())

    def rf():
        return RandomForestRegressor(n_estimators=1000, min_samples_leaf=2, random_state=SEED, n_jobs=4)

    rf_pred, rf_cv_r2, rf_cv_rmse = grouped_cv(
        df, lambda tr, te: rf().fit(tr[X_COLS], tr["logy"]).predict(te[X_COLS]))

    rf_full = rf().fit(df[X_COLS], df["logy"])
    perm = permutation_importance(rf_full, df[X_COLS], df["logy"], n_repeats=50, random_state=SEED)
    # held-out permutation importance, grouped folds
    from sklearn.model_selection import GroupKFold
    perm_cv = np.zeros(len(X_COLS))
    for tr, te in GroupKFold(5).split(df, groups=df["design_point"]):
        mdl = rf().fit(df.iloc[tr][X_COLS], df.iloc[tr]["logy"])
        perm_cv += permutation_importance(mdl, df.iloc[te][X_COLS], df.iloc[te]["logy"], n_repeats=50,
                                          random_state=SEED).importances_mean / 5
    import shap
    sv = shap.TreeExplainer(rf_full).shap_values(df[X_COLS])
    imp = pd.DataFrame({"parameter": M.PARAMS, "perm_importance_train": perm.importances_mean,
                        "perm_importance_train_sd": perm.importances_std,
                        "perm_importance_heldout": perm_cv, "mean_abs_shap": np.abs(sv).mean(0)})
    imp.to_csv(OUT / "importance_rf.csv", index=False)
    np.save(OUT / "shap_values.npy", sv)

    fit_tab = {"n_runs": len(df), "n_params": int(res.df_model + 1), "R2": res.rsquared,
               "adj_R2": res.rsquared_adj, "cv5_grouped_R2_ols": ols_cv_r2, "cv5_grouped_RMSE_log_ols": ols_cv_rmse,
               "cv5_grouped_R2_rf": rf_cv_r2, "cv5_grouped_RMSE_log_rf": rf_cv_rmse,
               "shapiro_W": sw.statistic, "shapiro_p": sw.pvalue, "breusch_pagan_p": bp[1],
               "max_vif": np.nanmax(list(vif.values())), "n_cooks_gt_4_over_n": int((cooks > 4 / len(df)).sum()),
               "max_cooks": float(cooks.max())}
    pd.DataFrame([fit_tab]).to_csv(OUT / "model_fit.csv", index=False)
    txt += ["Fit, cross-validation and diagnostics:", pd.Series(fit_tab).round(4).to_string(), "",
            "Random forest importance (log response):", imp.round(4).to_string(index=False), ""]

    diag = df[["slug", "design_point", "anchor_id", "replicate", "concentration", "chain_length",
               "spacing_pct", "ionic_strength", "y", "logy"]].copy()
    diag["fitted"] = res.fittedvalues
    diag["resid"] = res.resid
    diag["resid_studentized"] = infl.resid_studentized_external
    diag["leverage"] = infl.hat_matrix_diag
    diag["cooks_d"] = cooks
    diag["qc_cooks_gt_4_over_n"] = cooks > 4 / len(df)
    diag["cv_pred_ols"] = ols_pred
    diag["cv_pred_rf"] = rf_pred
    diag.to_csv(OUT / "run_diagnostics.csv", index=False)
    txt += ["Runs with Cook's D > 4/n:",
            diag.loc[diag.qc_cooks_gt_4_over_n, ["slug", "y", "resid_studentized", "leverage", "cooks_d"]]
            .round(3).to_string(index=False), ""]

    # ---- Gamma GLM (log link): variance grows with the mean -----------------------------
    glm = smf.glm(M.formula(y="y"), data=df, family=sm.families.Gamma(sm.families.links.Log())).fit()
    g = pd.DataFrame({"term": glm.params.index, "coef_log": glm.params.values, "se": glm.bse.values,
                      "p": glm.pvalues.values, "ci_low": glm.conf_int()[0].values,
                      "ci_high": glm.conf_int()[1].values, "ols_coef_log": res.params.reindex(glm.params.index).values})
    g.to_csv(OUT / "glm_gamma.csv", index=False)
    txt += ["Gamma GLM (log link) vs OLS on log y - coefficients:", g.round(4).to_string(index=False), ""]

    # ---- sensitivity responses (frozen before the fit) -----------------------------------
    sens = []
    for r, own in [("inter_free_12A", True), ("inter_free_8A", True), ("inter_free_16A", True),
                   ("inter_teamA_12A", True), ("inter_teamA_12A", False), ("inter_free_12A_from200ns", True)]:
        d = M.load(r, own_only=own)
        rr = M.fit(d)
        p2 = param_tests(d, rr)
        row = {"response": r, "runs": "own" if own else "own + teammate", "n": len(d),
               "R2": rr.rsquared, "adj_R2": rr.rsquared_adj}
        for _, q in p2.iterrows():
            row[f"lmg_{M.SHORT[q.parameter]}"] = q.lmg_share_R2
            row[f"p_{M.SHORT[q.parameter]}"] = q.p
        for t in X_COLS:
            row[f"coef_{t}"] = rr.params[t]
        sens.append(row)
    sens = pd.DataFrame(sens)
    sens.to_csv(OUT / "sensitivity.csv", index=False)
    txt += ["Sensitivity (same model, other frozen responses):", sens.round(4).to_string(index=False), ""]

    # ---- teammate runs (team definition A only) -------------------------------------------
    own_a = M.load("inter_teamA_12A", own_only=True)
    all_a = M.load("inter_teamA_12A")
    res_a = M.fit(own_a)
    tm = all_a[all_a.data_source == "teammate_summary"].copy()
    pr = res_a.get_prediction(tm).summary_frame(alpha=0.05)
    tm["pred_log"] = pr["mean"].to_numpy()
    tm["pi_low_log"], tm["pi_high_log"] = pr["obs_ci_lower"].to_numpy(), pr["obs_ci_upper"].to_numpy()
    tm["resid_log"] = tm.logy - tm.pred_log
    tm["outside_95pi"] = (tm.logy < tm.pi_low_log) | (tm.logy > tm.pi_high_log)
    tm[["slug", "concentration", "chain_length", "spacing_pct", "ionic_strength", "y", "pred_log",
        "pi_low_log", "pi_high_log", "resid_log", "outside_95pi"]].assign(pred=np.exp(tm.pred_log)) \
        .to_csv(OUT / "teammate_check.csv", index=False)
    rss = float((tm.resid_log ** 2).sum())
    tcheck = {"n_teammate": len(tm),
              "heldout_R2": 1 - rss / float(((tm.logy - tm.logy.mean()) ** 2).sum()),
              "heldout_RMSE_log": np.sqrt(rss / len(tm)),
              "own_cv5_RMSE_log": grouped_cv(own_a, lambda tr, te: M.fit(tr).predict(te).to_numpy())[2],
              "mean_resid_log": tm.resid_log.mean(),
              "mean_resid_p": stats.ttest_1samp(tm.resid_log, 0).pvalue,
              "n_outside_95pi": int(tm.outside_95pi.sum())}
    res_all = M.fit(all_a)
    lof_a = {f"own_{k}": v for k, v in lack_of_fit(own_a, res_a).items() if k in ("F", "p", "sd_residual_log")}
    lof_a.update({f"all_{k}": v for k, v in lack_of_fit(all_a, res_all).items() if k in ("F", "p", "sd_residual_log")})
    txt += ["Teammate runs (team definition A, summary only) predicted by the team-A model fitted on our",
            "runs alone (held out), then lack of fit of the team-A model without / with them:",
            pd.Series({**tcheck, **lof_a}).round(4).to_string(),
            tm.loc[tm.outside_95pi | (tm.resid_log.abs() > 0.5), ["slug", "y", "resid_log"]]
            .assign(pred=np.exp(tm.pred_log)).round(3).to_string(index=False), ""]

    # ---- reactive anchors A10-A19 (excluded): how far crosslinking moves the count ------------
    rx = pd.read_csv(ANALYSIS / "reactive_anchors.csv")
    rx = pd.concat([rx, M.coded(rx).rename(columns=M.SHORT)], axis=1)
    rx["pred_nonreactive"] = np.exp(res_all.predict(rx))
    rx["ratio_obs_over_pred"] = rx.inter_teamA_12A / rx.pred_nonreactive
    rx.drop(columns=X_COLS).to_csv(OUT / "reactive_anchors.csv", index=False)
    g = rx.groupby("anchor_id").agg(obs_mean=("inter_teamA_12A", "mean"), obs_sd=("inter_teamA_12A", "std"),
                                    pred_nonreactive=("pred_nonreactive", "first"),
                                    ratio=("ratio_obs_over_pred", "mean"))
    txt += ["Reactive anchors A10-A19 (crosslinking on, excluded from every fit) vs the non-reactive",
            "team-A model on all runs:",
            g.round(3).to_string(), ""]

    # ---- secondary metrics (same model, descriptive) -------------------------------------
    sec = []
    for r, logit in [("intra_free_12A", True), ("bridge_loop_ratio_12A", True), ("largest_cluster_frac", False),
                     ("h95_nm", True), ("rg_nm", True), ("ree_nm", True), ("inter_free_pairs_ever", False)]:
        d = M.load(r)
        if logit:
            if (d["y"] <= 0).any():
                continue
        else:
            d["logy"] = d["y"]
        rr = M.fit(d)
        sh = lmg_shares(d)
        sec.append({"metric": r, "scale": "log" if logit else "linear", "R2": rr.rsquared,
                    **{f"lmg_{M.SHORT[p]}": sh[p] for p in M.PARAMS}})
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / "secondary.csv", index=False)
    txt += ["Secondary metrics (same quadratic model; LMG R2 shares):", sec.round(3).to_string(index=False)]

    (OUT / "summary.txt").write_text("\n".join(txt))
    print("\n".join(txt[-60:]))


if __name__ == "__main__":
    main()
