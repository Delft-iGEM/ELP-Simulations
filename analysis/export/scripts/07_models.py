"""Step 7 - the frozen primary model next to a log-coded refit, plus two extra responses.

Writes analysis/export/models/<response>_<coding>/ for

  responses  inter      log(inter_free_12A), OLS - the frozen primary response
             intra      log(intra_free_12A), OLS
             interfrac  inter / (inter + intra) free contacts at 12 A: quasi-binomial GLM, logit
                        link, on the production pair-frame counts (inter + intra) x 4286 frames.
                        The dispersion is estimated (Pearson X2 / df), so the inference does not
                        treat autocorrelated frames as independent trials.
  codings    linear     frozen (settings.yaml): each parameter coded to [-1, +1] over its range
             log        natural log of each parameter, coded to [-1, +1] over the log of its
                        range. NOT frozen: chosen on 2026-10-16, after the lack of fit of the
                        frozen model had been seen. The frozen model stays the baseline.

Every folder has the stats/ file schema (05_figures.py can draw the OLS ones):
  coefficients.csv     coefficient, SE, t (OLS) / z (GLM), p, 95 % CI, VIF
  parameter_tests.csv  per-parameter F-test (drop every term with the parameter) and LMG shares
  lack_of_fit.csv      lack-of-fit F-test against the pure error of the triplicate anchors
  model_fit.csv        R2, adjusted R2 (deviance-based for the GLM), 5-fold grouped CV,
                       residual normality / heteroscedasticity, max VIF, Cook's distance
  run_diagnostics.csv  per run: fitted, residual, studentised residual, leverage, Cook's D,
                       held-out prediction
  effects_oat.csv      one parameter at a time, others at the base case: predictions at the
                       low end, base and high end of the range with 95 % CIs, and high/low ratio
plus models/comparison.csv and models/lmg_shares.csv across all six.

    analysis/scripts/run.sh analysis/scripts/07_models.py
"""
from __future__ import annotations

import importlib
import itertools
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.special import expit, logit
from statsmodels.stats.outliers_influence import variance_inflation_factor

sys.path.insert(0, str(Path(__file__).parent))
import model as M  # noqa: E402
from common import ANALYSIS, BASE, RANGES  # noqa: E402

ST = importlib.import_module("04_stats")      # lack_of_fit, param_tests, grouped_cv: the frozen code
OUT = ANALYSIS / "export" / "models"
RESPONSES = {"inter": "inter_free_12A", "intra": "intra_free_12A", "interfrac": None}
CODING_NOTE = {"linear": "frozen coding (settings.yaml)",
               "log": "log coding, chosen after the frozen model's lack of fit was seen (not frozen)"}


# ---- data -------------------------------------------------------------------------------
def load(resp: str, coding: str) -> pd.DataFrame:
    if RESPONSES[resp]:
        return M.load(RESPONSES[resp], coding=coding)
    df = M.load("intra_free_12A", coding=coding)            # trajectory runs: inter and intra both present
    tot = df.inter_free_12A + df.intra_free_12A
    df["n_hits"] = (tot * df.production_frames).round()     # pair-frames in production, exact integers
    df["y"] = df.inter_free_12A / tot
    df["logy"] = logit(df["y"])                             # link scale, for CV and plots
    return df


# ---- the quasi-binomial GLM, mirroring the OLS helpers in 04_stats -----------------------
def glm_fit(df, terms=M.TERMS):
    return smf.glm(M.formula(terms, "y"), data=df, family=sm.families.Binomial(),
                   var_weights=df["n_hits"]).fit(scale="X2")


def pseudo_r2(res):
    return 1 - res.deviance / res.null_deviance


def glm_lmg(df):
    def r2(subset):
        if not subset:
            return 0.0
        allowed = {M.SHORT[p] for p in subset}
        return pseudo_r2(glm_fit(df, [t for t in M.TERMS if M.term_vars(t) <= allowed]))

    cache, share = {}, {}
    for p in M.PARAMS:
        others = [q for q in M.PARAMS if q != p]
        tot = 0.0
        for k in range(len(others) + 1):
            subs = list(itertools.combinations(others, k))
            for sub in subs:
                a = cache.setdefault(frozenset(sub), r2(list(sub)))
                b = cache.setdefault(frozenset(sub + (p,)), r2(list(sub) + [p]))
                tot += (b - a) / (len(M.PARAMS) * len(subs))
        share[p] = tot
    return share


def glm_param_tests(df, res):
    shares = glm_lmg(df)
    rows = []
    for p in M.PARAMS:
        red = glm_fit(df, M.terms_without(p))
        q = int(red.df_resid - res.df_resid)
        F = (red.deviance - res.deviance) / q / res.scale
        rows.append({"parameter": p, "terms_dropped": q, "F": F, "df_num": q, "df_den": res.df_resid,
                     "p": float(stats.f.sf(F, q, res.df_resid)),
                     "unique_dR2": pseudo_r2(res) - pseudo_r2(red), "lmg_share_R2": shares[p],
                     "lmg_share_of_explained": shares[p] / pseudo_r2(res)})
    return pd.DataFrame(rows)


def glm_lack_of_fit(df, res):
    pe = smf.glm("y ~ design_point", data=df, family=sm.families.Binomial(), var_weights=df["n_hits"]).fit()
    df_pe, df_res = int(pe.df_resid), int(res.df_resid)
    F = ((res.deviance - pe.deviance) / (df_res - df_pe)) / (pe.deviance / df_pe)
    return {"ss_pure_error": pe.deviance, "df_pure_error": df_pe,
            "dispersion_pure_error": pe.pearson_chi2 / df_pe,
            "ss_lack_of_fit": res.deviance - pe.deviance, "df_lack_of_fit": df_res - df_pe,
            "dispersion_model": res.scale, "F": F, "p": float(stats.f.sf(F, df_res - df_pe, df_pe)),
            "n_design_points": int(df.design_point.nunique()), "scale": "deviance"}


# ---- shared pieces -------------------------------------------------------------------------
def effects_oat(res, coding, inv_link, ratio_kind, use_t):
    """Predictions one parameter at a time (others at the base case), 95 % CIs from the
    coefficient covariance; the high/low ratio comes from the contrast of the two rows."""
    di = res.model.data.design_info
    cov = np.asarray(res.cov_params())
    beta = np.asarray(res.params)
    q = stats.t.ppf(0.975, res.df_resid) if use_t else stats.norm.ppf(0.975)

    def row(values):
        d = pd.DataFrame([{M.SHORT[p]: float(M.to_coded(p, values[p], coding)) for p in M.PARAMS}])
        return np.asarray(patsy.build_design_matrices([di], d)[0])[0]

    def pred(x):
        eta, se = x @ beta, np.sqrt(x @ cov @ x)
        return inv_link(eta), inv_link(eta - q * se), inv_link(eta + q * se)

    xb = row(BASE)
    b = pred(xb)
    out = []
    for p in M.PARAMS:
        lo, hi = RANGES[p]
        xl, xh = row({**BASE, p: lo}), row({**BASE, p: hi})
        c = xh - xl
        dlt, se = c @ beta, np.sqrt(c @ cov @ c)
        out.append({"parameter": p, "base_value": BASE[p], "low_value": lo, "high_value": hi,
                    "base_pred": b[0], "base_ci_low": b[1], "base_ci_high": b[2],
                    **dict(zip(["low_pred", "low_ci_low", "low_ci_high"], pred(xl))),
                    **dict(zip(["high_pred", "high_ci_low", "high_ci_high"], pred(xh))),
                    "ratio_kind": ratio_kind, "ratio_high_over_low": np.exp(dlt),
                    "ratio_ci_low": np.exp(dlt - q * se), "ratio_ci_high": np.exp(dlt + q * se)})
    return pd.DataFrame(out)


def vif_table(res):
    X = res.model.exog
    return {n: (variance_inflation_factor(X, i) if n != "Intercept" else np.nan)
            for i, n in enumerate(res.model.exog_names)}


def run_one(resp: str, coding: str):
    out = OUT / f"{resp}_{coding}"
    out.mkdir(parents=True, exist_ok=True)
    df = load(resp, coding)
    is_glm = resp == "interfrac"
    res = glm_fit(df) if is_glm else M.fit(df)
    vif = vif_table(res)
    ci = res.conf_int()

    # coefficients
    coef = pd.DataFrame({"term": res.params.index, ("coef_logit" if is_glm else "coef_log"): res.params.values,
                         "se": res.bse.values, ("z" if is_glm else "t"): res.tvalues.values,
                         "p": res.pvalues.values, "ci_low": ci[0].values, "ci_high": ci[1].values,
                         "vif": [vif[n] for n in res.params.index]})
    coef["odds_ratio_per_unit_coded" if is_glm else "fold_change_per_unit_coded"] = np.exp(res.params.values)
    coef.to_csv(out / "coefficients.csv", index=False)

    # per-parameter tests + LMG, lack of fit
    pt = glm_param_tests(df, res) if is_glm else ST.param_tests(df, res)
    pt.to_csv(out / "parameter_tests.csv", index=False)
    lof = glm_lack_of_fit(df, res) if is_glm else {**ST.lack_of_fit(df, res), "scale": "log"}
    pd.DataFrame([lof]).to_csv(out / "lack_of_fit.csv", index=False)

    # diagnostics and grouped CV (folds as in 04_stats: GroupKFold by design point, no shuffling)
    infl = res.get_influence()
    cooks = np.asarray(infl.cooks_distance[0])
    if is_glm:
        resid = np.asarray(res.resid_deviance)
        bp_resid = np.asarray(res.resid_pearson)
        cv_pred, cv_r2, cv_rmse = ST.grouped_cv(df, lambda tr, te: logit(glm_fit(tr).predict(te).to_numpy()))
        r2 = pseudo_r2(res)
        adj = 1 - (res.deviance / res.df_resid) / (res.null_deviance / (len(df) - 1))
    else:
        resid = bp_resid = np.asarray(res.resid)
        cv_pred, cv_r2, cv_rmse = ST.grouped_cv(df, lambda tr, te: M.fit(tr).predict(te).to_numpy())
        r2, adj = res.rsquared, res.rsquared_adj
    sw = stats.shapiro(resid)
    bp = sm.stats.het_breuschpagan(bp_resid, res.model.exog)
    fit_tab = {"response": resp, "coding": coding, "coding_note": CODING_NOTE[coding],
               "model": "quasi-binomial GLM, logit link, var_weights = pair-frames" if is_glm else "OLS on log response",
               "n_runs": len(df), "n_design_points": int(df.design_point.nunique()),
               "n_params": int(res.df_model + 1),
               "R2": r2, "adj_R2": adj, "R2_kind": "deviance pseudo-R2" if is_glm else "R2",
               "cv5_grouped_R2_ols" if not is_glm else "cv5_grouped_R2": cv_r2,
               "cv5_grouped_RMSE_log_ols" if not is_glm else "cv5_grouped_RMSE_logit": cv_rmse,
               "lof_F": lof["F"], "lof_df_num": lof["df_lack_of_fit"], "lof_df_den": lof["df_pure_error"],
               "lof_p": lof["p"], "dispersion": res.scale if is_glm else np.nan,
               "shapiro_W": sw.statistic, "shapiro_p": sw.pvalue,
               "shapiro_on": "deviance residuals" if is_glm else "residuals",
               "breusch_pagan_p": bp[1], "breusch_pagan_on": "Pearson residuals" if is_glm else "residuals",
               "max_vif": np.nanmax(list(vif.values())),
               "n_cooks_gt_4_over_n": int((cooks > 4 / len(df)).sum()), "max_cooks": float(cooks.max())}
    pd.DataFrame([fit_tab]).to_csv(out / "model_fit.csv", index=False)

    diag = df[["slug", "design_point", "anchor_id", "replicate", "concentration", "chain_length",
               "spacing_pct", "ionic_strength", "y", "logy"] + (["n_hits"] if is_glm else [])].copy()
    diag["fitted"] = np.asarray(res.fittedvalues) if not is_glm else logit(np.asarray(res.fittedvalues))
    diag["resid"] = resid
    diag["resid_studentized"] = (np.asarray(infl.resid_studentized) if is_glm
                                 else np.asarray(infl.resid_studentized_external))
    diag["leverage"] = np.asarray(infl.hat_matrix_diag)
    diag["cooks_d"] = cooks
    diag["qc_cooks_gt_4_over_n"] = cooks > 4 / len(df)
    diag["cv_pred_ols" if not is_glm else "cv_pred_logit"] = cv_pred
    diag.to_csv(out / "run_diagnostics.csv", index=False)

    oat = effects_oat(res, coding, expit if is_glm else np.exp,
                      "odds ratio" if is_glm else "response ratio", use_t=not is_glm)
    oat.to_csv(out / "effects_oat.csv", index=False)

    if resp == "inter":
        # The random forest uses the raw parameters; trees split on order, so a monotone
        # recoding gives the same partitions. Its importances are copied for the comparison figure.
        shutil.copy(ANALYSIS / "stats" / "importance_rf.csv", out / "importance_rf.csv")
    return fit_tab, pt.assign(response=resp, coding=coding)


def main():
    fits, shares = [], []
    for resp in RESPONSES:
        for coding in ("linear", "log"):
            f, s = run_one(resp, coding)
            fits.append(f)
            shares.append(s)
            print(f"{resp:9s} {coding:6s} n={f['n_runs']} R2={f['R2']:.3f} "
                  f"CV R2={f.get('cv5_grouped_R2_ols', f.get('cv5_grouped_R2')):.3f} "
                  f"LOF F={f['lof_F']:.2f} p={f['lof_p']:.2g}")
    comp = pd.DataFrame(fits)
    comp["cv5_R2"] = comp.pop("cv5_grouped_R2_ols").fillna(comp.pop("cv5_grouped_R2"))
    comp["cv5_RMSE_link"] = comp.pop("cv5_grouped_RMSE_log_ols").fillna(comp.pop("cv5_grouped_RMSE_logit"))
    comp.to_csv(OUT / "comparison.csv", index=False)
    sh = pd.concat(shares)[["response", "coding", "parameter", "F", "p", "lmg_share_R2", "lmg_share_of_explained"]]
    sh.to_csv(OUT / "lmg_shares.csv", index=False)

    # the frozen primary model must come out exactly as in stats/
    a = pd.read_csv(ANALYSIS / "stats" / "coefficients.csv")
    b = pd.read_csv(OUT / "inter_linear" / "coefficients.csv")
    assert np.allclose(a.coef_log, b.coef_log) and np.allclose(a.ci_low, b.ci_low), "frozen model not reproduced"
    print("frozen primary model reproduced (coefficients and CIs match stats/)")


if __name__ == "__main__":
    main()
