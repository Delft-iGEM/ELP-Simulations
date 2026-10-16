"""Step 6 - the teammate's reactive anchors A10-A19 (crosslinking on), analysed on their own.

These 30 runs (10 design points x 3) cannot join the contact fits: bonded lysine pairs are
held at 0.6 nm, so the contact count measures the bonds, not the opportunity. What they do
measure is the crosslinked state itself. Conversion is saturated (20-22 of 22 possible bonds),
so the informative response is the bridge fraction: inter-chain bonds / all bonds.

  reactive_design_points.csv  per anchor: bridge fraction mean/SD, and the non-reactive models'
                              predictions at the same design point
  reactive_glm.csv            quasi-binomial GLM (logit) of bridges / formed on the coded parameters
  reactive_summary.txt        everything above in plain text

Ten design points support main effects only (5 parameters); there is no room for curvature
or interactions, so this is a direction-of-effect check, not a response surface.

    analysis/scripts/run.sh analysis/scripts/06_reactive.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor

sys.path.insert(0, str(Path(__file__).parent))
import model as M  # noqa: E402
from common import ANALYSIS  # noqa: E402

OUT = ANALYSIS / "stats"
X_COLS = [M.SHORT[p] for p in M.PARAMS]


def glm(df, rhs):
    return smf.glm(f"xl_bridge_frac ~ {rhs}", data=df, family=sm.families.Binomial(),
                   var_weights=df["xl_formed"]).fit(scale="X2")


def main():
    rx = pd.read_csv(ANALYSIS / "reactive_anchors.csv")
    rx = pd.concat([rx, M.coded(rx).rename(columns=M.SHORT)], axis=1)
    txt = [f"SubM4 reactive anchors (teammate, crosslinking on) - n = {len(rx)} runs, "
           f"{rx.anchor_id.nunique()} design points", "",
           f"Crosslinks formed: {rx.xl_formed.min()}-{rx.xl_formed.max()} of {rx.xl_max.iloc[0]} possible "
           f"(conversion {rx.xl_conversion.min():.2f}-{rx.xl_conversion.max():.2f}); "
           f"bridge fraction {rx.xl_bridge_frac.min():.2f}-{rx.xl_bridge_frac.max():.2f}", ""]

    # ---- main-effects quasi-binomial GLM, lack of fit against the replicates ------------------
    lin = glm(rx, " + ".join(X_COLS))
    sat = glm(rx, "anchor_id")                          # one mean per design point = pure error (string column -> factor)
    df_pe = int(sat.df_resid)
    pe_scale = sat.pearson_chi2 / df_pe
    F_lof = ((lin.deviance - sat.deviance) / (lin.df_resid - df_pe)) / (sat.deviance / df_pe)
    p_lof = float(stats.f.sf(F_lof, lin.df_resid - df_pe, df_pe))
    X = lin.model.exog
    vif = [np.nan] + [variance_inflation_factor(X, i) for i in range(1, X.shape[1])]
    ci = lin.conf_int()
    coef = pd.DataFrame({"term": lin.params.index, "coef_logit": lin.params.values, "se": lin.bse.values,
                         "p": lin.pvalues.values, "ci_low": ci[0].values, "ci_high": ci[1].values,
                         "odds_ratio_per_unit_coded": np.exp(lin.params.values), "vif": vif})
    coef.to_csv(OUT / "reactive_glm.csv", index=False)
    txt += ["Quasi-binomial GLM, logit(bridges / formed) ~ C + L + S + ion (coded [-1, +1]):",
            coef.round(4).to_string(index=False), "",
            f"dispersion (Pearson X2/df): main effects {lin.scale:.2f}, replicates only {pe_scale:.2f} "
            f"(1 = binomial; > 1 = replicates differ more than coin-flipping would)",
            f"lack of fit vs replicate means: F({int(lin.df_resid - df_pe)}, {df_pe}) = {F_lof:.2f}, p = {p_lof:.3g}",
            ""]

    # ---- does the non-reactive proxy rank these design points like the crosslinks do? ----------
    dp = rx.groupby("anchor_id").agg(
        concentration=("concentration", "first"), chain_length=("chain_length", "first"),
        spacing_pct=("spacing_pct", "first"), ionic_strength=("ionic_strength", "first"),
        bridge_frac_mean=("xl_bridge_frac", "mean"), bridge_frac_sd=("xl_bridge_frac", "std"),
        formed_mean=("xl_formed", "mean"), **{c: (c, "first") for c in X_COLS})
    prim = M.fit(M.load())                                        # free-only contacts, our 104 runs
    blr = M.fit(M.load("bridge_loop_ratio_12A"))                  # inter / intra free contacts
    dp["pred_inter_free_12A"] = np.exp(prim.predict(dp))
    r = np.exp(blr.predict(dp))
    dp["pred_bridge_share"] = r / (1 + r)                         # inter / (inter + intra) contacts
    dp["glm_bridge_frac"] = lin.predict(dp)
    dp.drop(columns=X_COLS).to_csv(OUT / "reactive_design_points.csv")
    rows = []
    for col, what in [("pred_inter_free_12A", "inter-chain contacts (primary proxy)"),
                      ("pred_bridge_share", "inter / (inter + intra) contacts")]:
        rho, p = stats.spearmanr(dp[col], dp.bridge_frac_mean)
        pr, pp = stats.pearsonr(dp[col], dp.bridge_frac_mean)
        rows.append({"non_reactive_prediction": what, "spearman": rho, "spearman_p": p,
                     "pearson": pr, "pearson_p": pp})
    proxy = pd.DataFrame(rows)
    txt += ["Per design point (bridge fraction observed in the reactive runs; non-reactive model predictions",
            "at the same parameters, from our 104 trajectories):",
            dp.drop(columns=X_COLS).round(3).to_string(), "",
            "Rank agreement over the 10 design points:", proxy.round(3).to_string(index=False)]

    (OUT / "reactive_summary.txt").write_text("\n".join(txt))
    print("\n".join(txt))


if __name__ == "__main__":
    main()
