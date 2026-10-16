"""Step 3 - compare the lysine arrangements and apply the frozen optimality rule.

Writes analysis_lysplace/stats/
  class_summary.csv     per class and response: n, median [bootstrap 95 % CI], IQR, mean, SD
  class_tests.csv       per response: Kruskal-Wallis H, p
  pairwise.csv          per response and class pair: Mann-Whitney U, p, Holm p, Cliff's delta [95 % CI]
  descriptor_glm.csv    quasi-binomial GLM of n_bridges / 22 on the design descriptors (all runs and
                        random only), with 5-fold CV and the test of class beyond the descriptors
  optimal.txt           the optimality rule of settings.yaml applied, step by step
  summary.txt           everything above in plain text

    analysis/scripts/run.sh analysis_lysplace/scripts/03_compare.py
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
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).parent))
from common import CLASSES, MOD, STATS  # noqa: E402

S = yaml.safe_load((MOD / "settings.yaml").read_text())
SEED = int(S["seed"])
NBOOT = int(S["comparison"]["bootstrap"])
ALPHA = float(S["comparison"]["alpha"])
MAXB = int(S["network"]["max_bonds"])
PRIMARY = "n_bridges"
RESPONSES = {   # name: (group, plain label)
    "n_bridges": ("network", "inter-chain crosslinks (bridges)"),
    "n_loops": ("network", "intra-chain crosslinks (loops)"),
    "n_formed": ("network", "crosslinks formed"),
    "bridge_frac": ("network", "bridge fraction"),
    "n_pairs_bridged": ("network", "distinct chain pairs bridged"),
    "secondary_loops": ("network", "extra bridges between already-joined chains"),
    "largest_cluster_frac": ("network", "largest connected cluster / 16 chains"),
    "cycle_rank": ("network", "network cycle rank"),
    "n_unreacted": ("network", "free lysines left unreacted"),
    "n_bridges_100ns": ("network", "bridges at 100 ns"),
    "t50_bridges_ns": ("timing", "time to half the final bridges (ns)"),
    "t90_bridges_ns": ("timing", "time to 90 % of the final bridges (ns)"),
    "h50_nm": ("height", "brush h50 (nm)"),
    "h95_nm": ("height", "brush h95 (nm)"),
    "z_bridge_nm": ("height", "bridge height (nm)"),
    "z_bridge_rel_h95": ("height", "bridge height / h95"),
    "bridge_frac_above_h50": ("height", "bridges above h50 (fraction of bond-frames)"),
    "z_loop_nm": ("height", "loop height (nm)"),
    "z_lys_unreacted_nm": ("height", "unreacted lysine height (nm)"),
    "z_rgd_nm": ("height", "RGD height (nm)"),
    "rgd_frac_above_h50": ("height", "RGD above h50 (fraction)"),
    "rg_nm": ("size", "Rg (nm)"),
    "ree_nm": ("size", "end-to-end distance (nm)"),
}


def boot_median_ci(x, rng):
    b = rng.choice(x, size=(NBOOT, len(x)), replace=True)
    return np.percentile(np.median(b, 1), [2.5, 97.5])


def cliffs_delta(x, y):
    d = np.sign(np.asarray(x)[:, None] - np.asarray(y)[None, :])
    return d.mean()


def boot_delta_ci(x, y, rng):
    bx = rng.choice(x, size=(NBOOT, len(x)), replace=True)
    by = rng.choice(y, size=(NBOOT, len(y)), replace=True)
    d = np.sign(bx[:, :, None] - by[:, None, :]).mean((1, 2))
    return np.percentile(d, [2.5, 97.5])


def compare(m):
    rng = np.random.default_rng(SEED)
    summ, tests, pw = [], [], []
    for r, (grp, lab) in RESPONSES.items():
        groups = {c: m.loc[m.design_class == c, r].dropna().to_numpy(float) for c in CLASSES}
        for c, x in groups.items():
            lo, hi = boot_median_ci(x, rng) if len(x) > 1 else (np.nan, np.nan)
            summ.append({"response": r, "group": grp, "label": lab, "class": c, "n": len(x), "median": np.median(x),
                         "median_ci_low": lo, "median_ci_high": hi, "q25": np.percentile(x, 25),
                         "q75": np.percentile(x, 75), "mean": x.mean(), "sd": x.std(ddof=1)})
        H, p = stats.kruskal(*[x for x in groups.values() if len(x)])
        tests.append({"response": r, "group": grp, "label": lab, "kruskal_H": H, "p": p})
        rows = []
        for a, b in itertools.combinations(CLASSES, 2):
            x, y = groups[a], groups[b]
            U, pu = stats.mannwhitneyu(x, y, alternative="two-sided")
            lo, hi = boot_delta_ci(x, y, rng)
            rows.append({"response": r, "class_a": a, "class_b": b, "median_a": np.median(x), "median_b": np.median(y),
                         "U": U, "p": pu, "cliffs_delta": cliffs_delta(x, y), "delta_ci_low": lo, "delta_ci_high": hi})
        padj = multipletests([x["p"] for x in rows], method="holm")[1]
        for x, q in zip(rows, padj):
            x["p_holm"] = q
        pw += rows
    return pd.DataFrame(summ), pd.DataFrame(tests), pd.DataFrame(pw)


def glm(df, rhs):
    d = df.assign(prop=df.n_bridges / MAXB, trials=MAXB)
    return smf.glm(f"prop ~ {rhs}", data=d, family=sm.families.Binomial(), var_weights=d.trials).fit(scale="X2")


def descriptor_models(m):
    from sklearn.model_selection import KFold
    m = m.assign(log_min_sep=np.log(m.min_sep_res))
    rhs = "log_min_sep + n_blocks_with_K + centroid_rel"
    rows, txt = [], []
    for label, d in [("all runs", m), ("random only", m[m.design_class == "random"])]:
        d = d.reset_index(drop=True)
        res = glm(d, rhs)
        pred = np.full(len(d), np.nan)
        for tr, te in KFold(5, shuffle=True, random_state=SEED).split(d):
            pred[te] = glm(d.iloc[tr], rhs).predict(d.iloc[te]) * MAXB
        cv_r2 = 1 - ((d.n_bridges - pred) ** 2).sum() / ((d.n_bridges - d.n_bridges.mean()) ** 2).sum()
        ci = res.conf_int()
        for t in res.params.index:
            rows.append({"data": label, "n": len(d), "term": t, "coef_logit": res.params[t], "se": res.bse[t],
                         "p": res.pvalues[t], "ci_low": ci.loc[t, 0], "ci_high": ci.loc[t, 1],
                         "odds_ratio": np.exp(res.params[t]), "pseudo_R2": 1 - res.deviance / res.null_deviance,
                         "cv5_R2_counts": cv_r2, "dispersion": res.scale})
        txt.append(f"{label}: n = {len(d)}, deviance pseudo-R2 {1 - res.deviance / res.null_deviance:.3f}, "
                   f"5-fold CV R2 (bridge counts) {cv_r2:.3f}, dispersion {res.scale:.2f}")
    # does the class label add anything once the descriptors are in? (quasi-F on deviance)
    full, red = glm(m, rhs + " + C(design_class)"), glm(m, rhs)
    q = int(red.df_resid - full.df_resid)
    F = (red.deviance - full.deviance) / q / full.scale
    p = float(stats.f.sf(F, q, full.df_resid))
    txt.append(f"class label beyond the descriptors (all runs): F({q}, {int(full.df_resid)}) = {F:.2f}, p = {p:.3g}")
    rows.append({"data": "all runs", "term": "C(design_class) | descriptors", "F": F, "p": p, "n": len(m)})
    return pd.DataFrame(rows), txt


def optimal(m, pw):
    lines = ["Optimality rule (settings.yaml):", S["optimal"]["rule"], ""]
    g = m.groupby("design_class")
    tab = pd.DataFrame({"median_bridges": g.n_bridges.median(),
                        "iqr_bridges": g.n_bridges.quantile(0.75) - g.n_bridges.quantile(0.25),
                        "frac_connected": g.largest_cluster_frac.apply(lambda s: (s >= 1).mean()),
                        "median_pairs_bridged": g.n_pairs_bridged.median(), "n": g.size()}).reindex(CLASSES)
    tab["viable"] = tab.frac_connected >= 0.9
    lines += ["Per class:", tab.round(3).to_string(), ""]
    viable = tab[tab.viable].sort_values("median_bridges", ascending=False)
    if viable.empty:
        lines.append("No class has a connected network in >= 90 % of its runs.")
        return "\n".join(lines), None
    top = viable.index[0]
    p = pw[pw.response == PRIMARY].set_index(["class_a", "class_b"]).p_holm

    def ph(a, b):
        return p.get((a, b), p.get((b, a)))

    tied = [top] + [c for c in viable.index[1:] if ph(top, c) >= ALPHA]
    lines.append(f"Step 1-2: viable classes by median bridges: {list(viable.index)}; top = {top}.")
    lines.append("Step 3: not separable from the top on n_bridges (Holm p >= 0.05): " + ", ".join(
        f"{c} (p = {ph(top, c):.3g})" for c in tied[1:]) if len(tied) > 1 else "Step 3: no class ties with the top.")
    best = tab.loc[tied].sort_values(["median_pairs_bridged", "iqr_bridges"], ascending=[False, True])
    choice = best.index[0]
    if len(tied) > 1:
        lines.append(f"Tie-break (distinct pairs bridged, then IQR): {best[['median_pairs_bridged', 'iqr_bridges']].to_dict('index')}")
    lines += ["", f"OPTIMAL CLASS: {choice}" + (f" (statistically tied with {', '.join(c for c in tied if c != choice)})"
                                                 if len(tied) > 1 else "")]
    return "\n".join(lines), choice


def main():
    STATS.mkdir(exist_ok=True)
    m = pd.read_csv(MOD / "runs_metrics.csv")
    summ, tests, pw = compare(m)
    summ.to_csv(STATS / "class_summary.csv", index=False)
    tests.to_csv(STATS / "class_tests.csv", index=False)
    pw.to_csv(STATS / "pairwise.csv", index=False)
    dg, dtxt = descriptor_models(m)
    dg.to_csv(STATS / "descriptor_glm.csv", index=False)
    otxt, _ = optimal(m, pw)
    (STATS / "optimal.txt").write_text(otxt + "\n")

    med = summ.pivot(index="response", columns="class", values="median").reindex(list(RESPONSES))[CLASSES]
    sc = pw[(pw.class_a == "uniform") & (pw.class_b == "scramble")].set_index("response")
    txt = ["RGD+ lysine placement - class comparison (n = 130; medians per class)", "",
           med.round(3).to_string(), "", "Kruskal-Wallis:", tests.set_index("response")[["kruskal_H", "p"]]
           .round(4).to_string(), "",
           f"Primary ({PRIMARY}) pairwise, Holm-adjusted:",
           pw[pw.response == PRIMARY][["class_a", "class_b", "median_a", "median_b", "p_holm", "cliffs_delta",
                                       "delta_ci_low", "delta_ci_high"]].round(4).to_string(index=False), "",
           "uniform vs scramble (Holm p, Cliff's delta) per response:",
           sc[["p_holm", "cliffs_delta"]].round(3).to_string(), "",
           "Descriptor GLM (quasi-binomial, n_bridges / 22):", *dtxt,
           dg.drop(columns=[c for c in ("F",) if c in dg]).round(4).to_string(index=False), "", otxt]
    (STATS / "summary.txt").write_text("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
