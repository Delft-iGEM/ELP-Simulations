"""Step 4 - what makes one lysine form a bridge, and which placement maximises bridges.

Unit: every lysine that is not surface-pinned (130 runs x 45). Logistic GEE (exchangeable within
run) for P(bridge) and P(loop), predictors fixed in settings.yaml:
  log(sep_nearest_res)  sequence distance to the nearest other K on the same chain
  n_k_within_30         other K within 30 residues on the same chain
  nearest_pinned        whether that nearest K is surface-pinned
  rel_pos, rel_pos^2    position along the chain (0 = N-terminus)

Validation: expected bridges per run = sum P(bridge) / 2 over its free lysines, 10-fold CV by run.
Optimisation: every design with 4 K on the 49 VPGXG positions that are not M (C(49, 4) = 211 876),
pinning averaged analytically (each K pinned with p = 19/64), expected bridges per 16 chains.
Recommendations are flagged when a predictor leaves the range seen in the data.

Writes analysis_lysplace/stats/lysine_gee_bridge.csv, lysine_gee_loop.csv, lysine_model_cv.csv,
design_ranking.csv (top 200 + every observed design), design_top.csv (with 95 % intervals),
design_rule.csv (best / median expected bridges per minimum K-K separation),
lysine_model.txt

    analysis/scripts/run.sh analysis_lysplace/scripts/04_lysine_model.py
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
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).parent))
from common import MOD, SIMS, STATS, parse_fasta  # noqa: E402

S = yaml.safe_load((MOD / "settings.yaml").read_text())
SEED = int(S["seed"])
Q_PIN = 19 / 64
N_CH = 16
FORMULA = "y ~ np.log(sep_nearest_res) + n_k_within_30 + nearest_pinned + rel_pos + I(rel_pos ** 2)"


def fit(df, outcome):
    d = df.assign(y=(df.state == outcome).astype(float), nearest_pinned=df.nearest_pinned.astype(float))
    return smf.gee(FORMULA, groups="slug", data=d, family=sm.families.Binomial(),
                   cov_struct=sm.cov_struct.Exchangeable()).fit()


def coef_table(res):
    ci = res.conf_int()
    return pd.DataFrame({"term": res.params.index, "coef_logit": res.params.values, "se": res.bse.values,
                         "p": res.pvalues.values, "ci_low": ci[0].values, "ci_high": ci[1].values,
                         "odds_ratio": np.exp(res.params.values)})


def design_features(k_res: np.ndarray, length: int):
    """Per-K predictors for designs given as (n_designs, 4) residue indices (sorted)."""
    d = np.abs(k_res[:, :, None] - k_res[:, None, :]).astype(float)
    idx = np.arange(4)
    d[:, idx, idx] = np.inf
    near = d.min(2)
    n_tie = np.isclose(d, near[:, :, None]).sum(2)                     # nearest neighbours at the same distance
    p_np = 1 - (1 - Q_PIN) ** n_tie                                    # P(at least one nearest K pinned)
    return near, (d < 30).sum(2), p_np, k_res / (length - 1)


def expected_bridges(params: np.ndarray, near, n30, p_np, rel):
    """Expected bridges per 16 chains, pinning averaged; params in FORMULA order."""
    b0, bs, b30, bnp, br, br2 = params
    eta = b0 + bs * np.log(near) + b30 * n30 + br * rel + br2 * rel ** 2
    p = (1 - p_np) * expit(eta) + p_np * expit(eta + bnp)
    return N_CH * ((1 - Q_PIN) * p).sum(-1) / 2


def main():
    L = pd.read_csv(MOD / "per_lysine.csv")
    m = pd.read_csv(MOD / "runs_metrics.csv")
    free = L[~L.pinned].reset_index(drop=True)
    rb, rl = fit(free, "bridge"), fit(free, "loop")
    cb, cl = coef_table(rb), coef_table(rl)
    cb.to_csv(STATS / "lysine_gee_bridge.csv", index=False)
    cl.to_csv(STATS / "lysine_gee_loop.csv", index=False)

    # run-level validation, 10-fold CV by run (GroupKFold has no shuffle; folds are deterministic)
    from sklearn.model_selection import GroupKFold
    pred = pd.Series(np.nan, index=free.index)
    for tr, te in GroupKFold(10).split(free, groups=free.slug):
        r = fit(free.iloc[tr], "bridge")
        pred.iloc[te] = r.predict(free.iloc[te].assign(nearest_pinned=free.nearest_pinned.iloc[te].astype(float)))
    run = (free.assign(p=pred).groupby("slug").p.sum() / 2).rename("pred_bridges_cv").reset_index()
    run = run.merge(m[["slug", "Name", "design_class", "n_bridges"]], on="slug")
    r2 = 1 - ((run.n_bridges - run.pred_bridges_cv) ** 2).sum() / ((run.n_bridges - run.n_bridges.mean()) ** 2).sum()
    rmse = float(np.sqrt(((run.n_bridges - run.pred_bridges_cv) ** 2).mean()))
    run.to_csv(STATS / "lysine_model_cv.csv", index=False)
    cls_cal = run.groupby("design_class")[["n_bridges", "pred_bridges_cv"]].mean()

    # design space: 4 K on the VPGXG positions of the observed construct that are not M
    fa = parse_fasta(SIMS / m.slug.iloc[0] / "runtime" / "molecules.fasta")
    slots = [(bi + 1, p, b["start"] + (p - 1) * 5 + 3) for bi, b in enumerate(fa["blocks"])
             for p in range(1, b["n_repeats"] + 1) if p not in b["M"]]
    combos = np.array(list(itertools.combinations(range(len(slots)), 4)))
    kres = np.array([s[2] for s in slots])[combos]
    near, n30, p_np, rel = design_features(kres, fa["length"])
    params = rb.params.to_numpy()
    eb = expected_bridges(params, near, n30, p_np, rel)
    obs_near = free.sep_nearest_res
    in_range = (near.min(1) >= obs_near.min()) & (near.max(1) <= obs_near.max())

    def label(c):
        """Design name in the run-name convention, e.g. '13R1G13D4' or '0R0G0D3+4+5+7'."""
        per = [[] for _ in range(4)]
        for i in c:
            per[slots[i][0] - 1].append(slots[i][1])
        b = ["+".join(map(str, x)) or "0" for x in per]
        return f"{b[0]}R{b[1]}G{b[2]}D{b[3]}"

    order = np.argsort(-eb)
    rows = pd.DataFrame({"design": [label(combos[i]) for i in order[:200]], "expected_bridges": eb[order[:200]],
                         "min_sep_res": near[order[:200]].min(1), "in_data_range": in_range[order[:200]],
                         "k_per_block": ["-".join(str(sum(slots[j][0] == b for j in combos[i])) for b in (1, 2, 3, 4))
                                         for i in order[:200]], "rank": np.arange(1, 201)})
    # observed designs, scored the same way
    obs = m[["Name", "design_class", "n_bridges", "K_residues"]].copy()
    ok = np.array([sorted(map(int, s.split("+"))) for s in obs.K_residues])
    f = design_features(ok, fa["length"])
    obs["expected_bridges"] = expected_bridges(params, *f)
    obs["rank_in_design_space"] = [int(np.searchsorted(-np.sort(eb)[::-1], -e, side="left")) + 1
                                   for e in obs.expected_bridges]
    pd.concat([rows.assign(source="design space"),
               obs.rename(columns={"Name": "design"}).assign(source="simulated")], ignore_index=True) \
        .to_csv(STATS / "design_ranking.csv", index=False)

    # uncertainty for the top designs and class means: parameter draws from the GEE covariance
    rng = np.random.default_rng(SEED)
    draws = rng.multivariate_normal(params, rb.cov_params().to_numpy(), size=2000)
    top = order[:10]
    sims = np.stack([expected_bridges(dr, near[top], n30[top], p_np[top], rel[top]) for dr in draws])
    best_any = np.stack([expected_bridges(dr, near, n30, p_np, rel).max() for dr in draws[:200]])
    dtop = pd.DataFrame({"design": [label(combos[i]) for i in top], "expected_bridges": eb[top],
                         "ci_low": np.percentile(sims, 2.5, 0), "ci_high": np.percentile(sims, 97.5, 0),
                         "min_sep_res": near[top].min(1), "in_data_range": in_range[top]})
    dtop.to_csv(STATS / "design_top.csv", index=False)

    # design rule: best achievable expected bridges by minimum K-K separation
    mins = near.min(1)
    rule = pd.DataFrame({"min_sep_res": mins, "eb": eb}).groupby("min_sep_res").eb.agg(["max", "median", "size"])
    rule.reset_index().to_csv(STATS / "design_rule.csv", index=False)

    txt = ["RGD+ lysine placement - per-lysine model (free lysines, n = %d in %d runs)" % (len(free), free.slug.nunique()),
           "", "P(bridge), logistic GEE (exchangeable within run):", cb.round(4).to_string(index=False), "",
           "P(loop):", cl.round(4).to_string(index=False), "",
           f"Run-level check, 10-fold CV by run: expected bridges = sum P / 2; R2 = {r2:.3f}, RMSE = {rmse:.2f} bridges",
           "class means, observed vs CV-predicted:", cls_cal.round(2).to_string(), "",
           f"Design space: {len(eb)} designs of 4 K on {len(slots)} positions; expected bridges per 16 chains "
           f"{eb.min():.1f}-{eb.max():.1f}; {in_range.mean():.1%} of designs inside the data range of sep_nearest.",
           f"Best design in every one of 200 parameter draws reaches {np.percentile(best_any, 2.5):.1f}-"
           f"{np.percentile(best_any, 97.5):.1f} expected bridges (95 %).", "",
           "Top 10 designs (95 % interval from the GEE covariance):", dtop.round(2).to_string(index=False), "",
           "Best and median expected bridges by minimum K-K separation (residues):",
           rule.loc[rule.index <= 120].round(2).to_string(), "",
           "Simulated designs: expected (model) vs observed bridges, best 10 observed:",
           obs.sort_values("n_bridges", ascending=False).head(10).round(2).to_string(index=False)]
    (STATS / "lysine_model.txt").write_text("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
