"""Step 6 - draft figures (Plotly, Poppins, MOSAIC palette). HTML + PNG in analysis/figures/.

    analysis/scripts/run.sh analysis/scripts/05_figures.py
    # the same figures for another model folder written with the stats/ schema, e.g. log coding:
    analysis/scripts/run.sh analysis/scripts/05_figures.py --stats DIR --out DIR --coding log
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import model as M  # noqa: E402
import style as St  # noqa: E402
from common import ANALYSIS, BASE, FIGURES, RANGES  # noqa: E402

STATS = ANALYSIS / "stats"
CODING = "linear"          # design-variable coding of the model in STATS (frozen: linear)
NOTE = ""                  # appended to figure titles for a non-frozen model
YLAB = "inter-chain K–K contacts per frame"
TERM_LABEL = {"C": "C", "L": "L", "S": "S", "ion": "I"}


def nice_term(t):
    if t.startswith("I("):
        v = t[2:].split("**")[0].strip()
        return f"{TERM_LABEL[v]}²"
    return " × ".join(TERM_LABEL[v] for v in t.split(":"))


def coefficient_plot():
    c = pd.read_csv(STATS / "coefficients.csv")
    c = c[c.term != "Intercept"].copy()
    c["label"] = c.term.map(nice_term)
    c = c.reindex(c.coef_log.abs().sort_values().index)
    sig = (c.ci_low > 0) | (c.ci_high < 0)
    fig = go.Figure()
    for mask, color, name in [(sig, St.BLUE, "95 % CI excludes 0"), (~sig, St.MUTED, "95 % CI includes 0")]:
        d = c[mask]
        fig.add_trace(go.Scatter(
            x=d.coef_log, y=d.label, mode="markers", name=name,
            marker=dict(color=color, size=10, line=dict(color="white", width=2)),
            error_x=dict(type="data", symmetric=False, array=d.ci_high - d.coef_log,
                         arrayminus=d.coef_log - d.ci_low, color=color, thickness=2, width=0),
            customdata=np.c_[d.ci_low, d.ci_high, d.p],
            hovertemplate="%{y}: %{x:.3f} [%{customdata[0]:.3f}, %{customdata[1]:.3f}], p = %{customdata[2]:.2g}<extra></extra>"))
    fig.add_vline(x=0, line=dict(color=St.MUTED, width=1))
    St.layout(fig, "Standardised coefficients, log response" + NOTE, width=760, height=600)
    fig.update_xaxes(title_text="coefficient (change in log response per coded unit)")
    fig.update_yaxes(categoryorder="array", categoryarray=list(c.label))
    fig.update_layout(legend=dict(orientation="h", y=1.06, x=0))
    St.save(fig, FIGURES / "coefficients")


def predicted_vs_heldout():
    d = pd.read_csv(STATS / "run_diagnostics.csv")
    fit = pd.read_csv(STATS / "model_fit.csv").iloc[0]
    fig = go.Figure()
    lo, hi = d.y.min() * 0.8, d.y.max() * 1.25
    fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=St.MUTED, width=1, dash="dot"),
                             showlegend=False, hoverinfo="skip"))
    for mask, color, sym, name in [(d.anchor_id.isna(), St.BLUE, "circle", "single run"),
                                   (d.anchor_id.notna(), St.ORANGE, "square", "anchor replicate")]:
        e = d[mask]
        fig.add_trace(go.Scatter(x=np.exp(e.cv_pred_ols), y=e.y, mode="markers", name=name,
                                 marker=dict(color=color, symbol=sym, size=9, line=dict(color="white", width=1.5)),
                                 text=e.slug, hovertemplate="%{text}<br>held-out %{x:.3f}, observed %{y:.3f}<extra></extra>"))
    St.layout(fig, f"Held-out prediction (5-fold grouped CV, R² = {fit.cv5_grouped_R2_ols:.2f}){NOTE}", width=640, height=600)
    fig.update_xaxes(type="log", title_text=f"predicted, held out ({YLAB})", range=np.log10([lo, hi]))
    fig.update_yaxes(type="log", title_text=f"observed ({YLAB})", range=np.log10([lo, hi]), scaleanchor="x")
    fig.update_layout(legend=dict(x=0.02, y=0.98))
    St.save(fig, FIGURES / "predicted_vs_heldout")


def _surface_traces(res, px, py, fixed, zmin, zmax, show_scale, n_arrows=8, arrow_frac=0.06):
    xs, ys, logy = M.predict_grid(res, px, py, fixed, coding=CODING)
    z = np.exp(logy)
    cont = go.Contour(x=xs, y=ys, z=z, colorscale=St.SURFACE_SCALE, zmin=zmin, zmax=zmax,
                      contours=dict(showlabels=True, labelfont=dict(family=St.FONT, size=11, color=St.INK),
                                    start=0.1, end=1.5, size=0.1),
                      line=dict(color="rgba(29,29,31,0.45)", width=1), showscale=show_scale,
                      colorbar=dict(title=dict(text=YLAB, side="right"), thickness=14),
                      hovertemplate=f"{M.LABEL[px]} %{{x:.3g}}<br>{M.LABEL[py]} %{{y:.3g}}<br>predicted %{{z:.3f}}<extra></extra>")
    # gradient of log(response) in coded units, drawn in natural units of the axes
    fixed_full = {**BASE, **fixed}
    ax_pts = np.linspace(*RANGES[px], n_arrows + 2)[1:-1]
    ay_pts = np.linspace(*RANGES[py], n_arrows + 2)[1:-1]
    grads = []
    for gx in ax_pts:
        for gy in ay_pts:
            pt = {M.SHORT[p]: float(M.to_coded(p, fixed_full[p], CODING)) for p in M.PARAMS}
            pt[M.SHORT[px]] = float(M.to_coded(px, gx, CODING))
            pt[M.SHORT[py]] = float(M.to_coded(py, gy, CODING))
            g = M.gradient_coded(res, pt)
            grads.append((gx, gy, g[M.SHORT[px]], g[M.SHORT[py]]))
    # Fixed-length arrows show the direction of steepest ascent (in coded units); how steep
    # it is shows in the contour spacing. Arrows where the surface is essentially flat
    # (< 5 % of the largest gradient on the panel) are left out rather than drawn as stubs.
    gmax = max(np.hypot(a, b) for *_, a, b in grads) or 1.0
    wx, wy = np.ptp(RANGES[px]), np.ptp(RANGES[py])
    arrows = []
    for gx, gy, a, b in grads:
        norm = np.hypot(a, b)
        if norm < 0.05 * gmax:
            continue
        # coded-unit direction -> fraction of each axis (the Jacobian is constant for linear coding)
        fx = a * float(M.natural_per_coded(px, gx, CODING)) / wx
        fy = b * float(M.natural_per_coded(py, gy, CODING)) / wy
        nf = np.hypot(fx, fy)
        dx, dy = fx / nf * arrow_frac * wx, fy / nf * arrow_frac * wy
        arrows.append(dict(x=gx + dx / 2, y=gy + dy / 2, ax=gx - dx / 2, ay=gy - dy / 2,
                           showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1.6, arrowcolor=St.INK))
    return cont, arrows


def response_surfaces():
    res = M.fit(M.load(coding=CODING))
    pt = pd.read_csv(STATS / "parameter_tests.csv").sort_values("lmg_share_R2", ascending=False)
    px, py = pt.parameter.iloc[0], pt.parameter.iloc[1]
    others = [p for p in M.PARAMS if p not in (px, py)]
    zmin, zmax = 0.0, 1.4

    cont, arrows = _surface_traces(res, px, py, {}, zmin, zmax, True)
    fig = go.Figure(cont)
    for a in arrows:
        a.update(xref="x", yref="y", axref="x", ayref="y")
    St.layout(fig, "Predicted stiffness proxy, other parameters at base case" + NOTE, width=820, height=640)
    fig.update_layout(annotations=arrows)
    fig.update_xaxes(title_text=M.LABEL[px], showgrid=False)
    fig.update_yaxes(title_text=M.LABEL[py], showgrid=False)
    St.save(fig, FIGURES / "response_surface_base")

    # 2x2: the other two parameters at low / high (coded -0.8 / +0.8, i.e. 10th / 90th % of range)
    lvl = {p: {"low": RANGES[p][0] + 0.1 * np.ptp(RANGES[p]), "high": RANGES[p][0] + 0.9 * np.ptp(RANGES[p])}
           for p in others}
    combos = [(a, b) for b in ("high", "low") for a in ("low", "high")]
    fmt = {"chain_length": lambda v: f"L = {v:.0f}", "ionic_strength": lambda v: f"I = {v:.3g} M",
           "concentration": lambda v: f"C = {v:.3g}", "spacing_pct": lambda v: f"S = {v:.3g} %"}
    titles = [f"{fmt[others[0]](lvl[others[0]][a])}, {fmt[others[1]](lvl[others[1]][b])}" for a, b in combos]
    fig = make_subplots(rows=2, cols=2, subplot_titles=titles, shared_xaxes=True, shared_yaxes=True,
                        horizontal_spacing=0.06, vertical_spacing=0.1)
    all_arrows = []
    for i, (a, b) in enumerate(combos):
        r, c = i // 2 + 1, i % 2 + 1
        fixed = {others[0]: lvl[others[0]][a], others[1]: lvl[others[1]][b]}
        cont, arrows = _surface_traces(res, px, py, fixed, zmin, zmax, i == 0, n_arrows=6)
        fig.add_trace(cont, r, c)
        k = "" if i == 0 else str(i + 1)
        for ar in arrows:
            ar.update(xref=f"x{k}", yref=f"y{k}", axref=f"x{k}", ayref=f"y{k}")
        all_arrows += arrows
    St.layout(fig, None, width=1100, height=900)
    fig.update_layout(annotations=list(fig.layout.annotations) + all_arrows)
    fig.update_annotations(selector=dict(showarrow=False), font=dict(family=St.FONT, size=13, color=St.INK))
    fig.update_xaxes(title_text=M.LABEL[px], row=2, showgrid=False)
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(title_text=M.LABEL[py], col=1, showgrid=False)
    fig.update_yaxes(showgrid=False)
    St.save(fig, FIGURES / "response_surface_2x2")


def residual_diagnostics():
    d = pd.read_csv(STATS / "run_diagnostics.csv")
    osm, osr = stats.probplot(d.resid_studentized, dist="norm", fit=False)
    order = np.argsort(d.resid_studentized.to_numpy())
    fig = make_subplots(rows=1, cols=3, subplot_titles=["residuals vs fitted", "normal Q–Q", "Cook's distance"],
                        horizontal_spacing=0.08)
    fig.add_trace(go.Scatter(x=d.fitted, y=d.resid_studentized, mode="markers", text=d.slug,
                             marker=dict(color=St.BLUE, size=8, line=dict(color="white", width=1)),
                             hovertemplate="%{text}<br>fitted %{x:.2f}, stud. resid %{y:.2f}<extra></extra>",
                             showlegend=False), 1, 1)
    fig.add_hline(y=0, line=dict(color=St.MUTED, width=1), row=1, col=1)
    fig.add_trace(go.Scatter(x=osm, y=osr, mode="markers", text=d.slug.to_numpy()[order],
                             marker=dict(color=St.BLUE, size=8, line=dict(color="white", width=1)),
                             hovertemplate="%{text}<extra></extra>", showlegend=False), 1, 2)
    fig.add_trace(go.Scatter(x=[-3, 3], y=[-3, 3], mode="lines", line=dict(color=St.MUTED, dash="dot", width=1),
                             showlegend=False, hoverinfo="skip"), 1, 2)
    thr = 4 / len(d)
    fig.add_trace(go.Bar(x=np.arange(len(d)), y=d.cooks_d, text=d.slug, showlegend=False,
                         marker=dict(color=np.where(d.cooks_d > thr, St.ORANGE, St.BLUE)),
                         hovertemplate="%{text}<br>Cook's D %{y:.3f}<extra></extra>"), 1, 3)
    fig.add_hline(y=thr, line=dict(color=St.ORANGE, width=1, dash="dash"), row=1, col=3)
    St.layout(fig, None, width=1300, height=430)
    fig.update_xaxes(title_text="fitted log response", row=1, col=1)
    fig.update_yaxes(title_text="studentised residual", row=1, col=1)
    fig.update_xaxes(title_text="theoretical quantile", row=1, col=2)
    fig.update_yaxes(title_text="studentised residual", row=1, col=2)
    fig.update_xaxes(title_text="run index", row=1, col=3)
    fig.update_yaxes(title_text="Cook's D", row=1, col=3)
    St.save(fig, FIGURES / "residual_diagnostics")


def noise_check():
    m = pd.read_csv(ANALYSIS / "runs_metrics.csv")
    a = m[m.anchor_id.notna()].groupby("anchor_id").inter_free_12A.agg(["mean", "std", "count"]).reset_index()
    a["cv"] = a["std"] / a["mean"]
    fig = make_subplots(rows=1, cols=2, subplot_titles=["replicate SD vs mean", "replicate CV vs mean"],
                        horizontal_spacing=0.12)
    for col, y, yt in [(1, a["std"], "replicate SD (contacts per frame)"), (2, a.cv, "coefficient of variation")]:
        fig.add_trace(go.Scatter(x=a["mean"], y=y, mode="markers+text", text=a.anchor_id, textposition="top center",
                                 textfont=dict(size=11, color=St.MUTED),
                                 marker=dict(color=St.ORANGE, symbol="square", size=10, line=dict(color="white", width=1.5)),
                                 hovertemplate="%{text}: mean %{x:.3f}, %{y:.3f}<extra></extra>", showlegend=False), 1, col)
        fig.update_yaxes(title_text=yt, row=1, col=col, rangemode="tozero")
        fig.update_xaxes(title_text=f"anchor mean ({YLAB})", row=1, col=col)
    St.layout(fig, None, width=1100, height=440)
    St.save(fig, FIGURES / "noise_check_anchors")


def importance_comparison():
    pt = pd.read_csv(STATS / "parameter_tests.csv").set_index("parameter")
    imp = pd.read_csv(STATS / "importance_rf.csv").set_index("parameter")
    share = pd.DataFrame({
        "quadratic model, LMG share of R²": pt.lmg_share_of_explained,
        "random forest, held-out permutation": imp.perm_importance_heldout / imp.perm_importance_heldout.sum(),
        "random forest, mean |SHAP|": imp.mean_abs_shap / imp.mean_abs_shap.sum()})
    share = share.loc[share.iloc[:, 0].sort_values().index]
    fig = go.Figure()
    for i, col in enumerate(share.columns):
        fig.add_trace(go.Bar(y=[M.LABEL[p] for p in share.index], x=share[col], name=col, orientation="h",
                             marker=dict(color=[St.BLUE, St.PURPLE, St.GREEN][i]),
                             hovertemplate="%{y}: %{x:.1%}<extra>" + col + "</extra>"))
    St.layout(fig, "Parameter importance, normalised to sum to 1" + NOTE, width=900, height=460)
    fig.update_layout(barmode="group", bargap=0.25, bargroupgap=0.08, legend=dict(orientation="h", y=-0.2, x=0))
    fig.update_xaxes(tickformat=".0%", title_text="share")
    St.save(fig, FIGURES / "importance_comparison")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", type=Path, default=STATS, help="model folder with the stats/ file schema")
    ap.add_argument("--out", type=Path, default=FIGURES)
    ap.add_argument("--coding", choices=["linear", "log"], default="linear")
    a = ap.parse_args()
    STATS, FIGURES, CODING = a.stats, a.out, a.coding
    NOTE = "" if CODING == "linear" else " — log-coded parameters"
    for f in (coefficient_plot, predicted_vs_heldout, response_surfaces, residual_diagnostics,
              noise_check, importance_comparison):
        f()
        print("done", f.__name__)
