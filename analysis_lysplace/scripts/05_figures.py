"""Step 5 - figures (Plotly HTML with hover + PNG), MOSAIC style from analysis/scripts/style.py.

Colour: the classes are ordinal in lysine spread (clustered < focused < random < uniform), so they
get one blue ramp, light to dark, plus a distinct symbol each. Figures are transparent (St.wiki_mode():
mid-grey text and axes, no white outlines) so they read on white and on near-black (#1a1a1a); no ramp
of 4 ordinal steps clears 3:1 on both, so the ramp is kept in the middle (>= 2.2:1 on both) and the
symbol is the second cue. Classes are always also named on an axis or by a direct label. No figure
titles (they go in the caption). Within the z-profile panels, series differ by line style as well as colour (orange and
green are close for deuteranopes).

    analysis/scripts/run.sh analysis_lysplace/scripts/05_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(1, str(Path(__file__).resolve().parents[2] / "analysis" / "scripts"))
import style as St  # noqa: E402
from common import CLASSES, FIGURES, MOD, STATS  # noqa: E402

St.wiki_mode()
RAMP = {"clustered": "#7C9CF5", "focused": "#4F7FF5", "random": "#2A5BE0", "uniform": "#1F48B8"}
SYM = {"clustered": "circle", "focused": "square", "random": "triangle-up", "uniform": "diamond"}
RNG = np.random.default_rng(20261016)


def strip(fig, m, col, row, c, ylab):
    for k in CLASSES:
        y = m.loc[m.design_class == k, col].to_numpy(float)
        x = np.full(len(y), CLASSES.index(k)) + RNG.uniform(-0.18, 0.18, len(y))
        fig.add_trace(go.Box(x=np.full(len(y), CLASSES.index(k)), y=y, boxpoints=False, fillcolor="rgba(0,0,0,0)",
                             line=dict(color=St.MUTED, width=1.5), width=0.6, showlegend=False, hoverinfo="skip"),
                      row, c)
        fig.add_trace(go.Scatter(x=x, y=y, mode="markers", showlegend=False,
                                 marker=dict(color=RAMP[k], symbol=SYM[k], size=8, line=dict(color=St.EDGE, width=1)),
                                 text=m.loc[m.design_class == k, "Name"],
                                 hovertemplate=f"%{{text}}<br>{ylab}: %{{y:.3g}}<extra>{k}</extra>"), row, c)
    fig.update_xaxes(tickvals=list(range(len(CLASSES))), ticktext=CLASSES, row=row, col=c, showgrid=False)
    fig.update_yaxes(title_text=ylab, row=row, col=c)


def class_panels(m):
    panels = [("n_bridges", "inter-chain crosslinks (bridges)"), ("n_loops", "intra-chain crosslinks (loops)"),
              ("n_pairs_bridged", "distinct chain pairs bridged"), ("largest_cluster_frac", "largest connected cluster / 16"),
              ("z_bridge_nm", "bridge height above wall (nm)"), ("bridge_frac_above_h50", "bridges above h50 (fraction)")]
    fig = make_subplots(rows=3, cols=2, horizontal_spacing=0.1, vertical_spacing=0.09)
    for i, (col, lab) in enumerate(panels):
        strip(fig, m, col, i // 2 + 1, i % 2 + 1, lab)
    St.layout(fig, None, width=1100, height=1150)
    St.save(fig, FIGURES / "class_comparison")


def z_profiles(m):
    p = np.load(MOD / "profiles.npz", allow_pickle=True)
    z = 0.5 * (p["edges_A"][1:] + p["edges_A"][:-1]) / 10
    keep = z <= 25
    slugs = list(p["slugs"])
    fig = make_subplots(rows=1, cols=len(CLASSES), shared_yaxes=True, subplot_titles=CLASSES, horizontal_spacing=0.025)
    series = [("prof_mass", "all beads (mass)", St.MUTED, "solid", True), ("prof_bridge", "bridges", St.BLUE, "solid", False),
              ("prof_loop", "loops", St.ORANGE, "dash", False), ("prof_rgd", "RGD motifs", St.GREEN, "dot", False)]
    for j, k in enumerate(CLASSES, start=1):
        idx = [slugs.index(s) for s in m.loc[m.design_class == k, "slug"]]
        for key, lab, color, dash, fill in series:
            y = p[key][idx].sum(0)
            y = y / (y.sum() * (z[1] - z[0])) if y.sum() else y
            fig.add_trace(go.Scatter(x=y[keep], y=z[keep], mode="lines", name=lab, legendgroup=lab, showlegend=j == 1,
                                     line=dict(color=color, width=2, dash=dash),
                                     fill="tozerox" if fill else None,
                                     fillcolor="rgba(110,110,115,0.15)" if fill else None,
                                     hovertemplate=f"{lab}<br>z %{{y:.1f}} nm: %{{x:.3f}} /nm<extra>{k}</extra>"), 1, j)
        h95 = m.loc[m.design_class == k, "h95_nm"].median()
        fig.add_hline(y=h95, line=dict(color=St.MUTED, width=1, dash="dot"), row=1, col=j)
        fig.update_xaxes(title_text="density (1/nm)", row=1, col=j)
    fig.update_yaxes(title_text="height above wall (nm)", row=1, col=1)
    St.layout(fig, None, width=1150, height=520)
    fig.update_layout(legend=dict(orientation="h", y=1.14, x=0))
    St.save(fig, FIGURES / "z_profiles")


def lysine_curves():
    L = pd.read_csv(MOD / "per_lysine.csv")
    free = L[~L.pinned]
    cb = pd.read_csv(STATS / "lysine_gee_bridge.csv").set_index("term").coef_logit
    cl = pd.read_csv(STATS / "lysine_gee_loop.csv").set_index("term").coef_logit
    bins = np.array([4, 7, 12, 17, 27, 45, 70, 110, 160, 250])
    free = free.assign(b=pd.cut(free.sep_nearest_res, bins))
    g = free.groupby("b", observed=True)
    xs = np.exp(np.linspace(np.log(5), np.log(240), 100))
    fig = go.Figure()
    for st, coef, color, sym, dash in [("bridge", cb, St.BLUE, "circle", "solid"), ("loop", cl, St.ORANGE, "square", "dash")]:
        k = g.state.apply(lambda s: (s == st).sum())
        n = g.size()
        ph = k / n
        z = 1.96
        lo = (ph + z**2 / (2 * n) - z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)
        hi = (ph + z**2 / (2 * n) + z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)
        xm = g.sep_nearest_res.median()
        fig.add_trace(go.Scatter(x=xm, y=ph, mode="markers", name=f"{st}, observed (binned)",
                                 marker=dict(color=color, symbol=sym, size=9, line=dict(color=St.EDGE, width=1.5)),
                                 error_y=dict(type="data", symmetric=False, array=hi - ph, arrayminus=ph - lo,
                                              color=color, thickness=1.5, width=0),
                                 customdata=n, hovertemplate="sep %{x:.0f} res: %{y:.2f} (n = %{customdata})<extra>" + st + "</extra>"))
        eta = coef["Intercept"] + coef["np.log(sep_nearest_res)"] * np.log(xs) + coef["rel_pos"] * 0.5 + coef["I(rel_pos ** 2)"] * 0.25
        fig.add_trace(go.Scatter(x=xs, y=expit(eta), mode="lines", name=f"{st}, model (nearest K free, mid-chain)",
                                 line=dict(color=color, width=2, dash=dash), hoverinfo="skip"))
    St.layout(fig, None, width=820, height=560)
    fig.update_xaxes(type="log", tickvals=[5, 10, 20, 50, 100, 200], ticktext=["5", "10", "20", "50", "100", "200"],
                     title_text="sequence distance to the nearest other lysine on the same chain (residues)")
    fig.update_yaxes(title_text="probability", range=[0, 1])
    fig.update_layout(legend=dict(x=0.02, y=0.98))
    St.save(fig, FIGURES / "lysine_bridge_probability")


def model_check(m):
    cv = pd.read_csv(STATS / "lysine_model_cv.csv")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 22], y=[0, 22], mode="lines", line=dict(color=St.MUTED, dash="dot", width=1),
                             showlegend=False, hoverinfo="skip"))
    for k in CLASSES:
        d = cv[cv.design_class == k]
        fig.add_trace(go.Scatter(x=d.pred_bridges_cv, y=d.n_bridges + RNG.uniform(-0.15, 0.15, len(d)), mode="markers",
                                 name=k, marker=dict(color=RAMP[k], symbol=SYM[k], size=9, line=dict(color=St.EDGE, width=1)),
                                 text=d.Name, hovertemplate="%{text}<br>predicted %{x:.1f}, observed %{y:.0f}<extra></extra>"))
    St.layout(fig, None, width=640, height=600)
    fig.update_xaxes(title_text="predicted bridges (sum of P(bridge) / 2)", range=[0, 22])
    fig.update_yaxes(title_text="observed bridges", range=[0, 22], scaleanchor="x")
    fig.update_layout(legend=dict(x=0.02, y=0.98))
    St.save(fig, FIGURES / "lysine_model_check")


def design_rule(m):
    r = pd.read_csv(STATS / "design_rule.csv")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=r.min_sep_res, y=r["max"], mode="lines", name="best design at this spacing (model)",
                             line=dict(color=St.INK, width=2), hovertemplate="min sep %{x}: best %{y:.1f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=r.min_sep_res, y=r["median"], mode="lines", name="median design at this spacing (model)",
                             line=dict(color=St.MUTED, width=2, dash="dash"),
                             hovertemplate="min sep %{x}: median %{y:.1f}<extra></extra>"))
    for k in CLASSES:
        d = m[m.design_class == k]
        fig.add_trace(go.Scatter(x=d.min_sep_res + RNG.uniform(-1.2, 1.2, len(d)), y=d.n_bridges, mode="markers", name=f"{k} (simulated)",
                                 marker=dict(color=RAMP[k], symbol=SYM[k], size=8, line=dict(color=St.EDGE, width=1)),
                                 text=d.Name, hovertemplate="%{text}<br>min sep %{x:.0f}, bridges %{y}<extra></extra>"))
    fig.add_hline(y=22, line=dict(color=St.MUTED, width=1, dash="dot"))
    St.layout(fig, None, width=900, height=560)
    fig.update_xaxes(title_text="minimum sequence distance between lysines on a chain (residues)")
    fig.update_yaxes(title_text="inter-chain crosslinks per 16 chains", range=[0, 23])
    fig.update_layout(legend=dict(x=0.6, y=0.32))
    St.save(fig, FIGURES / "design_rule")


def timecourse(m):
    B = pd.read_csv(MOD / "per_bond.csv")
    t = np.linspace(0, 390, 400)
    fig = go.Figure()
    for k in CLASSES:
        curves = []
        for s in m.loc[m.design_class == k, "slug"]:
            tb = np.sort(B.loc[(B.slug == s) & (B.kind == "bridge"), "t_ns"].to_numpy())
            curves.append(np.searchsorted(tb, t, side="right"))
        med = np.median(curves, 0)
        fig.add_trace(go.Scatter(x=t, y=med, mode="lines", name=k,
                                 line=dict(color=RAMP[k], width=2),
                                 hovertemplate="%{x:.0f} ns: %{y:.1f} bridges<extra>" + k + "</extra>"))
        fig.add_annotation(x=395, y=med[-1], text=k, showarrow=False, xanchor="left", font=dict(size=12, color=St.INK))
    St.layout(fig, None, width=860, height=520)
    fig.update_xaxes(title_text="time after crosslinking starts (ns)", range=[0, 450])
    fig.update_yaxes(title_text="inter-chain crosslinks")
    fig.update_layout(showlegend=False)
    St.save(fig, FIGURES / "bridge_timecourse")


if __name__ == "__main__":
    FIGURES.mkdir(exist_ok=True)
    m = pd.read_csv(MOD / "runs_metrics.csv")
    for f in (class_panels, z_profiles, model_check, design_rule, timecourse):
        f(m)
        print("done", f.__name__)
    lysine_curves()
    print("done lysine_curves")
