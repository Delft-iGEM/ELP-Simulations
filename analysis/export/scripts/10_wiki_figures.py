"""Step 10 - wiki figures: transparent background, log-coded model, Poppins, PNG (2x) + SVG.

Reads the exported models (analysis/export/models/*_log) and refits the same log-coded models
with the export's own code (07_models.py / model.py) - no change to data, design or model; the
refit is asserted to reproduce the exported coefficients. Writes analysis/export/figures_wiki/.

    analysis/scripts/run.sh analysis/scripts/10_wiki_figures.py
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).parent))
import model as M  # noqa: E402
import style as St  # noqa: E402
from common import ANALYSIS, BASE, RANGES, TIMESERIES  # noqa: E402

St.wiki_mode()                       # before anything draws
F5 = importlib.import_module("05_figures")
M7 = importlib.import_module("07_models")
ST4 = importlib.import_module("04_stats")

EXPORT = ANALYSIS / "export"
MODELS = EXPORT / "models"
OUT = EXPORT / "figures_wiki"
DATA = OUT / "data"

# ---- calibration for the tissue bands ---------------------------------------------------
G0_KPA = 3.0              # assumed storage modulus G'0 of the base construct
E_OVER_G = 3.0            # E = 3 G (incompressible)
# Target tissues, Young's modulus E in kPa (low, high).
# UNVERIFIED: placeholders written from memory of typical reported ranges, NOT taken from a
# source; each must be checked against the cited literature before the wiki uses it. No
# citation is given on purpose - add the reference next to each value once checked.
TISSUES = {
    "brain": (0.1, 1.0),     # UNVERIFIED, check against the cited literature (source: TODO)
    "fat": (1.0, 4.0),       # UNVERIFIED, check against the cited literature (source: TODO)
    "liver": (0.5, 2.0),     # UNVERIFIED, check against the cited literature (source: TODO)
    "muscle": (8.0, 17.0),   # UNVERIFIED, check against the cited literature (source: TODO)
}
NEAR_BASE = 0.5           # "within one design step of base": |coded log value - base| <= 0.5 on the other three
SEED_NOTE = "Illustrative: assumes G is proportional to the proxy and G'0 = 3 kPa. Needs one rheology measurement to calibrate."
EXCLUDE_RUN = "sub4-extra2-l58-s4-c0p01-i0p43-05"

TICKS = {"concentration": [0.01, 0.02, 0.04, 0.06], "chain_length": [50, 75, 100, 150],
         "spacing_pct": [5, 10, 15, 25], "ionic_strength": [0.05, 0.1, 0.2, 0.5]}
SHORTLAB = {"concentration": "grafting density", "chain_length": "chain length",
            "spacing_pct": "lysine spacing", "ionic_strength": "ionic strength"}
UNITS = {"concentration": "chains/nm²", "chain_length": "repeats", "spacing_pct": "% of length",
         "ionic_strength": "M"}


# ---- output -----------------------------------------------------------------------------
def _clean(fig):
    """Make a figure safe for both backgrounds and for SVG: WebGL traces -> SVG traces,
    white marker outlines and dark contour lines -> transparent / mid-grey."""
    if any(tr.type == "scattergl" for tr in fig.data):
        data = []
        for tr in fig.data:
            d = tr.to_plotly_json()
            data.append(go.Scatter(**{k: v for k, v in d.items() if k != "type"}) if tr.type == "scattergl" else tr)
        fig = go.Figure(data=data, layout=fig.layout)
    for tr in fig.data:
        if tr.type == "scatter" and tr.marker is not None and tr.marker.line is not None \
                and tr.marker.line.color == "white":
            tr.marker.line.color = St.CLEAR
        if tr.type in ("bar",) and tr.marker is not None and tr.marker.line is not None \
                and tr.marker.line.color == "white":
            tr.marker.line.color = St.CLEAR
        if tr.type == "contour" and tr.line is not None and tr.line.color and str(tr.line.color).startswith("rgba(29"):
            tr.line.color = "rgba(124,128,137,0.75)"
    return fig


def _arr(v):
    a = np.asarray(v)
    if a.dtype.kind == "f":
        return np.round(a, 6).tolist()
    return a.tolist()


def dump(fig, name, extra=None):
    """The plotted numbers, one compact JSON line per figure."""
    traces = []
    for tr in fig.data:
        d = {"name": tr.name, "type": tr.type}
        for k in ("x", "y", "z", "text", "base"):
            v = getattr(tr, k, None)
            if v is not None:
                d[k] = _arr(v)
        ey = getattr(tr, "error_y", None)
        ex = getattr(tr, "error_x", None)
        for lab, e in (("error_x", ex), ("error_y", ey)):
            if e is not None and e.array is not None:
                d[lab] = {"array": _arr(e.array), "arrayminus": _arr(e.arrayminus) if e.arrayminus is not None else None}
        traces.append(d)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / f"{name}.json").write_text(json.dumps({"figure": name, "traces": traces, **(extra or {})},
                                                  separators=(",", ":")))


def save(fig, stem, extra=None):
    name = Path(stem).name
    fig = _clean(fig)
    fig.update_layout(paper_bgcolor=St.CLEAR, plot_bgcolor=St.CLEAR)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.write_image(str(OUT / f"{name}.png"), scale=2)
    fig.write_image(str(OUT / f"{name}.svg"))
    dump(fig, name, extra)
    print("  wrote", name)


St.save = save                      # the 05_figures drawing functions call St.save(fig, FIGURES / name)
F5.STATS, F5.FIGURES, F5.CODING, F5.NOTE = MODELS / "inter_log", OUT, "log", ""


# ---- models (the export's own code; asserted equal to the exported coefficients) ----------
def fit_models():
    out = {}
    for resp in ("inter", "intra", "interfrac"):
        df = M7.load(resp, "log")
        res = M7.glm_fit(df) if resp == "interfrac" else M.fit(df)
        exp = pd.read_csv(MODELS / f"{resp}_log" / "coefficients.csv")
        col = "coef_logit" if resp == "interfrac" else "coef_log"
        assert np.allclose(res.params.values, exp[col].values, rtol=1e-8, atol=1e-10), f"{resp}_log not reproduced"
        out[resp] = (df, res)
    print("refit models reproduce models/{inter,intra,interfrac}_log/coefficients.csv")
    return out


MODELS_FIT = None


def coded_row(values, coding="log"):
    return {M.SHORT[p]: float(M.to_coded(p, values[p], coding)) for p in M.PARAMS}


def predict(resp, rows, coding="log"):
    """Mean prediction with 95 % CI (response scale) for rows of natural-unit parameter dicts.
    OLS: exp(eta +- t se); GLM: expit(eta +- z se)."""
    df, res = MODELS_FIT[resp]
    di = res.model.data.design_info
    X = np.asarray(patsy.build_design_matrices([di], pd.DataFrame([coded_row(r, coding) for r in rows]))[0])
    cov, beta = np.asarray(res.cov_params()), np.asarray(res.params)
    eta, se = X @ beta, np.sqrt(np.einsum("ij,jk,ik->i", X, cov, X))
    glm = resp == "interfrac"
    q = stats.norm.ppf(0.975) if glm else stats.t.ppf(0.975, res.df_resid)
    inv = expit if glm else np.exp
    return inv(eta), inv(eta - q * se), inv(eta + q * se)


def contrast_ratio(resp, rows_num, rows_den):
    """exp(eta_num - eta_den) with 95 % CI from the coefficient covariance (OLS: response ratio)."""
    df, res = MODELS_FIT[resp]
    di = res.model.data.design_info

    def X(rows):
        return np.asarray(patsy.build_design_matrices([di], pd.DataFrame([coded_row(r) for r in rows]))[0])

    c = X(rows_num) - X(rows_den)
    cov, beta = np.asarray(res.cov_params()), np.asarray(res.params)
    d, se = c @ beta, np.sqrt(np.einsum("ij,jk,ik->i", c, cov, c))
    q = stats.t.ppf(0.975, res.df_resid)
    return np.exp(d), np.exp(d - q * se), np.exp(d + q * se)


def fold_text(r):
    if 0.97 <= r <= 1.03:
        return f"×{r:.2f}"
    return f"×{r:.2g}" if r < 10 else f"×{r:.0f}"


def log_axis(fig, p, **kw):
    fig.update_xaxes(type="log", tickvals=TICKS[p], ticktext=[f"{v:g}" for v in TICKS[p]],
                     range=np.log10(RANGES[p]) + np.array([-0.04, 0.04]), minor=dict(showgrid=False),
                     title_text=f"{SHORTLAB[p]} ({UNITS[p]}, log)", **kw)


# ---- 1. re-exports ------------------------------------------------------------------------
def reexports():
    for f in (F5.coefficient_plot, F5.predicted_vs_heldout, F5.response_surfaces, F5.residual_diagnostics,
              F5.noise_check, F5.importance_comparison):
        f()
    # QC time series (frozen figure, not model-dependent): same runs, seed and equilibration frames
    eq = importlib.import_module("02_equilibration")
    inv = pd.read_csv(ANALYSIS / "run_inventory.csv")
    slugs = inv.loc[inv.status == "finished", "slug"].tolist()
    rng = np.random.default_rng(eq.SEED)
    pick = sorted(rng.choice(slugs, 12, replace=False))
    t0 = pd.read_csv(ANALYSIS / "qc" / "equilibration_per_run.csv").query("observable == 'inter_free_12A'") \
        .set_index("slug")["t0_frame"]
    fig = St.small_multiples(pick, lambda s: eq.series(np.load(TIMESERIES / f"{s}.npz"))["inter_free_12A"], t0,
                             title="Inter-chain free K–K contacts", ytitle="contacts / frame")
    save(fig, "qc_timeseries_sample_inter_free_12A")


# ---- 2. effect_curves ---------------------------------------------------------------------
def effect_curves():
    df, _ = MODELS_FIT["inter"]
    cb = coded_row(BASE)
    fig = make_subplots(rows=1, cols=4, shared_yaxes=True, horizontal_spacing=0.025,
                        subplot_titles=[SHORTLAB[p] for p in M.PARAMS])
    base_pred = predict("inter", [BASE])[0][0]
    rows, ndots, ends = [], {}, {}
    ymin, ymax = np.inf, 0.0
    for i, p in enumerate(M.PARAMS, 1):
        v = np.geomspace(*RANGES[p], 80)
        mid, lo, hi = predict("inter", [{**BASE, p: x} for x in v])
        ymin, ymax = min(ymin, lo.min()), max(ymax, hi.max())
        fig.add_trace(go.Scatter(x=np.r_[v, v[::-1]], y=np.r_[hi, lo[::-1]], fill="toself", mode="lines",
                                 fillcolor="rgba(31,107,255,0.20)", line=dict(width=0), hoverinfo="skip",
                                 showlegend=False, name=f"95 % CI ({p})"), 1, i)
        fig.add_trace(go.Scatter(x=v, y=mid, mode="lines", line=dict(color=St.BLUE, width=3), showlegend=False,
                                 name=f"model ({p})"), 1, i)
        others = [q for q in M.PARAMS if q != p]
        near = np.all([np.abs(df[M.SHORT[q]] - cb[M.SHORT[q]]) <= NEAR_BASE for q in others], axis=0)
        ndots[p] = int(near.sum())
        fig.add_trace(go.Scatter(x=df.loc[near, p], y=df.loc[near, "y"], mode="markers", showlegend=False,
                                 name=f"simulated runs ({p})", text=df.loc[near, "slug"],
                                 marker=dict(color=St.BLUE, size=7, opacity=0.55, line=dict(width=0))), 1, i)
        fig.add_trace(go.Scatter(x=[BASE[p]], y=[base_pred], mode="markers", showlegend=False, name=f"base case ({p})",
                                 marker=dict(color=St.ORANGE, size=11, symbol="circle",
                                             line=dict(color="rgba(0,0,0,0)", width=0))), 1, i)
        rising = mid[-1] > mid[0]
        label_pos = {"low": "bottom right" if rising else "top right", "high": "top left" if rising else "bottom left"}
        for side, x, y in (("low", v[0], mid[0]), ("high", v[-1], mid[-1])):
            r = y / base_pred
            ends[(p, side)] = r
            fig.add_trace(go.Scatter(x=[x], y=[y], mode="text", text=[fold_text(r)], showlegend=False,
                                     textposition=label_pos[side],
                                     textfont=dict(color=St.INK, size=13), hoverinfo="skip",
                                     name=f"fold change vs base ({p}, {side})"), 1, i)
        rows.append({"parameter": p, "low_value": RANGES[p][0], "high_value": RANGES[p][1], "base_value": BASE[p],
                     "base_pred": base_pred, "low_pred": mid[0], "high_pred": mid[-1],
                     "fold_low_vs_base": ends[(p, "low")], "fold_high_vs_base": ends[(p, "high")],
                     "n_runs_near_base": ndots[p]})
    St.layout(fig, None, width=1200, height=470)
    fig.update_layout(margin=dict(l=80, r=20, t=50, b=80))
    for i, p in enumerate(M.PARAMS, 1):
        log_axis(fig, p, col=i)
    fig.update_yaxes(type="log", range=np.log10([ymin * 0.75, ymax * 1.5]),
                     tickvals=[0.05, 0.1, 0.2, 0.5, 1], ticktext=["0.05", "0.1", "0.2", "0.5", "1"], minor=dict(showgrid=False))
    fig.update_yaxes(title_text="stiffness proxy (K–K contacts / frame, log)", col=1)
    fig.update_annotations(selector=dict(showarrow=False), font=dict(family=St.FONT, size=14, color=St.INK))
    pd.DataFrame(rows).to_csv(DATA / "effect_curves_endpoints.csv", index=False)
    save(fig, "effect_curves", {"n_runs_near_base": ndots, "near_base_rule": f"|coded log value - base| <= {NEAR_BASE} on the other three"})
    return ndots


# ---- 3. knobs_fold_change -----------------------------------------------------------------
def knobs_fold_change():
    oat = pd.read_csv(MODELS / "inter_log" / "effects_oat.csv").set_index("parameter")
    _, res = MODELS_FIT["inter"]
    recs = []
    for p in M.PARAMS:
        lo_r = contrast_ratio("inter", [{**BASE, p: RANGES[p][0]}], [BASE])
        hi_r = contrast_ratio("inter", [{**BASE, p: RANGES[p][1]}], [BASE])
        # the same predictions as models/inter_log/effects_oat.csv
        assert np.isclose(lo_r[0][0] * oat.loc[p, "base_pred"], oat.loc[p, "low_pred"], rtol=1e-8)
        assert np.isclose(hi_r[0][0] * oat.loc[p, "base_pred"], oat.loc[p, "high_pred"], rtol=1e-8)
        recs.append({"parameter": p, "low_vs_base": lo_r[0][0], "low_ci_low": lo_r[1][0], "low_ci_high": lo_r[2][0],
                     "high_vs_base": hi_r[0][0], "high_ci_low": hi_r[1][0], "high_ci_high": hi_r[2][0],
                     "high_over_low": oat.loc[p, "ratio_high_over_low"]})
    d = pd.DataFrame(recs)
    d["span"] = np.abs(np.log(d.high_over_low))
    d = d.sort_values("span").reset_index(drop=True)
    y = np.arange(len(d))
    fig = go.Figure()
    fig.add_trace(go.Bar(y=y, x=(d[["low_vs_base", "high_vs_base"]].max(axis=1) - d[["low_vs_base", "high_vs_base"]].min(axis=1)),
                         base=d[["low_vs_base", "high_vs_base"]].min(axis=1), orientation="h", width=0.5,
                         marker=dict(color="rgba(31,107,255,0.55)", line=dict(width=0)), showlegend=False,
                         name="low end to high end of range"))
    for side, c in (("low", St.INK), ("high", St.INK)):
        fig.add_trace(go.Scatter(x=d[f"{side}_vs_base"], y=y, mode="markers", showlegend=False, name=f"{side} end",
                                 marker=dict(color=St.BLUE, size=9, symbol="line-ns", line=dict(color=St.BLUE, width=3)),
                                 error_x=dict(type="data", symmetric=False, color=St.BLUE, thickness=2, width=5,
                                              array=d[f"{side}_ci_high"] - d[f"{side}_vs_base"],
                                              arrayminus=d[f"{side}_vs_base"] - d[f"{side}_ci_low"])))
    for yi, r in zip(y, d.itertuples()):
        fig.add_annotation(x=1, xref="paper", y=yi, xanchor="left", xshift=8, showarrow=False,
                           text=f"{fold_text(r.high_over_low)} low → high", font=dict(color=St.INK, size=13))
    fig.add_vline(x=1, line=dict(color=St.ORANGE, width=2))
    St.layout(fig, None, width=1000, height=380)
    fig.update_layout(margin=dict(l=170, r=170, t=30, b=70))
    fig.update_yaxes(tickvals=y, ticktext=[f"{SHORTLAB[p]} ({RANGES[p][0]:g}–{RANGES[p][1]:g} {UNITS[p]})" for p in d.parameter],
                     showgrid=False)
    fig.update_xaxes(type="log", tickvals=[0.1, 0.2, 0.5, 1, 2, 5], ticktext=["0.1", "0.2", "0.5", "1", "2", "5"],
                     range=np.log10([0.08, 6]), title_text="stiffness proxy relative to the base case (log); orange line = base = 1")
    d.drop(columns="span").to_csv(DATA / "knobs_fold_change.csv", index=False)
    save(fig, "knobs_fold_change")


# ---- 4. response_surface_base with tissue bands -------------------------------------------
def tissue_bands(p0, map_lo, map_hi):
    rows = []
    for k in (1, 2):
        for t, (elo, ehi) in TISSUES.items():
            glo, ghi = elo / E_OVER_G, ehi / E_OVER_G
            plo, phi = (p0 * (g / G0_KPA) ** (1.0 / k) for g in (glo, ghi))
            rows.append({"tissue": t, "k": k, "E_low_kPa": elo, "E_high_kPa": ehi, "G_low_kPa": glo, "G_high_kPa": ghi,
                         "proxy_low": plo, "proxy_high": phi, "map_min": map_lo, "map_max": map_hi,
                         "reachable_on_map": "yes" if (phi >= map_lo and plo <= map_hi) else "no",
                         "fully_inside_map": plo >= map_lo and phi <= map_hi, "values_unverified": True,
                         "base_proxy": p0, "G0_kPa": G0_KPA, "E_over_G": E_OVER_G})
    return pd.DataFrame(rows)


def response_surface_base():
    df, res = MODELS_FIT["inter"]
    px, py = "concentration", "spacing_pct"
    zmax = 1.4
    cont, arrows = F5._surface_traces(res, px, py, {}, 0.0, zmax, True)
    p0 = float(predict("inter", [BASE])[0][0])
    map_lo, map_hi = float(np.min(cont.z)), float(np.max(cont.z))
    tb = tissue_bands(p0, map_lo, map_hi)
    tb.to_csv(DATA / "tissue_bands.csv", index=False)

    W, H, l, r, t, b = 1100, 700, 80, 290, 150, 140
    Wp, Hp = W - l - r, H - t - b
    fig = go.Figure(cont)
    for a in arrows:
        a.update(xref="x", yref="y", axref="x", ayref="y")
    St.layout(fig, None, width=W, height=H)
    fig.update_layout(margin=dict(l=l, r=r, t=t, b=b), annotations=arrows)
    fig.update_traces(colorbar=dict(title=dict(text=""), thickness=16, lenmode="fraction", len=1,
                                    yanchor="bottom", y=0, ypad=0, xanchor="left", x=1 + 14 / Wp, xpad=0,
                                    tickfont=dict(color=St.INK), outlinewidth=0), selector=dict(type="contour"))
    fig.update_xaxes(title_text=M.LABEL[px], showgrid=False)
    fig.update_yaxes(title_text=M.LABEL[py], showgrid=False)
    x0 = 1 + 78 / Wp                        # bands start right of the colour bar and its tick labels
    colw, gap = 16 / Wp, 12 / Wp
    cols = [St.BLUE, St.GREEN, St.PURPLE, St.ORANGE]
    k1 = tb[tb.k == 1].reset_index(drop=True)
    for j, row in k1.iterrows():
        xa = x0 + j * (colw + gap)
        lo, hi = max(row.proxy_low, 0.0), min(row.proxy_high, zmax)
        shape = dict(type="rect", xref="paper", yref="paper", x0=xa, x1=xa + colw, y0=lo / zmax, y1=max(hi / zmax, lo / zmax + 0.004),
                     fillcolor=cols[j], opacity=0.85, line=dict(width=0))
        fig.add_shape(**shape)
        off = row.reachable_on_map == "no"
        fig.add_annotation(xref="paper", yref="paper", x=xa + colw / 2, y=1.015, xanchor="center", yanchor="bottom",
                           textangle=-90, showarrow=False, font=dict(color=cols[j], size=13),
                           text=f"{row.tissue}, E {row.E_low_kPa:g}–{row.E_high_kPa:g} kPa" + (" (below map)" if off else ""))
    # the map's own range as a thin grey bracket in a fifth column
    xm = x0 + 4 * (colw + gap)
    fig.add_shape(type="rect", xref="paper", yref="paper", x0=xm + colw / 2 - 1.5 / Wp, x1=xm + colw / 2 + 1.5 / Wp,
                  y0=map_lo / zmax, y1=map_hi / zmax, fillcolor=St.MUTED, line=dict(width=0), opacity=0.9)
    fig.add_annotation(xref="paper", yref="paper", x=xm + colw / 2, y=1.015, xanchor="center", yanchor="bottom",
                       textangle=-90, showarrow=False, font=dict(color=St.MUTED, size=13), text="map range")
    fig.add_annotation(xref="paper", yref="paper", x=xm + colw + 36 / Wp, y=0.5, xanchor="center", showarrow=False,
                       textangle=-90, font=dict(color=St.INK, size=14), text=F5.YLAB)
    fig.add_annotation(xref="paper", yref="paper", x=0, y=-(b - 85) / Hp, xanchor="left", yanchor="top", showarrow=False,
                       align="left", font=dict(color=St.INK, size=13),
                       text="Illustrative: assumes G is proportional to the proxy and G'0 = 3 kPa.<br>Needs one rheology measurement to calibrate.")
    save(fig, "response_surface_base", {"tissue_bands_k1": k1.to_dict("records"), "calibration": f"G'0={G0_KPA} kPa at proxy {p0:.4f}, E=3G, k=1"})
    return tb


# ---- 5. bridges_vs_loops ------------------------------------------------------------------
def bridges_vs_loops():
    inter = pd.read_csv(MODELS / "inter_log" / "effects_oat.csv").set_index("parameter")
    intra = pd.read_csv(MODELS / "intra_log" / "effects_oat.csv").set_index("parameter")
    order = list(inter.ratio_high_over_low.abs().apply(np.log).abs().sort_values().index)
    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.17,
                        subplot_titles=["fold change across each parameter's range", "bridge share at base chain length and ionic strength<br>"
                                       "<span style='color:#1F6BFF'>blue line: 50 % (bridges = loops)</span>, × = base case"])
    ypos = np.arange(len(order))
    recs = []
    for off, (lab, tab, col) in zip((+0.19, -0.19), [("bridges (inter-chain contacts)", inter, St.BLUE),
                                                      ("loops (intra-chain contacts)", intra, St.ORANGE)]):
        v = np.array([tab.loc[p, "ratio_high_over_low"] for p in order])
        lo = np.array([tab.loc[p, "ratio_ci_low"] for p in order])
        hi = np.array([tab.loc[p, "ratio_ci_high"] for p in order])
        fig.add_trace(go.Bar(y=ypos + off, x=v - 1, base=1, orientation="h", width=0.34, name=lab,
                             marker=dict(color=col, line=dict(width=0)), showlegend=True), 1, 1)
        fig.add_trace(go.Scatter(x=v, y=ypos + off, mode="markers", showlegend=False, name=f"{lab} 95 % CI",
                                 marker=dict(color=col, size=1, opacity=0),
                                 error_x=dict(type="data", symmetric=False, color=col, thickness=1.6, width=4,
                                              array=hi - v, arrayminus=v - lo)), 1, 1)
        fig.add_trace(go.Scatter(x=np.where(v >= 1, hi, lo), y=ypos + off, mode="text", showlegend=False, name=f"{lab} fold change",
                                 hoverinfo="skip",
                                 text=[fold_text(x) for x in v], textposition=["middle right" if x >= 1 else "middle left" for x in v],
                                 textfont=dict(color=St.INK, size=12)), 1, 1)
        recs += [{"series": lab, "parameter": p, "ratio_high_over_low": a, "ci_low": b_, "ci_high": c_}
                 for p, a, b_, c_ in zip(order, v, lo, hi)]
    fig.add_vline(x=1, line=dict(color=St.MUTED, width=1.5), row=1, col=1)
    fig.update_yaxes(tickvals=ypos, ticktext=[SHORTLAB[p] for p in order], showgrid=False, row=1, col=1)
    fig.update_xaxes(type="log", range=np.log10([0.025, 30]), tickvals=[0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20],
                     ticktext=["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10", "20"], row=1, col=1,
                     title_text="fold change, high end / low end of range (log, 1 = no change)")

    # panel B: quasi-binomial (logit) GLM bridge share, density x spacing
    _, res = MODELS_FIT["interfrac"]
    xs = np.linspace(*RANGES["concentration"], 61)
    ys = np.linspace(*RANGES["spacing_pct"], 61)
    X, Y = np.meshgrid(xs, ys)
    pts = [{**BASE, "concentration": a, "spacing_pct": b_} for a, b_ in zip(X.ravel(), Y.ravel())]
    z = predict("interfrac", pts)[0].reshape(X.shape)
    fig.add_trace(go.Contour(x=xs, y=ys, z=z, colorscale=St.SURFACE_SCALE, zmin=0, zmax=1, showscale=True,
                             contours=dict(start=0.1, end=0.9, size=0.1, showlabels=True,
                                           labelfont=dict(family=St.FONT, size=11, color=St.INK), labelformat=".0%"),
                             line=dict(color="rgba(124,128,137,0.75)", width=1), name="bridge share",
                             colorbar=dict(title=dict(text="bridge share, inter / (inter + intra)", side="right"), thickness=14,
                                           tickformat=".0%", tickfont=dict(color=St.INK), outlinewidth=0, x=1.0, xpad=0, len=1,
                                           yanchor="bottom", y=0, ypad=0)), 1, 2)
    fig.add_trace(go.Contour(x=xs, y=ys, z=z, showscale=False, hoverinfo="skip", name="50 % line",
                             contours=dict(coloring="lines", start=0.5, end=0.5, size=1, showlabels=False),
                             line=dict(color=St.BLUE, width=4), colorscale=[[0, St.BLUE], [1, St.BLUE]]), 1, 2)
    fig.add_trace(go.Scatter(x=[BASE["concentration"]], y=[BASE["spacing_pct"]], mode="markers", showlegend=False,
                             name="base case", marker=dict(color=St.INK, size=9, symbol="x-thin", line=dict(color=St.INK, width=2.5))), 1, 2)
    fig.update_xaxes(title_text=M.LABEL["concentration"], showgrid=False, row=1, col=2)
    fig.update_yaxes(title_text=M.LABEL["spacing_pct"], showgrid=False, row=1, col=2)
    St.layout(fig, None, width=1300, height=520)
    fig.update_layout(margin=dict(l=130, r=20, t=80, b=130), legend=dict(orientation="h", x=0.0, y=-0.28, xanchor="left"),
                      bargap=0.1)
    fig.update_annotations(selector=dict(showarrow=False, xref="paper"), font=dict(family=St.FONT, size=14, color=St.INK))
    pd.DataFrame(recs).to_csv(DATA / "bridges_vs_loops_panelA.csv", index=False)
    inter_frac = pd.read_csv(MODELS / "interfrac_log" / "effects_oat.csv").set_index("parameter")
    save(fig, "bridges_vs_loops", {"bridge_share_base": float(inter_frac.base_pred.iloc[0]),
                                   "bridge_share_range": {p: [float(inter_frac.loc[p, "low_pred"]), float(inter_frac.loc[p, "high_pred"])] for p in M.PARAMS}})


# ---- 6. reactive_bridge_fraction ----------------------------------------------------------
LABEL_POS_REACTIVE = ["bottom right", "top left", "top left", "middle left", "top left", "top right",
                      "bottom right", "bottom right", "middle right", "bottom right"]     # A10..A19, hand-placed


def reactive_bridge_fraction():
    d = pd.read_csv(ANALYSIS / "stats" / "reactive_design_points.csv")
    rho, prho = stats.spearmanr(d.pred_bridge_share, d.bridge_frac_mean)
    r, pr = stats.pearsonr(d.pred_bridge_share, d.bridge_frac_mean)
    # side check: the log-coded quasi-binomial GLM used in bridges_vs_loops panel B
    g = predict("interfrac", [{p: row[p] for p in M.PARAMS} for _, row in d.iterrows()])[0]
    rho_g, _ = stats.spearmanr(g, d.bridge_frac_mean)
    r_g, _ = stats.pearsonr(g, d.bridge_frac_mean)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0.3, 1], y=[0.3, 1], mode="lines", line=dict(color=St.MUTED, width=1.5, dash="dot"),
                             showlegend=False, hoverinfo="skip", name="y = x"))
    fig.add_trace(go.Scatter(x=d.pred_bridge_share, y=d.bridge_frac_mean, mode="markers+text", text=d.anchor_id, showlegend=False,
                             textposition=LABEL_POS_REACTIVE, textfont=dict(color=St.MUTED, size=12), name="reactive design points",
                             marker=dict(color=St.BLUE, size=11, line=dict(width=0)),
                             error_y=dict(type="data", array=d.bridge_frac_sd, color=St.BLUE, thickness=1.8, width=4)))
    fig.add_annotation(xref="paper", yref="paper", x=0.02, y=0.98, xanchor="left", yanchor="top", showarrow=False, align="left",
                       text=f"Spearman ρ = {rho:.2f}<br>Pearson r = {r:.2f}", font=dict(color=St.INK, size=15))
    St.layout(fig, None, width=760, height=640)
    fig.update_layout(margin=dict(l=80, r=30, t=30, b=70))
    fig.update_xaxes(title_text="predicted bridge share, non-reactive models", range=[0.3, 1.0], tickformat=".0%")
    fig.update_yaxes(title_text="observed bridge fraction (mean ± SD of 3 runs)", range=[0.3, 1.05], tickformat=".0%")
    save(fig, "reactive_bridge_fraction", {"spearman": rho, "spearman_p": prho, "pearson": r, "pearson_p": pr,
                                           "n_design_points": len(d), "n_runs": 3 * len(d),
                                           "side_check_interfrac_log_glm": {"spearman": rho_g, "pearson": r_g}})
    return dict(spearman=rho, pearson=r, spearman_glm=rho_g, pearson_glm=r_g)


# ---- 7. model comparison ------------------------------------------------------------------
def model_comparison():
    rows = []
    for coding in ("linear", "log"):
        fit = pd.read_csv(MODELS / f"inter_{coding}" / "model_fit.csv").iloc[0]
        oat = pd.read_csv(MODELS / f"inter_{coding}" / "effects_oat.csv").set_index("parameter")
        pt = pd.read_csv(MODELS / f"inter_{coding}" / "parameter_tests.csv").set_index("parameter")
        df = M.load(coding=coding)
        sub = df[df.slug != EXCLUDE_RUN].reset_index(drop=True)
        lof_wo = ST4.lack_of_fit(sub, M.fit(sub))
        row = {"model": "frozen linear-coded baseline" if coding == "linear" else "log-coded, all four parameters (not frozen)",
               "coding": coding, "n_runs": int(fit.n_runs), "n_design_points": int(fit.n_design_points),
               "R2": fit.R2, "adj_R2": fit.adj_R2, "cv5_grouped_R2": fit.cv5_grouped_R2_ols,
               "lof_F": fit.lof_F, "lof_df_num": int(fit.lof_df_num), "lof_df_den": int(fit.lof_df_den), "lof_p": fit.lof_p,
               "n_runs_without_extra2_05": len(sub), "lof_F_without_extra2_05": lof_wo["F"],
               "lof_df_num_without_extra2_05": lof_wo["df_lack_of_fit"], "lof_p_without_extra2_05": lof_wo["p"]}
        for p in M.PARAMS:
            row[f"fold_{p}"] = oat.loc[p, "ratio_high_over_low"]
            row[f"fold_{p}_ci_low"] = oat.loc[p, "ratio_ci_low"]
            row[f"fold_{p}_ci_high"] = oat.loc[p, "ratio_ci_high"]
        for p in M.PARAMS:
            row[f"lmg_share_{p}"] = pt.loc[p, "lmg_share_of_explained"]
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(DATA / "model_comparison.csv", index=False)
    return out


# ---- main ---------------------------------------------------------------------------------
def main():
    global MODELS_FIT
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    MODELS_FIT = fit_models()
    reexports()
    ndots = effect_curves()
    knobs_fold_change()
    tb = response_surface_base()
    bridges_vs_loops()
    rx = reactive_bridge_fraction()
    mc = model_comparison()
    info = {"runs_in_fits": {k: len(v[0]) for k, v in MODELS_FIT.items()}, "near_base_dots": ndots, "reactive": rx}
    (DATA / "run_info.json").write_text(json.dumps(info, indent=1, default=float))
    print(json.dumps(info, indent=1, default=float))
    print(tb.round(4).to_string())
    print(mc.T.to_string())


if __name__ == "__main__":
    main()
