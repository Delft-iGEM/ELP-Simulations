"""MOSAIC figure style for Plotly: Poppins, team palette, no subtitles or footnotes in figures.

Palette (validated with the dataviz checker on white): blue, orange, green and purple are
used for marks; yellow #FFAF00 fails the lightness/contrast check for thin marks (1.8:1),
so it only appears inside continuous colour scales. Green and orange are close for
deuteranopes (dE 7.2), so wherever both are plotted they also differ by marker symbol.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

YELLOW, ORANGE, GREEN, BLUE, PURPLE = "#FFAF00", "#FF3A00", "#00A445", "#1F6BFF", "#A855F7"
MARKS = [BLUE, ORANGE, GREEN, PURPLE]          # fixed order, never cycled
SYMBOLS = ["circle", "square", "diamond", "triangle-up"]
INK, MUTED, GRID = "#1d1d1f", "#6e6e73", "#e6e6e6"
# Sequential scale for the response surface: white -> yellow -> orange (one warm family).
SURFACE_SCALE = [[0.0, "#FFF6DD"], [0.45, YELLOW], [1.0, ORANGE]]
FONT = "Poppins, Arial, sans-serif"


def layout(fig, title=None, width=900, height=560):
    import plotly.graph_objects as go  # noqa: F401
    fig.update_layout(
        template="simple_white", font=dict(family=FONT, size=14, color=INK),
        title=dict(text=title, x=0.0, xanchor="left", font=dict(size=18)) if title else None,
        width=width, height=height, margin=dict(l=70, r=30, t=60 if title else 30, b=60),
        legend=dict(bgcolor="rgba(0,0,0,0)"), hoverlabel=dict(font_family=FONT),
    )
    fig.update_xaxes(showgrid=True, gridcolor=GRID, linecolor=MUTED, ticks="outside", tickcolor=MUTED)
    fig.update_yaxes(showgrid=True, gridcolor=GRID, linecolor=MUTED, ticks="outside", tickcolor=MUTED)
    return fig


def save(fig, stem: Path):
    """HTML always (interactive, hover); PNG too when kaleido is available."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(stem) + ".html", include_plotlyjs="cdn")
    try:
        fig.write_image(str(stem) + ".png", scale=2)
    except Exception as e:  # noqa: BLE001
        print(f"  (no PNG for {stem.name}: {e})")


def running_mean(x, w=100):
    c = np.cumsum(np.insert(np.asarray(x, float), 0, 0.0))
    out = np.full(len(x), np.nan)
    out[w - 1:] = (c[w:] - c[:-w]) / w
    return out


def small_multiples(slugs, get, t0_by_slug, title, ytitle, ncols=4, ns_per_frame=0.07):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    nrows = int(np.ceil(len(slugs) / ncols))
    short = [s.replace("sub4-", "") for s in slugs]
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=short, shared_xaxes=True,
                        horizontal_spacing=0.06, vertical_spacing=0.09)
    for i, s in enumerate(slugs):
        r, c = i // ncols + 1, i % ncols + 1
        y = get(s)
        t = np.arange(len(y)) * ns_per_frame
        fig.add_trace(go.Scattergl(x=t, y=y, mode="lines", line=dict(color="#c9d8ff", width=1),
                                   name="per frame", showlegend=i == 0,
                                   hovertemplate="%{x:.0f} ns: %{y:.2f}<extra></extra>"), r, c)
        fig.add_trace(go.Scattergl(x=t, y=running_mean(y), mode="lines", line=dict(color=BLUE, width=2),
                                   name="100-frame mean", showlegend=i == 0,
                                   hovertemplate="%{x:.0f} ns: %{y:.2f}<extra></extra>"), r, c)
        fig.add_vline(x=float(t0_by_slug[s]) * ns_per_frame, line=dict(color=ORANGE, width=2, dash="dash"),
                      row=r, col=c)
    layout(fig, title, width=1300, height=260 * nrows + 80)
    fig.update_annotations(font=dict(family=FONT, size=11, color=MUTED))
    fig.update_xaxes(title_text="time (ns)", row=nrows)
    fig.update_yaxes(title_text=ytitle, col=1)
    return fig
