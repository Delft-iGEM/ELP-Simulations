"""SASA of the lysine ramp at 290 K only, one figure per quantity.

    .venv/bin/python tools/plot_lysramp_290k.py

Reads figures/g6-kspc-kend/g6-lysramp.csv (tools/tea_lysramp_analysis.py) and
writes, next to it,

  g6-lysramp-290k-chain-sasa.png          whole-chain SASA, for light slides
  g6-lysramp-290k-chain-sasa-darkbg.png   the same, light text, for dark slides
  g6-lysramp-290k-domain-sasa.png         RGD + linker (23-mer) SASA, per copy

The chain figures have absolute SASA against % of the 49 replaceable pentads
replaced, a transparent background and no caption (for slides); the domain
figure has % change from the original sequence (kend, 4 lysines) against
lysine count.

The combined figures share their y axis with 370 K, which squashes the 290 K
lines; here they get the whole axis. Runs: the 357-aa RGD + linker construct
(4 x LDASTVYAVTGRGDSPASSAASA + ELP), single free chains, TEA gamma 6, 400 ns,
50-400 ns analysed. Two series: lysines spread evenly (spaced) and lysines added
from both chain ends (kend).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

CSV = Path("figures/g6-kspc-kend/g6-lysramp.csv")
T = 290
# ink, muted, grid, surface; a transparent figure takes the slide's background,
# so the dark-background version needs its own light text
LIGHT = {"ink": "#000000", "muted": "#000000", "grid": "#dcdee2", "surface": "#fcfcfb",
         "line": "#2a6fd1"}
DARK = {"ink": "#f2f2f2", "muted": "#b4b8be", "grid": "#4a4e55", "surface": "#1e1e1e",
        "line": "#5b9bf0"}
STYLE = {"kspc": ("lysines spread evenly", "-", "o"),
         "kend": ("lysines added from both chain ends", (0, (4, 2)), "s")}
N_REPLACEABLE = 49                               # pentads that can become VPGKG
# tag, column, sem column, y label, title ("" = none), y as % change from the original,
# x as % of the replaceable pentads, slide version (transparent, no caption)
FIGS = [
    ("chain-sasa", "chain_sasa_nm2", "chain_sasa_sem_nm2", "SASA per chain (nm²)",
     "", False, True, True),
    ("domain-sasa", "domain_sasa_per_site_nm2", "domain_sasa_sem_per_site_nm2",
     "SASA per domain, change vs original", "RGD + linker (23-mer) SASA vs lysine count, 290 K",
     True, False, False),
]


def draw(plt, rows, original, fig_spec, pal: dict, out: Path) -> None:
    tag, key, sem, ylab, title, y_pct, x_pct, slide = fig_spec
    ink, muted, grid, surface, line = (pal[k] for k in ("ink", "muted", "grid", "surface",
                                                        "line"))
    ref = float(original[key])
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    fig.patch.set_facecolor(surface)
    ax.set_facecolor(surface)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(ink)
    ax.tick_params(colors=muted, labelsize=9)
    ax.grid(color=grid, linewidth=0.7)
    ax.set_axisbelow(True)
    # on a transparent figure a marker ring in the surface colour would show as a
    # halo, so the ring goes
    ring = dict(mec=line, mew=0) if slide else dict(mec=surface, mew=1.2)
    for series, (_, ls, marker) in STYLE.items():
        s = sorted((r for r in rows if r["series"] == series), key=lambda r: int(r["n_lys"]))
        n = np.array([int(r["n_lys"]) for r in s])
        y = np.array([float(r[key]) for r in s])
        e = np.array([float(r[sem]) for r in s])
        if y_pct:
            y, e = 100 * (y / ref - 1), 100 * e / ref
        if x_pct:
            n = 100 * n / N_REPLACEABLE
        ax.plot(n, y, color=line, lw=2, ls=ls, zorder=2)
        ax.errorbar(n, y, yerr=e, fmt=marker, color=line, ms=6.5, elinewidth=1.3, capsize=2.5, capthick=1.3,
                    zorder=3, **ring)
    if y_pct:
        ax.axhline(0, color=muted, lw=1, ls=(0, (4, 3)), zorder=1)
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:+.1f}%" if round(v, 2) else "0%")
    counts = sorted({int(r["n_lys"]) for r in rows})
    if x_pct:
        ax.set_xticks([100 * c / N_REPLACEABLE for c in counts])
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
        ax.set_xlabel("Percentage of replaced pentades (VPGIG → VPGKG)", fontsize=9.5, color=ink)
    else:
        ax.set_xticks(counts)
        ax.set_xlabel(f"Lysines per chain (VPGKG pentads, of {N_REPLACEABLE} replaceable)",
                      fontsize=9.5, color=ink)
    ax.set_ylabel(ylab, fontsize=9.5, color=ink)
    handles = [plt.Line2D([], [], color=line, lw=2, ls=ls, marker=m, ms=6.5, label=lab, **ring)
               for lab, ls, m in STYLE.values()]
    ax.legend(handles=handles, handlelength=3.2, loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=8.6, labelcolor=ink)
    if title:
        fig.suptitle(title, x=0.02, y=0.98, ha="left", fontsize=12.5, color=ink, fontweight="bold")
    if not slide:
        fig.text(0.02, 0.012,
                 "357-aa RGD + linker construct (4 × LDASTVYAVTGRGDSPASSAASA + ELP), one free "
                 "chain per point, 290 K, TEA γ 6, 400 ns;\nSASA every 10th frame over 50–400 ns, "
                 "bar = 20-block SEM. "
                 + ("8% (4 lysines)" if x_pct else "n = 4")
                 + " in the 'both ends' series is the original sequence"
                 + (f";\n% = change from the original's {ref:.2f} nm² per domain "
                    "(dashed line = 0)." if y_pct else "."),
                 fontsize=7.8, color=muted, linespacing=1.5)
    fig.tight_layout(rect=(0, 0 if slide else 0.11 if y_pct else 0.08, 1, 0.94 if title else 0.97))
    if slide:
        fig.savefig(out, dpi=200, transparent=True)
    else:
        fig.savefig(out, dpi=200, facecolor=surface)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [r for r in csv.DictReader(open(CSV)) if int(float(r["temperature_K"])) == T]
    # the original sequence: 4 lysines, kend series
    original = next(r for r in rows if r["series"] == "kend" and int(r["n_lys"]) == 4)
    for spec in FIGS:
        tag, slide = spec[0], spec[-1]
        draw(plt, rows, original, spec, LIGHT, CSV.parent / f"g6-lysramp-290k-{tag}.png")
        if slide:
            draw(plt, rows, original, spec, DARK,
                 CSV.parent / f"g6-lysramp-290k-{tag}-darkbg.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
