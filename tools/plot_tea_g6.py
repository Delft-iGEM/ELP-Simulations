"""Figures for the gamma-6, 285 / 370 K set: Rg and SASA, grouped by what varied.

    .venv/bin/python tools/plot_tea_g6.py                 # reads figures/g6/g6-285-370.csv
    .venv/bin/python tools/plot_tea_g6.py --csv other.csv --out-prefix g6

Writes <prefix>-rg.png and <prefix>-sasa.png. Each has three panels, one per
thing that was varied:

  lysine placement   the 298-aa lysine scaffold, best 4 then worst 4 by
                     inter-chain share of crosslink hits
  RGD placement      the 357-aa construct, best 4 then worst 4 by 23-mer SASA,
                     with the original as reference
  8 lysines          the original with 8 VPGKG, grouped by design family, with
                     the 4-lysine original as reference

Each row is one sequence: a dot at 285 K, a dot at 370 K, joined so the
temperature response reads as the length of the bar. Error bars are block
standard errors (often smaller than the dot). A sequence whose runs have not
finished is listed and marked "pending" rather than dropped.

The lysine scaffold has no 23-mer cell-binding domain, only RGD tripeptides, so
its SASA panel shows tripeptide SASA on its own axis.

Palette #4269d0 / #c1442e is the project's validated pair (tools/plot_tea_temps.py):
cool for the cold run, warm for the hot one.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLD = "#4269d0"
HOT = "#c1442e"
INK = "#1a1a1a"
MUTED = "#6b6f76"
GRID = "#dcdee2"
SURFACE = "#fcfcfb"
TEMPS = ((285, COLD), (370, HOT))

# (header, [(stem, label), ...]) per panel; stems are run names without -<T>k.
LYS = [
    ("best 4 — inter-chain share", [
        ("g6-lys12-even-4", "lys12  even-4"),
        ("g6-lys13-even-4-wide", "lys13  even-4-wide"),
        ("g6-lys19-blockb", "lys19  block b"),
        ("g6-lys20-blockc", "lys20  block c"),
    ]),
    ("worst 4", [
        ("g6-lys02-blk4-mid", "lys02  block of 4, middle"),
        ("g6-lys03-blk4-cterm", "lys03  block of 4, C-term"),
        ("g6-lys09-blk13-nc", "lys09  1 + 3, N/C"),
        ("g6-lys01-blk4-nterm", "lys01  block of 4, N-term"),
    ]),
]
RGD = [
    ("reference", [("g6-orig", "original construct")]),
    ("best 4 — 23-mer SASA", [
        ("g6-rgd19-3n-1c", "rgd19  3 N + 1 C"),
        ("g6-rgd15-pairs-ends", "rgd15  pairs at ends"),
        ("g6-rgd14-tight-c", "rgd14  tight, C-term"),
        ("g6-rgd22-c-heavy", "rgd22  C-heavy"),
    ]),
    ("worst 4", [
        ("g6-rgd16-pairs-mid", "rgd16  pairs, middle"),
        ("g6-rgd10-clust-q1", "rgd10  clustered, 1st quarter"),
        ("g6-rgd03-even-inset", "rgd03  even, inset"),
        ("g6-rgd23-staircase", "rgd23  staircase"),
    ]),
]
K8 = [
    ("reference", [("g6-orig", "original (4 K)")]),
    ("spread", [
        ("g6-k8-02-even", "02  even"),
        ("g6-k8-03-gap-thirds", "03  2 per gap, thirds"),
        ("g6-k8-04-rgd-flank", "04  flanking each RGD"),
    ]),
    ("paired", [
        ("g6-k8-05-pairs-centre", "05  KK pair per gap"),
        ("g6-k8-01-orig-dup", "01  original motif ×2"),
    ]),
    ("clustered", [
        ("g6-k8-09-alt-mid", "09  alternating, middle"),
        ("g6-k8-08-blk44-ends", "08  4 + 4 at the ends"),
        ("g6-k8-07-blk8-mid", "07  block of 8, middle"),
        ("g6-k8-06-blk8-cterm", "06  block of 8, C-term"),
    ]),
    ("baseline", [("g6-k8-10-random", "10  random")]),
]
PANELS = [("Lysine placement  (298-aa scaffold, 4 K)", LYS),
          ("RGD placement  (357-aa construct)", RGD),
          ("8 lysines instead of 4  (357-aa construct)", K8)]


def load(path: Path) -> dict[str, dict[int, dict]]:
    out: dict[str, dict[int, dict]] = {}
    for row in csv.DictReader(open(path)):
        stem, t = row["run"].rsplit("-", 1)
        out.setdefault(stem, {})[int(t.rstrip("k"))] = {
            k: float(v) for k, v in row.items() if k != "run" and v not in ("", "nan")}
    return out


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="x", colors=MUTED, labelsize=8)
    ax.tick_params(axis="y", length=0, labelsize=8.5, colors=INK)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def panel(ax, title, groups, data, key, sem_key, xlabel):
    _style(ax)
    y, ticks, labels = 0.0, [], []
    for header, rows in groups:
        ax.text(0.0, y, header, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7.5, color=MUTED, style="italic",
                fontweight="bold")
        y += 1
        for stem, label in rows:
            ticks.append(y)
            labels.append(label)
            at = data.get(stem, {})
            pts = [(T, c, at[T][key], at[T].get(sem_key, 0.0))
                   for T, c in TEMPS if T in at and key in at[T]]
            if len(pts) == 2:
                ax.plot([pts[0][2], pts[1][2]], [y, y], color=GRID, linewidth=2,
                        zorder=1, solid_capstyle="round")
            for T, c, v, e in pts:
                ax.errorbar(v, y, xerr=e, fmt="o", ms=6.5, color=c, ecolor=c,
                            elinewidth=1, capsize=0, mec=SURFACE, mew=1.2, zorder=3)
            missing = [f"{T} K" for T, _ in TEMPS if T not in at]
            if missing:
                # Left edge: the 370 K end of a row, empty when 370 K is the one
                # missing, so the note never sits on a dot.
                ax.text(0.01, y, "pending: " + ", ".join(missing),
                        transform=ax.get_yaxis_transform(), ha="left", va="center",
                        fontsize=7.5, color=MUTED)
            y += 1
        y += 0.4
    ax.set_yticks(ticks, labels)
    ax.set_ylim(y - 0.6, -0.8)
    ax.set_title(title, loc="left", fontsize=10, color=INK, pad=6, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=8.5, color=MUTED)
    return y


def figure(data, specs, suptitle, out: Path, shared: list[tuple[int, ...]]):
    heights = [sum(len(r) + 1.4 for _, r in g) for _, g in PANELS]
    height = 0.24 * sum(heights) + 1.8
    fig, axes = plt.subplots(3, 1, figsize=(8.2, height),
                             gridspec_kw={"height_ratios": heights, "hspace": 0.45})
    fig.subplots_adjust(top=1 - 0.75 / height, bottom=0.45 / height)
    fig.patch.set_facecolor(SURFACE)
    for ax, (title, groups), (key, sem_key, xlabel) in zip(axes, PANELS, specs):
        panel(ax, title, groups, data, key, sem_key, xlabel)
    for group in shared:          # panels measuring the same quantity share an axis
        lo = min(axes[i].get_xlim()[0] for i in group)
        hi = max(axes[i].get_xlim()[1] for i in group)
        for i in group:
            axes[i].set_xlim(lo, hi)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=6.5, color=c, mec=SURFACE,
                          label=f"{T} K") for T, c in TEMPS]
    fig.legend(handles=handles, loc="upper right", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.98, 0.995), labelcolor=INK)
    fig.suptitle(suptitle, x=0.02, y=0.995, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=Path("figures/g6/g6-285-370.csv"))
    parser.add_argument("--out-prefix", default="figures/g6/g6-285-370")
    args = parser.parse_args(argv)
    data = load(args.csv)

    rg = ("rg_nm", "rg_sem_nm", "Rg (nm)")
    figure(data, [rg, rg, rg],
           "Radius of gyration, single chain, TEA γ = 6",
           Path(f"{args.out_prefix}-rg.png"), shared=[(0, 1, 2)])
    figure(data, [("rgd_sasa_nm2", "rgd_sasa_sem_nm2", "SASA of the 4 RGD tripeptides (nm²)"),
                  ("domain_sasa_nm2", "domain_sasa_sem_nm2",
                   "SASA of the 4 23-mer cell-binding domains (nm²)"),
                  ("domain_sasa_nm2", "domain_sasa_sem_nm2",
                   "SASA of the 4 23-mer cell-binding domains (nm²)")],
           "RGD exposure, single chain, TEA γ = 6",
           Path(f"{args.out_prefix}-sasa.png"), shared=[(1, 2)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
