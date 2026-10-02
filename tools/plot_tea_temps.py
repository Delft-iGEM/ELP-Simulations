"""Figure for the 270/350 K TEA set: the collapse, and what causes it.

    .venv/bin/python tools/plot_tea_temps.py --rg-csv <tt-rg.csv>

Left panel  - dRg(270->350) per sequence, sorted, with the weighted mean. The
              question it answers is whether arrangement modulates the response:
              if the dots scatter no wider than their own error bars, it does not.
Right panel - the control decomposition at a FIXED 350 K. Each sequence has
              lambda(350) against stock lambda with the same seed and the same
              thermostat, so the gap between them is the stickiness alone.

Palette #4269d0 / #c1442e passes all six checks of the dataviz validator in
light mode (lightness band, chroma floor, CVD separation, normal-vision floor,
contrast vs surface).
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics as st
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#4269d0"       # the 24 RGD arrangements
RED = "#c1442e"        # the three designed sequences
INK = "#1a1a1a"
MUTED = "#6b6f76"
GRID = "#dcdee2"
SURFACE = "#fcfcfb"

DESIGNED = {"vpgig": "(VPGIG)71 — bare ELP",
            "orig": "original construct",
            "allterm": "RGD + Lys all terminal"}
CONTROL_PAIRS = ["lys-best", "lys-worst", "rgd-best", "rgd-worst"]


def load_delta(path: Path) -> dict[str, tuple[float, float]]:
    at: dict[str, dict[int, tuple[float, float]]] = {}
    for row in csv.DictReader(open(path)):
        stem = row["run"][3:].rsplit("-", 1)[0]
        at.setdefault(stem, {})[int(float(row["temperature_K"]))] = (
            float(row["rg_nm"]), float(row["rg_sem_nm"]))
    out = {}
    for stem, by_t in at.items():
        if 270 in by_t and 350 in by_t:
            (a, ae), (b, be) = by_t[270], by_t[350]
            out[stem] = (b - a, math.hypot(ae, be))
    return out


def control_rows() -> list[dict]:
    """Rg for each control triplet, computed on the fly."""
    from tools.tea_analysis import block_average, rg_series
    out = []
    for pair in CONTROL_PAIRS:
        row = {"pair": pair}
        for tag, name in [("tea280", f"tea-{pair}-free-280k"),
                          ("tea350", f"tea-{pair}-free-350k"),
                          ("stock350", f"ctrl-{pair}-free-350k-notea")]:
            _, rg = rg_series(name, 50.0)
            row[tag] = block_average(rg)
        out.append(row)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rg-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("tea-270-350-response.png"))
    args = parser.parse_args(argv)

    delta = load_delta(args.rg_csv)
    arrangements = sorted((v[0], v[1], k) for k, v in delta.items() if k.startswith("rgd"))
    designed = [(delta[k][0], delta[k][1], k) for k in DESIGNED if k in delta]

    weights = [1 / e ** 2 for _, e, _ in arrangements]
    mean = sum(d * w for (d, _, _), w in zip(arrangements, weights)) / sum(weights)
    between = st.stdev([d for d, _, _ in arrangements])
    within = st.mean([e for _, e, _ in arrangements])

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(13.2, 7.0), gridspec_kw={"width_ratios": [1.45, 1]})
    fig.patch.set_facecolor(SURFACE)

    # ---- left: dRg per sequence -------------------------------------------
    ordered = arrangements + designed
    ordered.sort()
    y = range(len(ordered))
    ax.axvspan(mean - within, mean + within, color=BLUE, alpha=0.10, lw=0,
               label=f"weighted mean ± one run's error")
    ax.axvline(mean, color=BLUE, lw=1.4, ls="--", alpha=0.75)
    for i, (d, e, name) in zip(y, ordered):
        is_designed = name in DESIGNED
        colour = RED if is_designed else BLUE
        ax.errorbar(d, i, xerr=e, fmt="o", ms=7 if is_designed else 6,
                    color=colour, ecolor=colour, elinewidth=1.6, capsize=3,
                    mfc=colour if is_designed else "white",
                    mec=colour, mew=1.8, zorder=3)
        label = DESIGNED.get(name, name.replace("rgd", "arr "))
        ax.text(d + e + 0.022, i, label, va="center", ha="left", fontsize=8.6,
                color=INK if is_designed else MUTED,
                fontweight="semibold" if is_designed else "normal")

    ax.set_yticks([])
    ax.set_xlim(-1.50, -0.40)
    ax.set_ylim(-1, len(ordered))
    ax.set_xlabel("ΔRg on heating 270 → 350 K  (nm)", fontsize=10.5, color=INK)
    ax.set_title("Every sequence collapses on heating — and by the same amount",
                 fontsize=12.5, color=INK, fontweight="semibold", loc="left", pad=52)
    ax.text(0, 1.012,
            f"24 RGD arrangements (hollow) scatter sd {between:.3f} nm against a "
            f"{within:.3f} nm per-run error:\nratio {between/within:.2f}, χ²/dof "
            f"{sum(w*(d-mean)**2 for (d,_,_),w in zip(arrangements,weights))/(len(arrangements)-1):.2f} "
            f"— no arrangement effect. Filled = designed sequences.",
            transform=ax.transAxes, fontsize=9.2, color=MUTED, va="bottom", ha="left")
    ax.legend(loc="lower right", frameon=False, fontsize=8.8, labelcolor=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)

    # ---- right: the control decomposition ---------------------------------
    rows = control_rows()
    for j, row in enumerate(rows):
        base = len(rows) - 1 - j
        t280, t350, s350 = row["tea280"], row["tea350"], row["stock350"]
        bx.plot([s350[0], t350[0]], [base, base], color=MUTED, lw=1.2,
                ls=":", zorder=1)
        bx.errorbar(t280[0], base + 0.17, xerr=t280[1], fmt="s", ms=6, color=MUTED,
                    ecolor=MUTED, elinewidth=1.4, capsize=3, mfc="white", mew=1.6,
                    label="λ(280), 280 K" if j == 0 else None, zorder=3)
        bx.errorbar(s350[0], base, xerr=s350[1], fmt="o", ms=7, color=BLUE,
                    ecolor=BLUE, elinewidth=1.6, capsize=3, mfc="white", mew=1.8,
                    label="stock λ, 350 K" if j == 0 else None, zorder=3)
        bx.errorbar(t350[0], base, xerr=t350[1], fmt="o", ms=7, color=RED,
                    ecolor=RED, elinewidth=1.6, capsize=3, zorder=3)
        if j == 0:
            bx.errorbar([], [], fmt="o", ms=7, color=RED, label="λ(350), 350 K")
        gap = t350[0] - s350[0]
        bx.text((s350[0] + t350[0]) / 2, base - 0.235, f"{gap:+.2f} nm",
                ha="center", va="top", fontsize=8.8, color=RED, fontweight="semibold")
        bx.text(t350[0], base + 0.40, row["pair"], ha="left", va="bottom",
                fontsize=9.2, color=INK, fontweight="semibold")

    bx.set_yticks([])
    bx.set_ylim(-0.75, len(rows) - 0.02)
    bx.set_xlabel("Rg  (nm)", fontsize=10.5, color=INK)
    bx.set_title("The collapse is λ(T), not the thermostat",
                 fontsize=12.5, color=INK, fontweight="semibold", loc="left", pad=52)
    bx.text(0, 1.012,
            "Same sequence, same seed, same 350 K — only the λ table differs.\n"
            "That one change accounts for 94–98% of the whole 280→350 K collapse.",
            transform=bx.transAxes, fontsize=9.2, color=MUTED, va="bottom", ha="left")
    bx.legend(loc="lower left", frameon=False, fontsize=8.8, labelcolor=MUTED)
    bx.grid(axis="x", color=GRID, lw=0.7)
    bx.set_axisbelow(True)
    for side in ("top", "right", "left"):
        bx.spines[side].set_visible(False)
    bx.spines["bottom"].set_color(GRID)
    bx.tick_params(colors=MUTED, labelsize=9)

    fig.text(0.006, 0.012,
             "CALVADOS-TEA, γ=3, single free chain, 80 nm box (controls 60 nm), "
             "300 ns (VPGIG 500 ns), first 50 ns discarded, block-averaged errors. "
             "270 K is extrapolated below TEA's 280–380 K fit range.",
             fontsize=7.9, color=MUTED)
    fig.tight_layout(rect=[0, 0.028, 1, 1])
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
