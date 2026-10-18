"""Slide versions of the collapsed-dwell figures: one panel, untitled, transparent.

    .venv/bin/python tools/dwell_slides.py lysramp        # from the lysramp CSVs
    .venv/bin/python tools/dwell_slides.py rgd            # analyse + plot the RGD runs
    .venv/bin/python tools/dwell_slides.py rgd --plot-only

lysramp  the evenly spaced lysine series (g6-kspc*, 290 / 370 K) from
         figures/g6-kspc-kend/g6-lysramp-dwell*.csv (tools/tea_lysramp_dwells.py);
         writes g6-lysramp-dwells-{strip,collapse}-spaced.png next to them.
rgd      g6-orig and the eight g6-rgd* arrangements at 285 / 370 K (the runs of
         tools/g6_rgd_states.py). Dwells with the definitions of
         tools/tea_lysramp_dwells.py (hysteresis 0.95 / 1.05 Rg_theta, minimum 3
         frames, censoring at the window edges); the collapse fraction is the
         plain Rg < Rg_theta fraction. Writes figures/g6/g6-rgd-dwells.csv,
         g6-rgd-dwell-summary.csv and g6-rgd-dwells-{strip,collapse,rg}.png
         (rg: mean Rg from figures/g6/g6-285-370.csv, dotted Rg_theta).
         Arrangements are ordered by the free-chain RGD + linker ranking, best
         left; the original, free winner and free loser are in bold, the others
         carry their free-chain place (of 24); a chain sketch (yellow bar, red
         dot per RGD, N terminus left) sits under each name.

Strip: one faint point per dwell (jittered), filled marker = median, bar = IQR,
open marker = fewer than 10 dwells. Collapse: fraction of frames with
Rg < Rg_theta.
"""

from __future__ import annotations

import argparse
import csv
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_lysramp_dwells import (COLD, GRID, HOT, MIN_FRAMES, SEED, SMALL_N, analyse,
                                      summarise)

INK = "#000000"
LYS_DIR = Path("figures/g6-kspc-kend")
N_REPLACEABLE = 49
RGD_DIR = Path("figures/g6")
RGD_LABELS = ["rgd19-3n-1c", "rgd15-pairs-ends", "rgd14-tight-c", "rgd22-c-heavy", "orig",
              "rgd23-staircase", "rgd03-even-inset", "rgd10-clust-q1", "rgd16-pairs-mid"]
FREE_CSV = Path("figures/free-rgdarr/sasa-24-arrangements.csv")
RG_CSV = RGD_DIR / "g6-285-370.csv"
STATES_CSV = RGD_DIR / "g6-rgd-states.csv"          # domain_starts (0-based) per arrangement
CHAIN_FILL, CHAIN_EDGE, RGD_DOT = "#f5c518", "#b8920c", "#d7263d"
RGD_ROLE = {"orig": "original", "rgd19-3n-1c": "free winner",
            "rgd16-pairs-mid": "free loser"}


def _axes(plt, width: float, height: float):
    fig, ax = plt.subplots(figsize=(width, height))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK)
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def _save(fig, plt, ax, temps, out: Path, strip: bool, half: bool = False,
          theta: bool = False) -> None:
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, ms=7, label=f"{T} K")
               for T, c in temps]
    if half:
        handles.append(plt.Line2D([], [], color=INK, lw=1.2, ls=":", label="50% of the time"))
    if theta:
        handles.append(plt.Line2D([], [], color=INK, lw=1.2, ls=":",
                                  label="Rg,θ = 0.21·N$^{0.5}$"))
    if strip:
        handles += [plt.Line2D([], [], marker="o", ls="", color="#6b6f76", ms=7,
                               label="median, bar = IQR"),
                    plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#6b6f76", mew=1.6,
                               ms=7, label=f"n < {SMALL_N} dwells")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(handles), frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(out, dpi=200, transparent=True)
    plt.close(fig)
    print(f"wrote {out}")


def strip(plt, ax, cells, events, rng, half_width: float) -> None:
    """cells: (x, colour, summary row)."""
    for x, colour, s in cells:
        d = np.array([e["duration_ns"] for e in events
                      if e["run"] == s["run"] and e["rule"] == "hyst"])
        if len(d):
            ax.scatter(x + rng.uniform(-half_width, half_width, len(d)), d, s=6,
                       color=colour, alpha=0.25, linewidths=0, zorder=2)
        if int(s["n_dwells"]):
            small = int(s["n_dwells"]) < SMALL_N
            ax.plot([x, x], [float(s["q25_ns"]), float(s["q75_ns"])], color=colour, lw=2.2,
                    zorder=3, solid_capstyle="butt")
            ax.plot(x, float(s["median_ns"]), "o", ms=7.5, zorder=4,
                    mfc="none" if small else colour, mec=colour, mew=1.8)
    ax.set_yscale("log")
    ax.set_ylabel("Collapsed dwell (ns)", fontsize=10, color=INK)


def collapse(ax, cells, half: bool = False) -> None:
    if half:                                    # dotted guide at half the time collapsed
        ax.axhline(0.5, color=INK, lw=1.2, ls=":", zorder=1)
    for x, colour, s in cells:
        ax.plot(x, float(s["collapse_fraction"]), "o", ms=8, color=colour, mew=0, zorder=3)
    ax.set_ylim(-0.03, 1.0)
    ax.set_ylabel("Fraction of time collapsed", fontsize=10, color=INK)


def _read(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path)))


def lysramp(plt) -> None:
    summary = [s for s in _read(LYS_DIR / "g6-lysramp-dwell-summary.csv")
               if s["rule"] == "hyst" and s["series"] == "kspc"]
    events = [dict(e, duration_ns=float(e["duration_ns"]))
              for e in _read(LYS_DIR / "g6-lysramp-dwells.csv")]
    temps = sorted({int(s["T"]) for s in summary})
    temps = list(zip(temps, (COLD, HOT)))
    off = {temps[0][0]: -1.8, temps[1][0]: 1.8}
    colour = dict(temps)
    cells = [(100 * int(s["n_lys"]) / N_REPLACEABLE + off[int(s["T"])], colour[int(s["T"])], s)
             for s in summary]
    counts = sorted({int(s["n_lys"]) for s in summary})
    for kind in ("strip", "collapse"):
        fig, ax = _axes(plt, 7.4, 4.8 if kind == "strip" else 4.4)
        if kind == "strip":
            strip(plt, ax, cells, events, np.random.default_rng(SEED), 1.1)
        else:
            collapse(ax, cells, half=True)
        ax.set_xticks([100 * c / N_REPLACEABLE for c in counts])
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
        ax.set_xlim(2, 106)
        ax.set_xlabel("Percentage of replaced pentades (VPGIG → VPGKG)", fontsize=10,
                      color=INK)
        _save(fig, plt, ax, temps, LYS_DIR / f"g6-lysramp-dwells-{kind}-spaced.png",
              kind == "strip", half=kind == "collapse")


def rgd(plt, plot_only: bool) -> None:
    ev_csv, sum_csv = RGD_DIR / "g6-rgd-dwells.csv", RGD_DIR / "g6-rgd-dwell-summary.csv"
    if not plot_only:
        cells = [(lab, T) for lab in RGD_LABELS for T in (285, 370)]
        with Pool(len(cells)) as pool:
            results = pool.map(analyse, [f"g6-{lab}-{T}k" for lab, T in cells])
        rng = np.random.default_rng(SEED)
        events, summary = [], []
        for (lab, T), r in zip(cells, results):
            for rule, mn in (("hyst", MIN_FRAMES), ("bare", 1)):
                d = np.array([L for L, _ in r[rule]], dtype=float) * r["dt_ns"]
                cens = np.array([c for _, c in r[rule]], dtype=bool)
                events += [{"run": r["run"], "rule": rule, "duration_ns": L,
                            "censored": int(c)} for L, c in zip(d, cens)]
                summary.append({"run": r["run"], "label": lab, "T": T, "rule": rule,
                                "rg_theta_nm": r["theta_nm"],
                                "collapse_fraction": r["collapse_fraction"]}
                               | summarise(d, cens, mn, r["dt_ns"], rng))
        for path, rows in ((ev_csv, events), (sum_csv, summary)):
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
        print(f"wrote {ev_csv}, {sum_csv}")
    summary = [s for s in _read(sum_csv) if s["rule"] == "hyst"]
    events = [dict(e, duration_ns=float(e["duration_ns"])) for e in _read(ev_csv)]
    temps = [(285, COLD), (370, HOT)]
    off = {285: -0.17, 370: 0.17}
    x_of = {lab: i for i, lab in enumerate(RGD_LABELS)}
    cells = [(x_of[s["label"]] + off[int(s["T"])], dict(temps)[int(s["T"])], s)
             for s in summary]
    free = sorted(_read(FREE_CSV), key=lambda r: -float(r["motif_23mer_nm2"]))
    place = {r["run"]: i + 1 for i, r in enumerate(free)}         # 1 = most exposed of 24
    rg = {r["run"]: r for r in _read(RG_CSV)}
    for kind in ("strip", "collapse", "rg"):
        fig, ax = _axes(plt, 10.0, 5.0 if kind == "strip" else 4.6)
        ax.grid(axis="x", visible=False)
        if kind == "strip":
            strip(plt, ax, cells, events, np.random.default_rng(SEED), 0.1)
        elif kind == "collapse":
            collapse(ax, cells, half=True)
        else:
            rg_panel(ax, x_of, temps, off, rg)
        for lab in RGD_ROLE:
            ax.axvspan(x_of[lab] - 0.42, x_of[lab] + 0.42, color=GRID, alpha=0.45, lw=0,
                       zorder=0)
        ax.set_xticks(range(len(RGD_LABELS)))
        ticks = []
        for lab in RGD_LABELS:
            name = "01-orig" if lab == "orig" else lab[3:]
            n = place[name]
            note = {"orig": f"(original,\n{n}{_suffix(n)} place)"}.get(
                lab, f"({RGD_ROLE[lab]})" if lab in RGD_ROLE else f"placed #{n}")
            # three lines each: the original keeps its name on one line
            ticks.append((name if lab == "orig" else name.replace("-", "\n", 1)) + "\n" + note)
        ax.set_xticklabels(ticks, fontsize=8.5)
        sketches(ax, fig, x_of)
        for t, lab in zip(ax.get_xticklabels(), RGD_LABELS):
            if lab in RGD_ROLE:
                t.set_fontweight("bold")
        ax.set_xlim(-0.6, len(RGD_LABELS) - 0.4)
        ax.set_xlabel("RGD arrangement (by free-chain RGD + linker SASA, highest left)",
                      fontsize=10, color=INK)
        _save(fig, plt, ax, temps, RGD_DIR / f"g6-rgd-dwells-{kind}.png", kind == "strip",
              half=kind == "collapse", theta=kind == "rg")


def sketches(ax, fig, x_of) -> None:
    """A small chain sketch per arrangement under its (three-line) name:
    yellow bar = chain (N terminus left), red dot = the G of each RGD."""
    from matplotlib.transforms import ScaledTranslation, blended_transform_factory
    starts = {r["label"]: (int(r["n_residues"]), [int(v) for v in r["domain_starts"].split()])
              for r in _read(STATES_CSV)}
    tr = (blended_transform_factory(ax.transData, ax.transAxes)
          + ScaledTranslation(0, -47 / 72, fig.dpi_scale_trans))
    half = 0.36
    for lab, x in x_of.items():
        n, dom = starts[lab]
        for lw, c, z in ((7, CHAIN_EDGE, 1), (5.8, CHAIN_FILL, 2)):
            ax.plot([x - half, x + half], [0, 0], color=c, lw=lw, solid_capstyle="round",
                    transform=tr, clip_on=False, zorder=z)
        g = [x - half + 2 * half * (d + 12) / (n - 1) for d in dom]
        ax.plot(g, [0] * len(g), "o", color=RGD_DOT, ms=3.6, mec=INK, mew=0.4,
                transform=tr, clip_on=False, zorder=3)
    ax.xaxis.labelpad = 22


def _suffix(n: int) -> str:
    return "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")


def rg_panel(ax, x_of, temps, off, rg) -> None:
    """Mean Rg +- block SEM per run; dotted Rg_theta for each construct's own length."""
    for lab, x in x_of.items():
        th = float(rg[f"g6-{lab}-285k"]["rg_theta_nm"])
        ax.hlines(th, x - 0.5, x + 0.5, color=INK, lw=1.2, ls=":", zorder=1)
        for T, c in temps:
            r = rg[f"g6-{lab}-{T}k"]
            ax.errorbar(x + off[T], float(r["rg_nm"]), yerr=float(r["rg_sem_nm"]), fmt="o",
                        ms=8, color=c, elinewidth=1.4, capsize=2.5, mew=0, zorder=3)
    ax.set_ylim(0, None)
    ax.set_ylabel("Radius of gyration (nm)", fontsize=10, color=INK)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("which", choices=["lysramp", "rgd"])
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if args.which == "lysramp":
        lysramp(plt)
    else:
        rgd(plt, args.plot_only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
