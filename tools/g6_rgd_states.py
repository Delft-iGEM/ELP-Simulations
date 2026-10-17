"""SASA in the collapsed and in the extended state, for the g6 RGD-placement runs.

    .venv/bin/python tools/g6_rgd_states.py               # analyse + plot
    .venv/bin/python tools/g6_rgd_states.py --plot-only   # re-plot from the CSV

Runs: g6-orig and the eight g6-rgd* arrangements (the best 4 and worst 4 of
the free-chain RGD sweep, tools/tea_g6.py), each at 285 and 370 K: single free
chains, CALVADOS-TEA gamma 6, 400 ns, 50-400 ns analysed.

Every frame is called collapsed (Rg below the theta-state Rg, 0.21 N^0.5 nm,
as in tools/tea_g6_states.py) or extended, and three SASAs are averaged
separately over the two sets of frames:

  rgd     the 4 RGD tripeptides, summed            every frame
  23mer   the 4 LDASTVYAVTGRGDSPASSAASA domains     every frame
  chain   the whole chain                          every 5th frame

Errors: 20-block SEM over the blocks that contain frames in that state
(tools/tea_g6_states._state_mean); fewer than 3 such blocks gives no error.
At 285 K a chain is collapsed only ~1.5 % of the time, so the 285 K collapsed
means rest on a few short visits.

g6-orig is 357 aa, the rgd* arrangements 352 aa (tools/tea_g6.py); the domain
and RGD SASAs compare directly, the whole-chain SASA carries the 5 extra
residues of the original.

Writes figures/g6/g6-rgd-states.csv and g6-rgd-states-{rgd,23mer,chain}.png:
two panels (285 K, 370 K), collapsed and extended side by side per row, rows
ordered by the free-chain RGD + linker ranking, a chain sketch per row, and the
original, the free-chain winner (19-3n-1c) and loser (16-pairs-mid) highlighted.
Also g6-rgd-states-{rgd,23mer,chain}-combined.png: both temperatures on one axis,
colour = temperature, filled circle = extended, open square = collapsed.
"""

from __future__ import annotations

import argparse
import csv
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_analysis import rg_theta
from tools.tea_g6_analysis import DISCARD_NS, DOMAIN, _whole
from tools.tea_g6_states import _state_mean

TEMPS = (285, 370)
CHAIN_STRIDE = 5
OUT_DIR = Path("figures/g6")
OUT_CSV = OUT_DIR / "g6-rgd-states.csv"
FREE_CSV = Path("figures/free-rgdarr/sasa-24-arrangements.csv")
# sequence label -> free-chain arrangement label (for the ranking); orig is 01-orig
# in free-chain ranking order (most exposed first)
LABELS = ["rgd19-3n-1c", "rgd15-pairs-ends", "rgd14-tight-c", "rgd22-c-heavy", "orig",
          "rgd23-staircase", "rgd03-even-inset", "rgd10-clust-q1", "rgd16-pairs-mid"]
ROLE = {"orig": "original", "rgd19-3n-1c": "free winner", "rgd16-pairs-mid": "free loser"}
QUANT = [("rgd", "RGD SASA (nm²)"), ("23mer", "RGD + linker SASA (nm²)"),
         ("chain", "Whole-chain SASA (nm²)")]
INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
STATE = {"extended": ("#2a6fd1", "o"), "collapsed": ("#e8642c", "s")}
CHAIN_FILL, CHAIN_EDGE, RGD_DOT = "#f5c518", "#b8920c", "#d7263d"
TEMP_COLD, TEMP_HOT = "#4269d0", "#c1442e"            # 285 / 370 K, as in the dwell figures


def analyse(name: str) -> dict:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads

    runtime = Path("simulations") / name / "runtime"
    top, dcd = trajectory_paths(name)
    traj = md.load(str(dcd), top=str(top))
    meta = read_metadata(runtime)
    t = np.arange(traj.n_frames) * float(meta["time_per_frame_ps"]) / 1000.0
    traj = traj[t >= DISCARD_NS]
    box = traj.unitcell_lengths[0]
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    n = len(seq)

    xyz = np.stack([_whole(f, box) for f in traj.xyz])
    rg = np.sqrt(((xyz - xyz.mean(axis=1, keepdims=True)) ** 2).sum(axis=2).mean(axis=1))
    collapsed = rg < rg_theta(n)

    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    dom_starts = [m.start() for m in re.finditer(DOMAIN, seq)]
    dom = np.array([i for s in dom_starts for i in range(s, s + len(DOMAIN))])
    in_rgd = np.isin(dom, [m.start() + k for m in re.finditer("RGD", seq) for k in range(3)])
    s_dom, s_rgd = np.empty(traj.n_frames), np.empty(traj.n_frames)
    for f in range(traj.n_frames):
        s = sasa_of_beads(traj.xyz[f], radii, which=dom, box=box)
        s_dom[f], s_rgd[f] = s.sum(), s[in_rgd].sum()
    every = np.arange(0, traj.n_frames, CHAIN_STRIDE)
    s_chain = np.array([sasa_of_beads(traj.xyz[f], radii, box=box).sum() for f in every])

    row = {"run": name, "label": name[3:].rsplit("-", 1)[0],
           "temperature_K": int(round(float(meta["temperature_K"]))), "n_residues": n,
           "domain_starts": " ".join(str(s) for s in dom_starts),
           "collapse_fraction": float(collapsed.mean()),
           "collapsed_frames": int(collapsed.sum()), "frames": int(len(collapsed))}
    for q, x, mask in (("rgd", s_rgd, collapsed), ("23mer", s_dom, collapsed),
                       ("chain", s_chain, collapsed[every])):
        for state, m in (("collapsed", mask), ("extended", ~mask)):
            row[f"{q}_{state}_nm2"], row[f"{q}_{state}_sem_nm2"] = _state_mean(x, m)
    print(f"  {name}: collapsed {row['collapse_fraction']:.1%}  23mer "
          f"{row['23mer_collapsed_nm2']:.1f} / {row['23mer_extended_nm2']:.1f}", flush=True)
    return row


def plot(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = LABELS[::-1]                                   # free-chain best at the top
    y = {lab: i for i, lab in enumerate(labels)}
    by = {(r["label"], int(r["temperature_K"])): r for r in rows}
    for q, xlabel in QUANT:
        fig, axes = plt.subplots(1, 3, figsize=(12.6, 6.4), sharey=True,
                                 gridspec_kw={"width_ratios": (3, 3, 1.4), "wspace": 0.08})
        for ax in axes:
            ax.set_facecolor(SURFACE)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            ax.tick_params(axis="y", length=0)
            ax.tick_params(colors=INK, labelsize=9)
            for lab in ROLE:
                ax.axhspan(y[lab] - 0.42, y[lab] + 0.42, color=GRID, alpha=0.45, lw=0, zorder=0)
        for ax, T in zip(axes[:2], TEMPS):
            ax.spines["bottom"].set_color(INK)
            ax.grid(color=GRID, linewidth=0.7, axis="x")
            ax.set_axisbelow(True)
            for lab in labels:
                r = by.get((lab, T))
                if not r:
                    continue
                for k, (state, (c, m)) in enumerate(STATE.items()):
                    v, e = float(r[f"{q}_{state}_nm2"]), float(r[f"{q}_{state}_sem_nm2"])
                    bold = lab in ROLE
                    ax.errorbar(v, y[lab] + (0.14 if k == 0 else -0.14),
                                xerr=None if np.isnan(e) else e, fmt=m, color=c,
                                ms=7.5 if bold else 6, elinewidth=1.4, capsize=2.5,
                                capthick=1.4, mec=INK if bold else c, mew=1.2 if bold else 0,
                                zorder=3)
            ax.set_title(f"{T} K", loc="left", fontsize=11, color=INK, fontweight="bold")
            ax.set_xlabel(xlabel, fontsize=10, color=INK)
        ax0 = axes[0]
        ax0.set_yticks(range(len(labels)))
        names = [lab.replace("rgd", "", 1) if lab != "orig" else "01-orig" for lab in labels]
        free = sorted(csv.DictReader(open(FREE_CSV)), key=lambda r: -float(r["motif_23mer_nm2"]))
        place = {r["run"]: i + 1 for i, r in enumerate(free)}        # 1 = most exposed of 24
        ticks = []
        for n, lab in zip(names, labels):
            k = place[n]
            if lab == "orig":
                sfx = "th" if 10 <= k % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
                ticks.append(f"{n}  (original, {k}{sfx} place)")
            else:
                ticks.append(f"{n}  ({ROLE[lab]})" if lab in ROLE else f"{n}  (placed #{k})")
        ax0.set_yticklabels(ticks)
        for t, lab in zip(ax0.get_yticklabels(), labels):
            t.set_color(INK)
            if lab in ROLE:
                t.set_fontweight("bold")
        ax0.set_ylim(-0.6, len(labels) - 0.4)

        sx = axes[2]
        sx.tick_params(bottom=False, labelbottom=False)
        sx.spines["bottom"].set_visible(False)
        for lab in labels:
            r = by.get((lab, TEMPS[0])) or by.get((lab, TEMPS[1]))
            n = int(r["n_residues"])
            sx.plot([0, 1], [y[lab]] * 2, color=CHAIN_EDGE, lw=9.5, solid_capstyle="round",
                    zorder=1)
            sx.plot([0, 1], [y[lab]] * 2, color=CHAIN_FILL, lw=8, solid_capstyle="round",
                    zorder=2)
            g = [(int(s) + 12) / (n - 1) for s in r["domain_starts"].split()]   # the G of RGD
            sx.plot(g, [y[lab]] * len(g), "o", color=RGD_DOT, ms=6, mec=INK, mew=0.6, zorder=3)
        sx.set_xlim(-0.06, 1.06)
        for x, t in ((0, "N"), (1, "C")):
            sx.text(x, len(labels) - 0.55, t, ha="center", va="bottom", fontsize=9,
                    color=INK, fontweight="bold")
        sx.set_xlabel("RGD positions", fontsize=10, color=INK)

        handles = [plt.Line2D([], [], marker=m, ls="", color=c, ms=7,
                              label=f"{state} (Rg {'<' if state == 'collapsed' else '>'} Rg_θ)")
                   for state, (c, m) in STATE.items()]
        handles.append(plt.Line2D([], [], marker="o", ls="", color=RGD_DOT, mec=INK, mew=0.6,
                                  ms=6, label="RGD position (yellow bar = chain)"))
        fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9,
                   labelcolor=INK, bbox_to_anchor=(0.5, 0.995))
        fig.subplots_adjust(left=0.21, right=0.985, bottom=0.09, top=0.88)
        out = OUT_DIR / f"g6-rgd-states-{q}.png"
        fig.savefig(out, dpi=200, transparent=True)
        plt.close(fig)
        print(f"wrote {out}")


def plot_combined(rows: list[dict]) -> None:
    """Both temperatures on one SASA axis: per row a 285 K and a 370 K dumbbell,
    filled circle = extended, open square = collapsed, colour = temperature.
    Writes g6-rgd-states-<q>-combined.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = LABELS[::-1]
    y = {lab: i for i, lab in enumerate(labels)}
    by = {(r["label"], int(r["temperature_K"])): r for r in rows}
    free = sorted(csv.DictReader(open(FREE_CSV)), key=lambda r: -float(r["motif_23mer_nm2"]))
    place = {r["run"]: i + 1 for i, r in enumerate(free)}
    temp_colour = {285: TEMP_COLD, 370: TEMP_HOT}
    dy = {285: 0.17, 370: -0.17}
    for q, xlabel in QUANT:
        fig, (ax, sx) = plt.subplots(1, 2, figsize=(10.4, 6.6), sharey=True,
                                     gridspec_kw={"width_ratios": (4.2, 1.3), "wspace": 0.05})
        for a in (ax, sx):
            for side in ("top", "right", "left"):
                a.spines[side].set_visible(False)
            a.tick_params(axis="y", length=0)
            a.tick_params(colors=INK, labelsize=9)
            for lab in ROLE:
                a.axhspan(y[lab] - 0.42, y[lab] + 0.42, color=GRID, alpha=0.45, lw=0, zorder=0)
        ax.spines["bottom"].set_color(INK)
        ax.grid(color=GRID, linewidth=0.7, axis="x")
        ax.set_axisbelow(True)
        for lab in labels:
            bold = lab in ROLE
            for T in TEMPS:
                r = by.get((lab, T))
                if not r:
                    continue
                c, yy = temp_colour[T], y[lab] + dy[T]
                ext, col = float(r[f"{q}_extended_nm2"]), float(r[f"{q}_collapsed_nm2"])
                ax.plot([col, ext], [yy, yy], color=c, lw=1.6 if bold else 1.1, alpha=0.6,
                        zorder=2)
                for v, state in ((ext, "extended"), (col, "collapsed")):
                    e = float(r[f"{q}_{state}_sem_nm2"])
                    filled = state == "extended"
                    ax.errorbar(v, yy, xerr=None if np.isnan(e) else e,
                                fmt="o" if filled else "s", color=c, mfc=c if filled else "white",
                                mec=c, mew=1.8, ms=7.5 if bold else 6.5, elinewidth=1.3,
                                capsize=2.5, capthick=1.3, zorder=3)
        names = [lab.replace("rgd", "", 1) if lab != "orig" else "01-orig" for lab in labels]
        ticks = []
        for n, lab in zip(names, labels):
            k = place[n]
            if lab == "orig":
                sfx = "th" if 10 <= k % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
                ticks.append(f"{n}  (original, {k}{sfx} place)")
            else:
                ticks.append(f"{n}  ({ROLE[lab]})" if lab in ROLE else f"{n}  (placed #{k})")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(ticks)
        for t, lab in zip(ax.get_yticklabels(), labels):
            t.set_color(INK)
            if lab in ROLE:
                t.set_fontweight("bold")
        ax.set_ylim(-0.6, len(labels) - 0.4)
        ax.set_xlabel(xlabel, fontsize=10, color=INK)

        sx.tick_params(bottom=False, labelbottom=False)
        sx.spines["bottom"].set_visible(False)
        for lab in labels:
            r = by.get((lab, TEMPS[0])) or by.get((lab, TEMPS[1]))
            n = int(r["n_residues"])
            sx.plot([0, 1], [y[lab]] * 2, color=CHAIN_EDGE, lw=9.5, solid_capstyle="round",
                    zorder=1)
            sx.plot([0, 1], [y[lab]] * 2, color=CHAIN_FILL, lw=8, solid_capstyle="round",
                    zorder=2)
            g = [(int(s) + 12) / (n - 1) for s in r["domain_starts"].split()]   # the G of RGD
            sx.plot(g, [y[lab]] * len(g), "o", color=RGD_DOT, ms=6, mec=INK, mew=0.6, zorder=3)
        sx.set_xlim(-0.06, 1.06)
        sx.set_xlabel("RGD positions", fontsize=10, color=INK)

        handles = [plt.Line2D([], [], color=temp_colour[T], lw=3, label=f"{T} K")
                   for T in TEMPS]
        handles += [plt.Line2D([], [], marker="o", ls="", color=MUTED, mec=MUTED, mew=1.8,
                               ms=7, label="extended (Rg > Rg_θ)"),
                    plt.Line2D([], [], marker="s", ls="", mfc="white", mec=MUTED, mew=1.8,
                               ms=7, label="collapsed (Rg < Rg_θ)"),
                    plt.Line2D([], [], marker="o", ls="", color=RGD_DOT, mec=INK, mew=0.6,
                               ms=6, label="RGD position (yellow bar = chain)")]
        fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9,
                   labelcolor=INK, bbox_to_anchor=(0.5, 0.995))
        fig.subplots_adjust(left=0.27, right=0.985, bottom=0.09, top=0.88)
        out = OUT_DIR / f"g6-rgd-states-{q}-combined.png"
        fig.savefig(out, dpi=200, transparent=True)
        plt.close(fig)
        print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args(argv)
    if not args.plot_only:
        names = [f"g6-{lab}-{T}k" for lab in LABELS for T in TEMPS]
        with Pool(len(names)) as pool:
            rows = pool.map(analyse, names)
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {OUT_CSV}")
    rows = list(csv.DictReader(open(OUT_CSV)))
    plot(rows)
    plot_combined(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
