"""Solvent exposure of the RGD construct across the g6-tsw temperature sweep.

    .venv/bin/python tools/tea_tsweep_sasa.py               # analyse + plot
    .venv/bin/python tools/tea_tsweep_sasa.py --plot-only   # re-plot from the CSVs

Runs: g6-tsw-275k .. 380k (tools/tea_tsweep.py): one free 357-aa chain, TEA
gamma 6, 400 ns, the same seed at every T; 275 K is TEA extrapolated below its
280 K fit range. Unfinished runs are skipped.

Coarse-grained Shrake-Rupley SASA (tools/sasa.py: bead radius sigma/2, probe
0.14 nm, minimum image). Every bead's SASA is computed on every 10th frame from
50 ns on (~500 frames, 0.7 ns apart), and grouped:

  chain         the whole chain
  rgd_site      the R, G and D beads of one RGD             (mean of the 4)
  domain_site   one 23-mer cell-binding domain              (mean of the 4)
  lys_site      one lysine bead, the crosslink site         (mean of the 4)
  elp           the ELP blocks (everything outside the 4 domains), per residue

`*_exposed` is the same SASA divided by what those beads would have with
nothing around them (4 pi (r + probe)^2 each): 1 = fully exposed. Error bars
are the SEM from 20 block averages of the per-frame values.

Writes, in figures/g6-tsw/:
  g6-tsw-sasa.csv           one row per temperature, every group, nm^2 and exposed
  g6-tsw-sasa-residues.csv  exposed fraction of every residue at every temperature
  g6-tsw-sasa-{chain,rgd,domain,lys,elp}.png   one quantity each, vs T
  g6-tsw-sasa-exposed.png   the five groups as exposed fraction, one shared scale
  g6-tsw-sasa-profile.png   exposed fraction along the sequence (residue x T)
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_tsweep import TEMPERATURES

STRIDE = 10
DISCARD_NS = 50.0
DOMAIN = "LDASTVYAVTGRGDSPASSAASA"
OUT_DIR = Path("figures/g6-tsw")
FIT_LO = 280                                     # TEA fit range starts here
RUNS = [(int(T), f"g6-tsw-{int(T)}k") for T in TEMPERATURES]
GROUPS = ["chain", "rgd_site", "domain_site", "lys_site", "elp"]


def groups_of(seq: str) -> dict[str, list[list[int]]]:
    dom = [list(range(m.start(), m.start() + len(DOMAIN))) for m in re.finditer(DOMAIN, seq)]
    in_dom = {i for d in dom for i in d}
    return {"chain": [list(range(len(seq)))],
            "rgd_site": [list(range(m.start(), m.start() + 3)) for m in re.finditer("RGD", seq)],
            "domain_site": dom,
            "lys_site": [[i] for i, c in enumerate(seq) if c == "K"],
            # per residue: one "site" per ELP bead
            "elp": [[i] for i in range(len(seq)) if i not in in_dom]}


def analyse(item: tuple[int, str]) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import PROBE_NM, bead_radii, sasa_of_beads
    from tools.tea_analysis import block_average

    T, name = item
    try:
        top, dcd = trajectory_paths(name)
    except FileNotFoundError:
        return None
    if not Path(dcd).is_file():
        return None
    runtime = Path("simulations") / name / "runtime"
    traj = md.load(str(dcd), top=str(top))
    ps = float(read_metadata(runtime)["time_per_frame_ps"])
    traj = traj[np.arange(traj.n_frames) * ps / 1000.0 >= DISCARD_NS][::STRIDE]
    box = traj.unitcell_lengths[0]
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    if traj.n_atoms != len(seq):
        raise ValueError(f"{name}: {traj.n_atoms} beads, {len(seq)} residues")
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    alone = 4.0 * np.pi * (radii + PROBE_NM) ** 2

    per_bead = np.array([sasa_of_beads(f, radii, box=box) for f in traj.xyz])   # frames x beads
    row = {"T": T, "run": name, "n_residues": len(seq), "frames": len(traj),
           "extrapolated": int(T < FIT_LO)}
    for g, sites in groups_of(seq).items():
        v = np.mean([per_bead[:, s].sum(axis=1) for s in sites], axis=0)      # per site, per frame
        ref = np.mean([alone[s].sum() for s in sites])
        row[f"{g}_nm2"], row[f"{g}_sem_nm2"] = block_average(v)
        row[f"{g}_exposed"], row[f"{g}_exposed_sem"] = block_average(v / ref)
    row["chain_per_residue_nm2"] = row["chain_nm2"] / len(seq)
    row["_residue_exposed"] = (per_bead.mean(axis=0) / alone).tolist()
    row["_seq"] = seq
    return row


# ---- figures ---------------------------------------------------------------------

INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
LINE = "#5a5f66"                                 # one series: neutral ink, T on x
QUANT = [  # tag, group, y label, title
    ("chain", "chain", "SASA per chain (nm²)", "Whole-chain SASA"),
    ("rgd", "rgd_site", "SASA per RGD (nm²)", "RGD exposure"),
    ("domain", "domain_site", "SASA per 23-mer domain (nm²)",
     "Cell-binding domain (RGD + linker) exposure"),
    ("lys", "lys_site", "SASA per lysine (nm²)", "Lysine (crosslink site) exposure"),
    ("elp", "elp", "SASA per ELP residue (nm²)", "ELP-block exposure"),
]
# five groups on one exposed-fraction axis: categorical, fixed order
GROUP_COLOUR = {"chain": "#5a5f66", "elp": "#2a6fd1", "domain_site": "#1aaa7a",
                "rgd_site": "#e8642c", "lys_site": "#8a5cc2"}
GROUP_NAME = {"chain": "whole chain", "elp": "ELP blocks", "domain_site": "23-mer domain",
              "rgd_site": "RGD", "lys_site": "lysine"}


def _axes(plt):
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("Temperature (K)", fontsize=9.5, color=INK)
    ax.set_xticks(range(280, 381, 10))
    ax.axvspan(272, FIT_LO - 2.5, color=GRID, alpha=0.5, lw=0, zorder=0)
    return fig, ax


def _marks(ax, m: dict) -> None:
    for key, label in (("T_theta_K", "T$_θ$"), ("T_half_K", "T$_{½}$")):
        if m.get(key) and np.isfinite(m[key]):
            ax.axvline(m[key], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
            ax.text(m[key], 1.0, f" {label} {m[key]:.0f} K", transform=ax.get_xaxis_transform(),
                    fontsize=8, color=MUTED, va="top", ha="left" if key == "T_theta_K" else "right")


def _finish(fig, plt, ax, title: str, caption: str, out: Path) -> None:
    ax.set_xlim(272, 383)
    fig.suptitle(title, x=0.02, y=0.98, ha="left", fontsize=12.5, color=INK, fontweight="bold")
    fig.text(0.02, 0.012, caption, fontsize=7.8, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.16, 1, 0.94))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


CAPTION = ("One free 357-aa chain per temperature (g6-tsw, TEA γ 6, 400 ns, same seed); every "
           "10th frame from 50 ns. Bar = 20-block SEM.\nGrey band: below TEA's 280 K fit range "
           "(extrapolated). Dashed: T$_θ$ (mean Rg = Rg$_θ$) and T$_{½}$ (half the frames "
           "collapsed),\nfrom g6-tsw.json.")


def plot(rows: list[dict], res: dict, m: dict, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T = np.array([r["T"] for r in rows])
    for tag, g, ylab, title in QUANT:
        fig, ax = _axes(plt)
        y = np.array([r[f"{g}_nm2"] for r in rows])
        e = np.array([r[f"{g}_sem_nm2"] for r in rows])
        ax.plot(T, y, color=LINE, lw=2, zorder=2)
        ax.errorbar(T, y, yerr=e, fmt="o", color=LINE, ms=6, mec=SURFACE, mew=1.2,
                    elinewidth=1.3, capsize=2.5, zorder=3)
        _marks(ax, m)
        ax.set_ylabel(ylab, fontsize=9.5, color=INK)
        drop = 100 * (1 - y[-1] / y[T == 290][0]) if 290 in T else float("nan")
        ax2 = ax.secondary_yaxis("right", functions=(
            lambda v, ref=rows[0][f"{g}_nm2"] / rows[0][f"{g}_exposed"]: v / ref,
            lambda f, ref=rows[0][f"{g}_nm2"] / rows[0][f"{g}_exposed"]: f * ref))
        ax2.set_ylabel("exposed fraction", fontsize=9, color=MUTED)
        ax2.tick_params(colors=MUTED, labelsize=8.5)
        ax2.spines["right"].set_color(GRID)
        _finish(fig, plt, ax, title,
                CAPTION + f"\nRight axis: the same value as a fraction of the beads' SASA with "
                f"nothing around them. 290 → {int(T[-1])} K: {drop:+.0f}% buried.".replace("+", "−")
                if drop > 0 else CAPTION, out_dir / f"g6-tsw-sasa-{tag}.png")

    # all groups on one scale: exposed fraction
    fig, ax = _axes(plt)
    ends = {}
    for g in ("lys_site", "rgd_site", "domain_site", "elp", "chain"):
        y = np.array([r[f"{g}_exposed"] for r in rows])
        e = np.array([r[f"{g}_exposed_sem"] for r in rows])
        c = GROUP_COLOUR[g]
        ax.plot(T, y, color=c, lw=2, zorder=2)
        ax.errorbar(T, y, yerr=e, fmt="o", color=c, ms=5.5, mec=SURFACE, mew=1.2,
                    elinewidth=1.2, capsize=2, zorder=3, label=GROUP_NAME[g])
        ends[g] = y[-1]
    # direct labels at the right end, pushed apart where lines end close together
    gap = 0.011
    placed = []
    for g, v in sorted(ends.items(), key=lambda kv: kv[1]):
        y_lab = max(v, placed[-1] + gap) if placed else v
        placed.append(y_lab)
        ax.text(T[-1] + 1.5, y_lab, GROUP_NAME[g], fontsize=8, color=INK, va="center")
    _marks(ax, m)
    ax.set_ylabel("Exposed fraction (SASA / SASA alone)", fontsize=9.5, color=INK)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.07), frameon=False, fontsize=8.3,
              labelcolor=INK, ncol=5, handletextpad=0.3, columnspacing=1.2)
    ax.set_xlim(272, 395)
    fig.suptitle("Which parts of the chain get buried on heating", x=0.02, y=0.98, ha="left",
                 fontsize=12.5, color=INK, fontweight="bold")
    fig.text(0.02, 0.012, CAPTION + "\nExposed fraction = SASA of the group divided by its "
             "beads' SASA with nothing around them (1 = fully exposed).",
             fontsize=7.8, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.17, 1, 0.9))
    fig.savefig(out_dir / "g6-tsw-sasa-exposed.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_dir / 'g6-tsw-sasa-exposed.png'}")

    # profile along the sequence
    seq = res["seq"]
    temps = sorted(int(t) for t in res["exposed"])
    from matplotlib.colors import TwoSlopeNorm
    z = np.array([res["exposed"][str(t)] for t in temps])
    z = z - z[0]                                       # change since the coldest run
    k = np.ones(5) / 5                                 # one pentad, removes bead-size stripes
    z = np.array([np.convolve(row, k, mode="same") for row in z])
    lim = np.abs(z).max()
    fig, ax = plt.subplots(figsize=(10.5, 4.9))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    im = ax.imshow(z, aspect="auto", origin="lower", cmap="RdBu",
                   norm=TwoSlopeNorm(0.0, -lim, lim),
                   extent=(-0.5, len(seq) - 0.5, temps[0] - 2.5, temps[-1] + 2.5))
    for mt in re.finditer(DOMAIN, seq):
        ax.axvspan(mt.start() - 0.5, mt.end() - 0.5, ymin=1.0, ymax=1.04, color="#1aaa7a",
                   clip_on=False, lw=0)
        r = mt.start() + DOMAIN.index("RGD")
        ax.axvspan(r - 0.5, r + 2.5, ymin=1.0, ymax=1.04, color="#e8642c", clip_on=False, lw=0)
    for i, c in enumerate(seq):
        if c == "K":
            ax.plot(i, temps[-1] + 2.5, marker="v", color="#8a5cc2", ms=5, clip_on=False,
                    transform=ax.transData, zorder=5)
    ax.set_xlabel("Residue", fontsize=9.5, color=INK)
    ax.set_ylabel("Temperature (K)", fontsize=9.5, color=INK)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label(f"Change in exposed fraction since {temps[0]} K", fontsize=9, color=INK)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("Where along the chain exposure is lost on heating", x=0.02, y=0.98, ha="left",
                 fontsize=12.5, color=INK, fontweight="bold")
    fig.text(0.02, 0.92, "Above the plot: green = 23-mer domain, orange = its RGD, "
             "purple ▼ = lysine. Red = more buried than at " + f"{temps[0]} K.", fontsize=8.5, color=MUTED)
    fig.text(0.02, 0.012, "Exposed fraction of each bead (mean SASA over every 10th frame from 50 ns / SASA "
             f"with nothing around it), minus the same bead's value at {temps[0]} K,\nthen a "
             "5-residue running mean. One free chain per temperature, TEA γ 6, same seed; "
             "rows below 280 K are TEA extrapolated. Per-residue values are in g6-tsw-sasa-residues.csv.", fontsize=7.8, color=MUTED,
             linespacing=1.5)
    fig.tight_layout(rect=(0, 0.08, 1, 0.9))
    fig.savefig(out_dir / "g6-tsw-sasa-profile.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out_dir / 'g6-tsw-sasa-profile.png'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=22)
    args = parser.parse_args(argv)
    out_csv, res_csv = OUT_DIR / "g6-tsw-sasa.csv", OUT_DIR / "g6-tsw-sasa-residues.csv"

    if args.plot_only:
        rows = [{k: (v if k == "run" else float(v)) for k, v in r.items()}
                for r in csv.DictReader(open(out_csv))]
        lines = list(csv.reader(open(res_csv)))
        seq = "".join(r[1] for r in lines[1:])
        res = {"seq": seq, "exposed": {h: [float(r[j]) for r in lines[1:]]
                                       for j, h in enumerate(lines[0]) if j >= 2}}
    else:
        with Pool(args.jobs) as pool:
            got = [r for r in pool.map(analyse, RUNS) if r]
        got.sort(key=lambda r: r["T"])
        seq = got[0]["_seq"]
        res = {"seq": seq, "exposed": {str(r["T"]): r.pop("_residue_exposed") for r in got}}
        rows = [{k: v for k, v in r.items() if k != "_seq"} for r in got]
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        with open(res_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["residue", "aa"] + list(res["exposed"]))
            for i, c in enumerate(seq):
                w.writerow([i, c] + [f"{res['exposed'][t][i]:.4f}" for t in res["exposed"]])
        print(f"wrote {out_csv}, {res_csv}  ({len(rows)} runs)")
        for r in rows:
            print(f"  {int(r['T'])} K  chain {r['chain_nm2']:6.1f}  RGD {r['rgd_site_nm2']:.3f}"
                  f"  domain {r['domain_site_nm2']:.2f}  K {r['lys_site_nm2']:.3f}"
                  f"  ELP/res {r['elp_nm2']:.3f}")
    m = json.loads((OUT_DIR / "g6-tsw.json").read_text())
    plot(rows, res, m, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
