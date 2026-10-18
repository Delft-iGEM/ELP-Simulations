"""Mean distance between the 4 RGD sites of a chain, for all 24 RGD arrangements.

    .venv/bin/python tools/rgd_distance.py               # analyse + plot
    .venv/bin/python tools/rgd_distance.py --plot-only   # re-plot from the CSV

For every chain and every analysed frame: unwrap the chain along its bonds,
take each RGD's centre (mean of its R, G and D beads), and average the 6
pairwise distances between the 4 centres. That is averaged over frames (and,
in a brush, over the 16 chains); the error is the 20-block SEM over frames.
The nearest-neighbour distance (each RGD to its closest other RGD, averaged)
is written to the CSV too.

Two sets of the same 24 sequences, both 293.15 K, no TEA:

  free     free-rgdarr<NN>-*  one free chain, 200 ns, 35 nm box (filled dots)
  grafted  rgdarr<NN>-*       16 chains, preattached, 0.04 chains/nm^2, 400 ns
                              (open dots)

Frames before 50 ns are discarded; every 8th frame after that is used, the
same sampling as figures/free-rgdarr/sasa-24-arrangements.csv. The free runs
sit in a 35 nm box that is too small for a 352-aa chain (it reaches its own
periodic image), which tends to compact it; compare arrangements within a set,
not absolute values between sets.

Writes figures/free-rgdarr/rgd-distance-24-arrangements.{csv,png}.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

DISCARD_NS = 50.0
STRIDE = 8
OUT_DIR = Path("figures/free-rgdarr")
SASA_CSV = OUT_DIR / "sasa-24-arrangements.csv"
INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
END_COLOUR = {0: "#2a6fd1", 1: "#e8642c", 2: "#1aaa7a"}
END_LABEL = {0: "RGD at neither chain end", 1: "RGD at one chain end",
             2: "RGD at both chain ends"}
ORIGINAL = "01-orig"


def arrangements() -> list[dict]:
    return [{"label": r["run"], "copies_at_chain_end": int(r["copies_at_chain_end"]),
             "motif_23mer_nm2": float(r["motif_23mer_nm2"])}
            for r in csv.DictReader(open(SASA_CSV))]


def analyse(name: str) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.tea_analysis import block_average
    from tools.tea_g6_analysis import _whole

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
    n_chain = traj.n_atoms // len(seq)
    starts = [m.start() for m in re.finditer("RGD", seq)]
    pairs = list(itertools.combinations(range(len(starts)), 2))

    mean_d, nn_d = [], []
    for frame in traj.xyz:
        m_f, n_f = [], []
        for c in range(n_chain):
            pos = _whole(frame[c * len(seq):(c + 1) * len(seq)], box)
            centres = np.array([pos[s:s + 3].mean(axis=0) for s in starts])
            d = np.linalg.norm(centres[:, None] - centres[None], axis=2)
            m_f.append(np.mean([d[i, j] for i, j in pairs]))
            np.fill_diagonal(d, np.inf)
            n_f.append(d.min(axis=1).mean())
        mean_d.append(np.mean(m_f))
        nn_d.append(np.mean(n_f))
    out = {"run": name, "n_chains": n_chain, "frames": len(traj)}
    out["mean_nm"], out["mean_sem_nm"] = block_average(np.array(mean_d))
    out["nn_nm"], out["nn_sem_nm"] = block_average(np.array(nn_d))
    return out


def plot(rows: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = sorted({r["label"] for r in rows if r["set"] == "free"},
                    key=lambda lab: next(r["mean_nm"] for r in rows
                                         if r["label"] == lab and r["set"] == "free"))
    meta = {a["label"]: a for a in arrangements()}
    by_23 = sorted(meta.values(), key=lambda a: a["motif_23mer_nm2"])
    role = {ORIGINAL: "original", by_23[-1]["label"]: "winner", by_23[0]["label"]: "loser"}
    get = {(r["label"], r["set"]): r for r in rows}

    fig, ax = plt.subplots(figsize=(8.4, 9.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.tick_params(axis="y", length=0)
    ax.grid(color=GRID, linewidth=0.7, axis="x")
    ax.set_axisbelow(True)

    for i, lab in enumerate(labels):
        c = END_COLOUR[meta[lab]["copies_at_chain_end"]]
        bold = lab in role
        if bold:
            ax.axhspan(i - 0.42, i + 0.42, color=GRID, alpha=0.45, lw=0, zorder=0)
        f, g = get.get((lab, "free")), get.get((lab, "grafted"))
        if f and g:
            ax.plot([f["mean_nm"], g["mean_nm"]], [i, i], color=c, lw=1, alpha=0.35, zorder=1)
        # grafted ring drawn first and a little larger, so a free dot at the same
        # value stays visible inside it
        for r, filled in ((g, False), (f, True)):
            if not r:
                continue
            ax.errorbar(r["mean_nm"], i, xerr=r["mean_sem_nm"], fmt="o", color=c,
                        ms=(8.5 if bold else 6.0) if filled else (11 if bold else 9.5),
                        mfc=c if filled else SURFACE,
                        mec=INK if bold else (SURFACE if filled else c),
                        mew=1.6 if bold or not filled else 1.0,
                        elinewidth=2 if bold else 1.3, capsize=2.5, zorder=4 if filled else 3)
    for s, ls in (("free", (0, (4, 3))), ("grafted", (0, (1, 2)))):
        v = [r["mean_nm"] for r in rows if r["set"] == s]
        if v:
            ax.axvline(np.mean(v), color=MUTED, lw=1, ls=ls, zorder=0)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([f"{lab}  ({role[lab]})" if lab in role else lab for lab in labels])
    for t, lab in zip(ax.get_yticklabels(), labels):
        t.set_fontweight("bold" if lab in role else "normal")
        t.set_color(INK if lab in role else MUTED)
    ax.set_ylim(-0.7, len(labels) - 0.3)
    ax.set_xlabel("Mean distance between the RGD sites of a chain (nm)", fontsize=10, color=INK)

    handles = [plt.Line2D([], [], marker="o", ls="", color=c, mec=SURFACE, ms=7.5,
                          label=END_LABEL[n]) for n, c in END_COLOUR.items()]
    handles += [plt.Line2D([], [], marker="o", ls="", color=MUTED, mec=SURFACE, ms=7.5,
                           label="single free chain"),
                plt.Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=MUTED, mew=1.6, ms=7.5,
                           label="grafted, 16 chains")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, 1.01), ncol=3,
              frameon=False, fontsize=8.6, labelcolor=INK, columnspacing=1.2)
    fig.suptitle("Distance between RGD sites across all 24 RGD arrangements", x=0.02, y=0.985,
                 ha="left", fontsize=13, color=INK, fontweight="bold")
    fig.text(0.02, 0.952, "Mean of the 6 pairwise distances between the 4 RGD centres "
             "of one chain", fontsize=9.5, color=MUTED, ha="left")
    fig.text(0.02, 0.012,
             "RGD centre = mean of its R, G, D beads; chains unwrapped along their bonds. Free: "
             "one chain, 200 ns, 35 nm box. Grafted: 16 chains,\npreattached, 0.04 chains/nm², "
             "400 ns, averaged over chains. Both 293 K, no TEA, every 8th frame from 50 ns; "
             "bar = 20-block SEM.\nDashed line = mean of the free set, dotted = grafted. Rows "
             "sorted by the free-chain distance;\nbold: the original, and the winner and loser "
             "by RGD + linker SASA.", fontsize=8.0, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.09, 1, 0.905))
    fig.subplots_adjust(top=0.865)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args(argv)
    out_csv = OUT_DIR / "rgd-distance-24-arrangements.csv"

    if args.plot_only:
        rows = [{k: (v if k in ("run", "label", "set") else float(v)) for k, v in r.items()}
                for r in csv.DictReader(open(out_csv))]
    else:
        cells = [(a["label"], s, f"{p}rgdarr{a['label']}")
                 for a in arrangements() for s, p in (("free", "free-"), ("grafted", ""))]
        with Pool(args.jobs) as pool:
            results = pool.map(analyse, [n for _, _, n in cells])
        rows = [{"label": lab, "set": s} | r for (lab, s, _), r in zip(cells, results) if r]
        missing = [n for (_, _, n), r in zip(cells, results) if not r]
        if missing:
            print(f"missing ({len(missing)}): {', '.join(missing)}")
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out_csv}")
        for r in sorted(rows, key=lambda r: (r["set"], r["mean_nm"])):
            print(f"  {r['set']:<8}{r['label']:<18}{r['mean_nm']:6.2f} ±{r['mean_sem_nm']:.2f}"
                  f"   nn {r['nn_nm']:5.2f}")
    plot(rows, OUT_DIR / "rgd-distance-24-arrangements.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
