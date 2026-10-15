"""Rg and SASA for the lysine-fraction series (tools/tea_lysfrac.py), and its figure.

    .venv/bin/python tools/tea_lysfrac_analysis.py                  # analyse + plot
    .venv/bin/python tools/tea_lysfrac_analysis.py --plot-only      # re-plot the CSV

Per run, over 50-400 ns, everything tools/tea_g6_analysis.py measures (Rg,
collapse fraction, extent, K-K contacts) plus two SASAs the g6 set did not need:

  chain_sasa   SASA of the whole chain, every 10th frame
  lys_sasa     SASA summed over the lysine beads, divided by the number of
               lysines, so 10 % and 100 % read on the same per-residue scale

These chains have no RGD, so the g6 script's RGD and 23-mer SASA columns are
dropped. Writes figures/g6-kfrac/g6-kfrac.{csv,png}.
"""

from __future__ import annotations

import argparse
import csv
import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np

FRACTIONS = list(range(0, 101, 10))
TEMPS = (285, 370)
SASA_STRIDE = 10
DISCARD_NS = 50.0


def run_name(pct: int, t: int) -> str:
    return f"g6-kfrac{pct:03d}-{t}k"


def analyse(name: str) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads
    from tools.tea_analysis import block_average
    from tools.tea_g6_analysis import analyse_run

    row = analyse_run(name)
    if row is None:
        return None
    for k in ("rgd_sasa_nm2", "rgd_sasa_sem_nm2", "domain_sasa_nm2", "domain_sasa_sem_nm2"):
        row.pop(k)

    top, dcd = trajectory_paths(name)
    traj = md.load(str(dcd), top=str(top))
    meta = read_metadata(Path("simulations") / name / "runtime")
    t = np.arange(traj.n_frames) * float(meta["time_per_frame_ps"]) / 1000.0
    traj = traj[t >= DISCARD_NS]
    box = traj.unitcell_lengths[0]
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    lys = np.array([a.index for a in traj.topology.atoms if a.residue.name == "LYS"])

    chain, per_lys = [], []
    for frame in traj.xyz[::SASA_STRIDE]:
        s = sasa_of_beads(frame, radii, box=box)
        chain.append(s.sum())
        if lys.size:
            per_lys.append(s[lys].mean())
    row["chain_sasa_nm2"], row["chain_sasa_sem_nm2"] = block_average(chain)
    row["lys_sasa_per_k_nm2"], row["lys_sasa_per_k_sem_nm2"] = (
        block_average(per_lys) if per_lys else (float("nan"), float("nan")))
    row["lys_pct"] = int(name.split("kfrac")[1][:3])
    return row


# ---- figure ------------------------------------------------------------------

COLD = "#4269d0"
HOT = "#c1442e"
INK = "#1a1a1a"
MUTED = "#6b6f76"
GRID = "#dcdee2"
SURFACE = "#fcfcfb"


def plot(rows: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [("rg_nm", "rg_sem_nm", "Radius of gyration (nm)"),
              ("chain_sasa_nm2", "chain_sasa_sem_nm2", "SASA, whole chain (nm²)"),
              ("lys_sasa_per_k_nm2", "lys_sasa_per_k_sem_nm2", "SASA per lysine (nm²)")]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    fig.patch.set_facecolor(SURFACE)
    for ax, (key, sem, label) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for t, c in ((285, COLD), (370, HOT)):
            pts = sorted((r["lys_pct"], r[key], r[sem]) for r in rows
                         if int(round(r["temperature_K"])) == t and not math.isnan(r[key]))
            if not pts:
                continue
            x, y, e = map(np.array, zip(*pts))
            ax.plot(x, y, color=c, linewidth=2, zorder=2)
            ax.errorbar(x, y, yerr=e, fmt="o", ms=7, color=c, ecolor=c, elinewidth=1,
                        capsize=0, mec=SURFACE, mew=1.5, zorder=3, label=f"{t} K")
            # direct label at the right end of each line
            ax.annotate(f"{t} K", (x[-1], y[-1]), xytext=(6, 0), textcoords="offset points",
                        va="center", fontsize=8.5, color=INK)
        if key == "rg_nm":
            theta = rows[0]["rg_theta_nm"]
            ax.axhline(theta, color=MUTED, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
            ax.text(1, theta, " θ-state Rg", transform=ax.get_yaxis_transform(),
                    va="bottom", ha="right", fontsize=8, color=MUTED)
        ax.set_xlabel("pentads that are VPGKG (%)", fontsize=9, color=MUTED)
        ax.set_title(label, loc="left", fontsize=10, color=INK, fontweight="bold")
        ax.set_xticks(FRACTIONS)
        ax.set_xlim(-4, 110)
    handles = [plt.Line2D([], [], marker="o", color=c, mec=SURFACE, ms=7, label=f"{t} K")
               for t, c in ((285, COLD), (370, HOT))]
    fig.legend(handles=handles, loc="upper right", ncol=2, frameon=False, fontsize=9,
               labelcolor=INK, bbox_to_anchor=(0.99, 1.0))
    fig.suptitle("Lysine fraction, (VPGIG)₉₆ → N-terminal VPGKG block, single chain, TEA γ = 6",
                 x=0.01, y=0.99, ha="left", fontsize=11.5, color=INK, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=Path("figures/g6-kfrac/g6-kfrac.csv"))
    parser.add_argument("--out", type=Path, default=Path("figures/g6-kfrac/g6-kfrac.png"))
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)

    if args.plot_only:
        rows = [{k: (v if k == "run" else float(v)) for k, v in r.items()}
                for r in csv.DictReader(open(args.csv))]
    else:
        names = [run_name(p, t) for p in FRACTIONS for t in TEMPS]
        with Pool(args.jobs) as pool:
            results = pool.map(analyse, names)
        rows = [r for r in results if r]
        missing = [n for n, r in zip(names, results) if not r]
        if missing:
            print(f"not finished ({len(missing)}): {', '.join(missing)}")
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
        for r in sorted(rows, key=lambda r: (r["lys_pct"], r["temperature_K"])):
            print(f"{r['run']:<22}{r['rg_nm']:>7.2f} ±{r['rg_sem_nm']:.2f}"
                  f"{r['collapse_fraction']:>6.0%}{r['span_p99_nm']:>7.1f}"
                  f"{r['chain_sasa_nm2']:>8.1f}{r['lys_sasa_per_k_nm2']:>7.3f}")
    plot(rows, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
