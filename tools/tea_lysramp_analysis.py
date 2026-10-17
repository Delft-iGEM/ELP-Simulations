"""Rg and SASA for the lysine-ramp series (tools/tea_lysramp.py), and its figure.

    .venv/bin/python tools/tea_lysramp_analysis.py                  # analyse + plot
    .venv/bin/python tools/tea_lysramp_analysis.py --plot-only      # re-plot the CSV

Per run, over 50-400 ns, everything tools/tea_g6_analysis.py measures (Rg
against the theta-state Rg, collapse fraction, extent, K-K contacts, RGD and
23-mer SASA) plus the whole-chain SASA. RGD and 23-mer SASA are divided by the
number of copies (4 in every sequence here), so they read per binding site.
Unfinished runs are skipped; re-run once they are in. Writes figures/g6-kspc-kend/g6-lysramp.csv
and figures/g6-kspc-kend/g6-lysramp-{rg,rg-both,rgd-sasa,chain-sasa}.png, with the lysine
count on the x axis as % of the 49 replaceable pentads.
"""

from __future__ import annotations

import argparse
import csv
import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_lysramp import COUNTS, SERIES, TEMPERATURES

SASA_STRIDE = 10
DISCARD_NS = 50.0
N_SITES = 4
N_REPLACEABLE = 49           # pentads that can be VPGKG; the x axis is n / 49 in %


def run_name(series: str, n: int, t: float) -> str:
    return f"g6-{series}{n:02d}-{int(t)}k"


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
        row[k.replace("_nm2", "_per_site_nm2")] = row.pop(k) / N_SITES

    top, dcd = trajectory_paths(name)
    traj = md.load(str(dcd), top=str(top))
    meta = read_metadata(Path("simulations") / name / "runtime")
    t = np.arange(traj.n_frames) * float(meta["time_per_frame_ps"]) / 1000.0
    traj = traj[t >= DISCARD_NS]
    box = traj.unitcell_lengths[0]
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    chain = [sasa_of_beads(f, radii, box=box).sum() for f in traj.xyz[::SASA_STRIDE]]
    row["chain_sasa_nm2"], row["chain_sasa_sem_nm2"] = block_average(chain)
    row["series"] = name[3:7]
    return row


# ---- figure ------------------------------------------------------------------

COLD = "#4269d0"
HOT = "#c1442e"
INK = "#1a1a1a"
MUTED = "#6b6f76"
GRID = "#dcdee2"
SURFACE = "#fcfcfb"
# Temperature is the colour (as in g6-kfrac.png); placement is line and marker.
STYLE = {"kspc": ("spaced", "-", "o"), "kend": ("from both ends", (0, (4, 2)), "s")}


PANELS = [("rg", "rg_nm", "rg_sem_nm", "Radius of gyration (nm)"),
          ("rg-both", "rg_nm", "rg_sem_nm", "Radius of gyration (nm)"),
          ("rgd-sasa", "rgd_sasa_per_site_nm2", "rgd_sasa_sem_per_site_nm2",
           "SASA per RGD (nm²)"),
          ("chain-sasa", "chain_sasa_nm2", "chain_sasa_sem_nm2", "SASA, whole chain (nm²)")]


def plot(rows: list[dict], out_dir: Path) -> None:
    """One figure per quantity: <out_dir>/g6-lysramp-{rg,rgd-sasa,chain-sasa}.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = [(int(T), c) for T, c in zip(TEMPERATURES, (COLD, HOT))]
    out_dir.mkdir(parents=True, exist_ok=True)
    for tag, key, sem, label in PANELS:
        # the Rg figures are slide versions: transparent, untitled; rg has the spaced
        # series only, rg-both the spaced and the from-both-ends series
        slide = tag.startswith("rg")
        styles = {"kspc": STYLE["kspc"]} if tag == "rg" else STYLE
        ring = dict(mec=SURFACE, mew=0 if slide else 1.3)   # no surface halo when transparent
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        fig.patch.set_facecolor(SURFACE)
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for T, colour in temps:
            for series, (_, ls, marker) in styles.items():
                pts = sorted((x["n_lys"], x[key], x[sem]) for x in rows
                             if x["series"] == series
                             and int(round(x["temperature_K"])) == T
                             and not math.isnan(x[key]))
                if not pts:
                    continue
                n, y, e = map(np.array, zip(*pts))
                n = 100 * n / N_REPLACEABLE
                ax.plot(n, y, color=colour, linewidth=2, linestyle=ls, zorder=2)
                ax.errorbar(n, y, yerr=e, fmt=marker, ms=6.5, color=colour, ecolor=colour,
                            elinewidth=1, capsize=0, zorder=3, **ring)
        if key == "rg_nm":
            theta = rows[0]["rg_theta_nm"]
            ax.axhline(theta, color=MUTED, linewidth=1, linestyle=(0, (1, 2)), zorder=1)
            ax.text(0.99, theta, f"θ-state Rg {theta:.2f} nm",
                    transform=ax.get_yaxis_transform(),
                    va="bottom", ha="right", fontsize=8, color=MUTED)
        ax.set_ylabel(label, fontsize=9.5, color=INK)
        if slide:
            ax.set_xlabel("Percentage of replaced pentades (VPGIG → VPGKG)", fontsize=9.5,
                          color=INK)
        else:
            ax.set_xlabel(f"Lysine pentads, % of the {N_REPLACEABLE} replaceable pentads "
                          f"(VPGKG; the rest VPGIG)", fontsize=9, color=MUTED)
        ax.set_xticks([100 * c / N_REPLACEABLE for c in COUNTS])
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
        if slide:
            handles = [plt.Line2D([], [], color=c, lw=2, marker="o", ms=6.5, label=f"{T} K",
                                  **ring) for T, c in temps]
            if tag == "rg-both":
                handles = ([plt.Line2D([], [], color=c, lw=2, label=f"{T} K") for T, c in temps]
                           + [plt.Line2D([], [], color=INK, lw=2, linestyle=ls, marker=m,
                                         ms=6.5, mew=0, label=f"lysines {lab}")
                              for lab, ls, m in STYLE.values()])
        else:
            handles = ([plt.Line2D([], [], color=c, lw=2, label=f"{T} K") for T, c in temps]
                       + [plt.Line2D([], [], color=MUTED, lw=2, linestyle=ls, marker=m, ms=6.5,
                                     mec=SURFACE, label=lab) for lab, ls, m in STYLE.values()])
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13),
                  ncol=2 if tag == "rg-both" else 4,
                  frameon=False, fontsize=8.5, labelcolor=INK, handlelength=3)
        if not slide:
            ax.set_title(f"{label.split(' (')[0]}, lysine ramp on the 357-aa RGD construct",
                         loc="left", fontsize=10.5, color=INK, fontweight="bold")
        fig.tight_layout()
        out = out_dir / f"g6-lysramp-{tag}.png"
        if slide:
            fig.savefig(out, dpi=200, transparent=True)
        else:
            fig.savefig(out, dpi=200, facecolor=SURFACE)
        plt.close(fig)
        print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=Path("figures/g6-kspc-kend/g6-lysramp.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("figures/g6-kspc-kend"))
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)

    if args.plot_only:
        rows = [{k: (v if k in ("run", "series") else float(v)) for k, v in r.items()}
                for r in csv.DictReader(open(args.csv))]
    else:
        names = [run_name(s, n, t) for s in SERIES for n in COUNTS for t in TEMPERATURES]
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
        print(f"{'run':<18}{'Rg':>7}{'':>6}{'theta':>7}{'coll':>6}{'span':>7}"
              f"{'RGD':>7}{'chain':>8}")
        for r in rows:
            print(f"{r['run']:<18}{r['rg_nm']:>7.2f} ±{r['rg_sem_nm']:.2f}"
                  f"{r['rg_theta_nm']:>7.2f}{r['collapse_fraction']:>6.0%}"
                  f"{r['span_p99_nm']:>7.1f}{r['rgd_sasa_per_site_nm2']:>7.2f}"
                  f"{r['chain_sasa_nm2']:>8.1f}")
    plot(rows, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
