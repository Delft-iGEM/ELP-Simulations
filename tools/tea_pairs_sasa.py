"""SASA of the tea_pairs runs: one sequence, free single chain vs 16-chain brush.

    .venv/bin/python tools/tea_pairs_sasa.py               # analyse + plot
    .venv/bin/python tools/tea_pairs_sasa.py --plot-only   # re-plot from the CSV

Runs (tools/tea_pairs.py): tea-{rgd,lys}-{best,worst}-{free,brush}-{280,350}k,
CALVADOS-TEA gamma 3, 300 ns, the first 50 ns discarded. Free: one chain in a
60 nm box. Brush: 16 chains grafted at 0.04 chains/nm^2 through residue 0 (the
"Z" anchor bead) in a 20 x 20 x 50 nm box. The two geometries share every
residue except residue 0 (L or H free, Z grafted).

Coarse-grained Shrake-Rupley SASA (tools/sasa.py, bead radius sigma/2, probe
0.14 nm, minimum image in x, y and z). In the brush every bead of all 16 chains
can occlude, so neighbouring chains count as burial; the grafting wall is not a
bead and does not. Values are per site, averaged over every site of every
chain in the frame, then over frames:

  rgd_site      R, G and D beads of one RGD                 (all runs)
  domain_site   the 23-mer cell-binding domain               (rgd runs)
  lys_site      one lysine bead, the crosslink site          (lys runs)
  chain         a whole chain                                 (all runs, every 5th sample)

Sampled every 10th frame (0.7 ns). Error bars are the SEM from 20 block
averages of the per-frame means.

Writes figures/tea-pairs/tea-pairs-sasa.csv and, one figure per set and
quantity, tea-pairs-{rgd,lys}-<quantity>.png.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

STRIDE = 10
CHAIN_EVERY = 5                  # whole-chain SASA on every 5th sampled frame
DISCARD_NS = 50.0
DOMAIN = "LDASTVYAVTGRGDSPASSAASA"
SETS = {"rgd": ("best", "worst"), "lys": ("best", "worst")}
TEMPS = (280, 350)
GEOMS = ("free", "brush")
OUT_DIR = Path("figures/tea-pairs")
SOURCE = {"rgd-best": "rgdarr19-3n-1c", "rgd-worst": "rgdarr16-pairs-mid",
          "lys-best": "lysarr12-even-4", "lys-worst": "lysarr02-blk4-mid"}


def run_name(s: str, rank: str, geom: str, t: int) -> str:
    return f"tea-{s}-{rank}-{geom}-{t}k"


def analyse(name: str) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads
    from tools.tea_analysis import block_average

    try:
        top, dcd = trajectory_paths(name)
    except FileNotFoundError:
        return None
    if not Path(dcd).is_file():
        return None
    traj = md.load(str(dcd), top=str(top))
    runtime = Path("simulations") / name / "runtime"
    meta = read_metadata(runtime)
    ps = float(meta["time_per_frame_ps"])
    traj = traj[np.arange(traj.n_frames) * ps / 1000.0 >= DISCARD_NS][::STRIDE]
    box = traj.unitcell_lengths[0]
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    n_chain = traj.n_atoms // len(seq)
    if n_chain * len(seq) != traj.n_atoms:
        raise ValueError(f"{name}: {traj.n_atoms} beads is not a whole number of chains")
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))

    groups = {"rgd_site": [list(range(m.start(), m.start() + 3))
                           for m in re.finditer("RGD", seq)]}
    if name.startswith("tea-rgd"):
        groups["domain_site"] = [list(range(m.start(), m.start() + len(DOMAIN)))
                                 for m in re.finditer(DOMAIN, seq)]
    else:
        groups["lys_site"] = [[i] for i, c in enumerate(seq) if c == "K"]
    # every site of every chain, as global bead indices
    sites = {g: [[c * len(seq) + i for i in s] for c in range(n_chain) for s in ss]
             for g, ss in groups.items()}
    want = sorted({i for ss in sites.values() for s in ss for i in s})

    per = {g: [] for g in sites}
    chain = []
    for k, frame in enumerate(traj.xyz):
        got = dict(zip(want, sasa_of_beads(frame, radii, which=want, box=box)))
        for g, ss in sites.items():
            per[g].append(np.mean([sum(got[i] for i in s) for s in ss]))
        if k % CHAIN_EVERY == 0:
            chain.append(sasa_of_beads(frame, radii, box=box).sum() / n_chain)

    _, rank, geom, t = name[4:].split("-")
    row = {"run": name, "set": name[4:7], "rank": rank, "geometry": geom,
           "T": int(t[:-1]), "n_chains": n_chain, "n_residues": len(seq),
           "frames": len(traj), "chain_frames": len(chain),
           "sites": {g: len(groups[g]) for g in groups}}
    for g, v in list(per.items()) + [("chain", chain)]:
        row[f"{g}_nm2"], row[f"{g}_sem_nm2"] = block_average(np.array(v))
    return row


# ---- figures ---------------------------------------------------------------------

INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
# geometry is the identity here, so not the temperature blue/red of the g6 figures
GEOM_STYLE = {"free": ("#2f8a68", "single free chain"),
              "brush": ("#8a5cc2", "grafted, 16-chain brush")}
RANK_STYLE = {"best": ("-", "o"), "worst": ((0, (4, 2)), "s")}
QUANTITIES = {
    "rgd": [("rgd-sasa", "rgd_site", "SASA per RGD (nm²)", "RGD solvent exposure"),
            ("domain-sasa", "domain_site", "SASA per 23-mer domain (nm²)",
             "Cell-binding domain (23-mer) exposure"),
            ("chain-sasa", "chain", "SASA per chain (nm²)", "Whole-chain SASA")],
    "lys": [("lys-sasa", "lys_site", "SASA per lysine (nm²)", "Lysine (crosslink site) exposure"),
            ("rgd-sasa", "rgd_site", "SASA per RGD (nm²)", "RGD solvent exposure"),
            ("chain-sasa", "chain", "SASA per chain (nm²)", "Whole-chain SASA")],
}
SET_NAME = {"rgd": "RGD-arrangement pair", "lys": "lysine-arrangement pair"}


def plot(rows: list[dict], out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    for s, quantities in QUANTITIES.items():
        for tag, key, ylabel, title in quantities:
            fig, ax = plt.subplots(figsize=(6.6, 4.8))
            fig.patch.set_facecolor(SURFACE)
            ax.set_facecolor(SURFACE)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color(GRID)
            ax.tick_params(colors=MUTED, labelsize=8.5)
            ax.grid(color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            for geom, (colour, _) in GEOM_STYLE.items():
                for rank, (ls, marker) in RANK_STYLE.items():
                    pts = sorted((r["T"], r[f"{key}_nm2"], r[f"{key}_sem_nm2"]) for r in rows
                                 if r["set"] == s and r["geometry"] == geom
                                 and r["rank"] == rank and not math.isnan(r[f"{key}_nm2"]))
                    if not pts:
                        continue
                    T, y, e = map(np.array, zip(*pts))
                    ax.plot(T, y, color=colour, lw=2, ls=ls, zorder=2)
                    ax.errorbar(T, y, yerr=e, fmt=marker, ms=7, color=colour, elinewidth=1.1,
                                capsize=0, mec=SURFACE, mew=1.3, zorder=3)
                    ax.annotate(f"{y[-1]:.2f}", (T[-1], y[-1]), xytext=(7, 0),
                                textcoords="offset points", va="center", fontsize=7.8,
                                color=INK)
            ax.set_xticks(TEMPS)
            ax.set_xticklabels([f"{t} K" for t in TEMPS])
            ax.set_xlim(TEMPS[0] - 12, TEMPS[1] + 18)
            ax.set_xlabel("Temperature", fontsize=9.5, color=INK)
            ax.set_ylabel(ylabel, fontsize=9.5, color=INK)
            ranks = {rk: SOURCE[f"{s}-{rk}"] for rk in RANK_STYLE}
            handles = ([plt.Line2D([], [], color=c, lw=2, label=lab)
                        for c, lab in GEOM_STYLE.values()]
                       + [plt.Line2D([], [], color=MUTED, lw=2, ls=ls, marker=m, ms=6.5,
                                     mec=SURFACE, label=f"{rk} ({ranks[rk]})")
                          for rk, (ls, m) in RANK_STYLE.items()])
            ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
                      ncol=2, frameon=False, fontsize=8.3, labelcolor=INK, handlelength=3)
            ax.set_title(f"{title}, {SET_NAME[s]}", loc="left", fontsize=10.5, color=INK,
                         fontweight="bold")
            r0 = next(r for r in rows if r["set"] == s)
            per = ("per site, averaged over every site of every chain" if key != "chain"
                   else "per chain, averaged over the chains")
            fig.text(0.01, 0.01,
                     f"Coarse-grained SASA (bead radius σ/2, probe 0.14 nm), {per};\n"
                     f"50–300 ns, every {STRIDE}th frame; bar = 20-block SEM. Brush: neighbouring "
                     f"chains occlude, the grafting wall does not.\nTEA γ 3, "
                     f"{int(r0['n_residues'])} aa, one run per point.", fontsize=7.6, color=MUTED, linespacing=1.5)
            fig.tight_layout(rect=(0, 0.1, 1, 1))
            out = out_dir / f"tea-pairs-{s}-{tag}.png"
            fig.savefig(out, dpi=200, facecolor=SURFACE)
            plt.close(fig)
            print(f"wrote {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args(argv)
    out_csv = OUT_DIR / "tea-pairs-sasa.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    text = ("run", "set", "rank", "geometry", "sites")

    if args.plot_only:
        rows = [{k: (v if k in text else int(v) if k == "T" else float(v))
                 for k, v in r.items()} for r in csv.DictReader(open(out_csv))]
        plot(rows, OUT_DIR)
        return 0

    names = [run_name(s, rk, g, t) for s, ranks in SETS.items() for rk in ranks
             for g in GEOMS for t in TEMPS]
    with Pool(args.jobs) as pool:
        results = pool.map(analyse, names)
    missing = [n for n, r in zip(names, results) if not r]
    if missing:
        print(f"missing ({len(missing)}): {', '.join(missing)}")
    rows = [r for r in results if r]
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in text, k))
    for r in rows:
        r["sites"] = ";".join(f"{g}={n}" for g, n in r["sites"].items())
        for k in keys:
            r.setdefault(k, float("nan"))
    order = ["run", "set", "rank", "geometry", "T", "n_chains", "n_residues", "frames",
             "chain_frames", "sites"]
    fields = order + [k for k in keys if k not in order]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}")
    print(f"{'run':<28}{'RGD':>13}{'domain':>14}{'Lys':>13}{'chain':>15}")
    for r in rows:
        cell = lambda k: (f"{r[k + '_nm2']:>7.2f}±{r[k + '_sem_nm2']:<5.2f}"
                          if not math.isnan(r[k + "_nm2"]) else f"{'-':>13}")
        print(f"{r['run']:<28}{cell('rgd_site')}{cell('domain_site')} {cell('lys_site')}"
              f"{r['chain_nm2']:>9.1f}±{r['chain_sem_nm2']:<5.1f}")
    plot(rows, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
