"""Collapsed vs extended states in the gamma-6 single-chain set.

    .venv/bin/python tools/tea_g6_states.py --csv figures/g6/g6-states.csv --plot figures/g6/g6-states.png

Per run, over 50-400 ns, every frame:

  collapsed     fraction of frames with Rg below the theta-state Rg
                (0.21 N^0.5 nm, tools/tea_analysis.py), the time that is in ns,
                and how it is spent: number of collapsed episodes, their mean
                and longest length
  RGD exposure  SASA per RGD tripeptide averaged separately over collapsed and
                extended frames, in nm^2 and as a percentage of the same three
                beads on a straight, fully stretched chain (only their bonded
                neighbours occlude them there - the most exposed RGD can be)

The comparison that isolates conformation is collapsed vs extended *at one
temperature*; across temperatures the bead sizes and contacts change too.
Errors are block-averaged over 20 consecutive blocks, using only the blocks
that contain frames in that state; with fewer than 3 such blocks there is no
error and the mean rests on a handful of correlated visits.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from tools.tea_analysis import rg_theta
from tools.tea_g6_analysis import DISCARD_NS, FRAMES_EXPECTED, _whole

N_BLOCKS = 20
STRETCHED_BOND_NM = 0.38


def _state_mean(x: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Mean of x over frames in `mask`, and the block s.e. over blocks that have any."""
    if not mask.any():
        return float("nan"), float("nan")
    size = len(x) // N_BLOCKS
    blocks = [x[i * size:(i + 1) * size][mask[i * size:(i + 1) * size]]
              for i in range(N_BLOCKS)]
    means = np.array([b.mean() for b in blocks if len(b)])
    sem = means.std(ddof=1) / math.sqrt(len(means)) if len(means) >= 3 else float("nan")
    return float(x[mask].mean()), float(sem)


def _episodes(mask: np.ndarray) -> np.ndarray:
    """Lengths, in frames, of the runs of consecutive True."""
    edges = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    return np.flatnonzero(edges == -1) - np.flatnonzero(edges == 1)


def analyse_run(name: str) -> tuple[dict, np.ndarray] | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads

    try:
        top, dcd = trajectory_paths(name)
    except FileNotFoundError:
        return None
    if not Path(dcd).is_file():
        return None
    traj = md.load(str(dcd), top=str(top))
    if traj.n_frames < FRAMES_EXPECTED:
        return None

    runtime = Path("simulations") / name / "runtime"
    meta = read_metadata(runtime)
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    ps = float(meta["time_per_frame_ps"])
    t = np.arange(traj.n_frames) * ps / 1000.0
    traj = traj[t >= DISCARD_NS]
    box = traj.unitcell_lengths[0]

    xyz = np.stack([_whole(f, box) for f in traj.xyz])
    dev = xyz - xyz.mean(axis=1, keepdims=True)
    rg = np.sqrt((dev ** 2).sum(axis=2).mean(axis=1))
    n = len(seq)
    theta = rg_theta(n)
    collapsed = rg < theta
    ep = _episodes(collapsed) * ps / 1000.0

    row = {
        "run": name,
        "temperature_K": float(meta["temperature_K"]),
        "n_residues": n,
        "rg_theta_nm": theta,
        "ns_used": len(rg) * ps / 1000.0,
        "collapse_fraction": float(collapsed.mean()),
        "collapsed_ns": float(collapsed.sum() * ps / 1000.0),
        "n_collapsed_episodes": len(ep),
        "mean_episode_ns": float(ep.mean()) if len(ep) else 0.0,
        "longest_episode_ns": float(ep.max()) if len(ep) else 0.0,
    }

    starts = [m.start() for m in re.finditer("RGD", seq)]
    row["n_rgd"] = len(starts)
    per_frame = np.column_stack([rg / theta, np.full(len(rg), np.nan)])
    if starts:
        radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
        idx = [i for s in starts for i in range(s, s + 3)]
        sasa = np.array([sasa_of_beads(f, radii, which=idx, box=box).sum()
                         for f in traj.xyz]) / len(starts)
        line = np.zeros((n, 3))
        line[:, 0] = np.arange(n) * STRETCHED_BOND_NM
        ref = sasa_of_beads(line, radii, which=idx).sum() / len(starts)
        per_frame[:, 1] = sasa
        row["rgd_sasa_stretched_nm2"] = ref
        for state, mask in (("collapsed", collapsed), ("extended", ~collapsed)):
            m, e = _state_mean(sasa, mask)
            row[f"rgd_sasa_{state}_nm2"] = m
            row[f"rgd_sasa_{state}_sem_nm2"] = e
            row[f"rgd_exposure_{state}_pct"] = 100 * m / ref
        row["rgd_sasa_diff_nm2"] = row["rgd_sasa_extended_nm2"] - row["rgd_sasa_collapsed_nm2"]
    return row, per_frame


def _plot(per_run: dict[str, tuple[dict, np.ndarray]], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    bins = np.linspace(0.5, 2.5, 21)
    centres = (bins[:-1] + bins[1:]) / 2
    for temp, colour in ((285.0, "#2a6fdb"), (370.0, "#d94a2b")):
        pts = np.concatenate([pf for row, pf in per_run.values()
                              if row["temperature_K"] == temp and row["n_rgd"]])
        if not len(pts):
            continue
        ref = np.mean([row["rgd_sasa_stretched_nm2"] for row, _ in per_run.values()
                       if row["temperature_K"] == temp and row["n_rgd"]])
        which = np.digitize(pts[:, 0], bins) - 1
        mean = np.array([pts[which == b, 1].mean() if (which == b).sum() >= 50 else np.nan
                         for b in range(len(centres))])
        ax.plot(centres, 100 * mean / ref, "o-", color=colour, label=f"{temp:.0f} K")
    ax.axvline(1.0, color="0.4", ls="--", lw=1)
    ax.text(0.98, ax.get_ylim()[0], "theta ", ha="right", va="bottom", color="0.4")
    ax.set_xlabel("Rg / Rg(theta)")
    ax.set_ylabel("RGD SASA, % of stretched chain")
    ax.set_title("RGD exposure against compaction, all g6 runs with RGD")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="*")
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args(argv)

    names = args.runs or sorted(p.name for p in Path("simulations").glob("g6-*")
                                if (p / "runtime").is_dir())
    with ProcessPoolExecutor(args.jobs) as pool:
        results = dict(zip(names, pool.map(analyse_run, names)))
    per_run = {k: v for k, v in results.items() if v}
    skipped = [k for k, v in results.items() if not v]
    if not per_run:
        print("no finished g6-* runs")
        return 1
    rows = [row for row, _ in per_run.values()]

    print(f"{'run':<32}{'T':>5}{'coll':>7}{'ns':>7}{'eps':>6}{'mean ns':>9}{'max ns':>8}"
          f"{'RGD coll %':>14}{'RGD ext %':>13}")
    for r in rows:
        line = (f"{r['run']:<32}{r['temperature_K']:>5.0f}{r['collapse_fraction']:>7.0%}"
                f"{r['collapsed_ns']:>7.0f}{r['n_collapsed_episodes']:>6}"
                f"{r['mean_episode_ns']:>9.2f}{r['longest_episode_ns']:>8.1f}")
        if r["n_rgd"]:
            line += (f"{r['rgd_exposure_collapsed_pct']:>14.1f}"
                     f"{r['rgd_exposure_extended_pct']:>13.1f}")
        print(line)
    if skipped:
        print(f"\nnot finished ({len(skipped)}): {', '.join(skipped)}")
    if args.csv:
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, restval="")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
    if args.plot:
        _plot(per_run, args.plot)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
