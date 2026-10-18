"""SASA and Rg of the grafted RGD arrangements (rgdg01..24).

    .venv/bin/python tools/rgdg_analysis.py               # analyse finished runs + plot
    .venv/bin/python tools/rgdg_analysis.py --plot-only   # re-plot from the CSV

The rgdg runs (tools/rgdarr_graft.py) put the 24 sequences of the free-chain
RGD sweep on a surface: 16 chains hanging from surface-bonded lysines (0.3 of
them bonded at t=0, fixed), 0.03 chains/nm^2, 500 ns, 293.15 K, no TEA, no
crosslinking. Runs whose metadata.csv does not say `completed` are skipped, so
re-running as more finish just adds them.

Per run, from 100 ns on:

  rgd     SASA of the 4 RGD tripeptides of a chain, summed, mean over the 16
          chains                                     every 10th frame
  23mer   the same for the 4 LDASTVYAVTGRGDSPASSAASA domains
                                                     every 10th frame
  chain   SASA of a whole chain, mean over the 16    every 80th frame
  rg      radius of gyration of a chain (unwrapped along its bonds), mean
          over the 16                                every frame

Every bead of every chain can occlude (minimum image); the surface itself is a
wall, not beads, and does not. Error = 20-block SEM over time.

Writes figures/rgdg/rgdg-analysis.csv and, per quantity, rgdg-<q>.png in the
style of figures/free-rgdarr/sasa-24-arrangements-*-schematic.png: sorted rows,
colour/shape = RGD domains at chain ends, a chain sketch per row, and the same
three runs highlighted as in the free-chain figures (01-orig, and the free-chain
RGD + linker winner 19-3n-1c and loser 16-pairs-mid).
"""

from __future__ import annotations

import argparse
import csv
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

START_NS = 100.0
SITE_STRIDE = 10
CHAIN_STRIDE = 80
DOMAIN = "LDASTVYAVTGRGDSPASSAASA"
OUT_DIR = Path("figures/rgdg")
OUT_CSV = OUT_DIR / "rgdg-analysis.csv"
FREE_CSV = Path("figures/free-rgdarr/sasa-24-arrangements.csv")
# highlighted as in the free-chain figures; winner/loser are by free-chain 23-mer SASA
ROLE = {"01-orig": "original", "19-3n-1c": "free winner", "16-pairs-mid": "free loser"}
FIELDS = ["run", "label", "frames_analysed", "window_ns",
          "rgd_nm2", "rgd_sem_nm2", "motif_23mer_nm2", "motif_23mer_sem_nm2",
          "chain_sasa_nm2", "chain_sasa_sem_nm2", "rg_nm", "rg_sem_nm"]
# tag, column, sem column, x label
PANELS = [
    ("rgd", "rgd_nm2", "rgd_sem_nm2", "RGD SASA per chain (nm²)"),
    ("23mer", "motif_23mer_nm2", "motif_23mer_sem_nm2", "RGD + linker SASA per chain (nm²)"),
    ("chain-sasa", "chain_sasa_nm2", "chain_sasa_sem_nm2", "Whole-chain SASA (nm²)"),
    ("rg", "rg_nm", "rg_sem_nm", "Radius of gyration (nm)"),
]


def finished_runs() -> list[str]:
    from tools.metadata import read_metadata
    out = []
    for d in sorted(Path("simulations").glob("rgdg[0-9][0-9]-*")):
        try:
            if read_metadata(d / "runtime").get("status") == "completed":
                out.append(d.name)
        except (FileNotFoundError, OSError):
            pass
    return out


def analyse(name: str) -> dict:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, find_motif, sasa_of_beads
    from tools.tea_analysis import block_average
    from tools.tea_g6_analysis import _whole

    runtime = Path("simulations") / name / "runtime"
    top, dcd = trajectory_paths(name)
    traj = md.load(str(dcd), top=str(top))
    ns_per_frame = float(read_metadata(runtime)["time_per_frame_ps"]) / 1000.0
    first = max(0, int(np.ceil(START_NS / ns_per_frame)) - 1)
    n_total = traj.n_frames
    traj = traj[first:]
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    n_res = len(seq)
    n_chain = traj.n_atoms // n_res
    box = traj.unitcell_lengths[0]
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))

    rgd_res = [i for s, e in find_motif(seq, "RGD") for i in range(s, e)]
    dom_res = [i for s, e in find_motif(seq, DOMAIN) for i in range(s, e)]
    assert len(rgd_res) == 12 and len(dom_res) == 92, name
    rgd = np.array([c * n_res + i for c in range(n_chain) for i in rgd_res])
    dom = np.array([c * n_res + i for c in range(n_chain) for i in dom_res])
    in_rgd = np.isin(dom, rgd)

    rgd_s, dom_s = [], []
    for f in range(0, traj.n_frames, SITE_STRIDE):
        s = sasa_of_beads(traj.xyz[f], radii, which=dom, box=box)
        dom_s.append(s.sum() / n_chain)
        rgd_s.append(s[in_rgd].sum() / n_chain)       # the RGD beads are inside the domain
    chain = [sasa_of_beads(traj.xyz[f], radii, box=box).sum() / n_chain
             for f in range(0, traj.n_frames, CHAIN_STRIDE)]
    rg = np.empty(traj.n_frames)
    for f in range(traj.n_frames):
        per = []
        for c in range(n_chain):
            w = _whole(traj.xyz[f, c * n_res:(c + 1) * n_res], box)
            per.append(np.sqrt(((w - w.mean(axis=0)) ** 2).sum(axis=1).mean()))
        rg[f] = np.mean(per)

    row = {"run": name, "label": name[len("rgdg"):], "frames_analysed": traj.n_frames,
           "window_ns": f"{first * ns_per_frame:.0f}-{n_total * ns_per_frame:.0f}"}
    row["rgd_nm2"], row["rgd_sem_nm2"] = block_average(np.array(rgd_s))
    row["motif_23mer_nm2"], row["motif_23mer_sem_nm2"] = block_average(np.array(dom_s))
    row["chain_sasa_nm2"], row["chain_sasa_sem_nm2"] = block_average(np.array(chain))
    row["rg_nm"], row["rg_sem_nm"] = block_average(rg)
    print(f"  {name}: rgd {row['rgd_nm2']:.2f}  23mer {row['motif_23mer_nm2']:.2f}  "
          f"chain {row['chain_sasa_nm2']:.1f}  rg {row['rg_nm']:.2f}", flush=True)
    return row


def plot(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from tools.plot_sasa_24 import draw

    free = {r["run"]: r for r in csv.DictReader(open(FREE_CSV))}
    for tag, key, sem, xlabel in PANELS:
        data = [{"run": r["label"], "copy_starts": free[r["label"]]["copy_starts"],
                 "copies_at_chain_end": float(free[r["label"]]["copies_at_chain_end"]),
                 key: float(r[key]), sem: float(r[sem])} for r in rows]
        role = {k: v for k, v in ROLE.items() if any(d["run"] == k for d in data)}
        draw(plt, data, (tag, key, sem, None, "", "", True), role, schematic=True,
             out=OUT_DIR / f"rgdg-{tag}.png", xlabel=xlabel)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args(argv)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not args.plot_only:
        names = finished_runs()
        print(f"{len(names)} finished rgdg runs")
        with Pool(min(args.jobs, len(names))) as pool:
            rows = pool.map(analyse, names)
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {OUT_CSV}")
    rows = list(csv.DictReader(open(OUT_CSV)))
    plot(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
