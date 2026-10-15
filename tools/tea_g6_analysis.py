"""Analysis for the gamma-6, 285 / 370 K single-chain set (tools/tea_g6.py).

    .venv/bin/python tools/tea_g6_analysis.py              # every finished g6-* run
    .venv/bin/python tools/tea_g6_analysis.py --csv figures/g6/g6-285-370.csv

Per run, over 50-400 ns:

  rg          block-averaged Rg, the same definition as tools/tea_analysis.py
  collapsed   fraction of frames below the theta-state Rg
  extent      p99 of the chain's largest span, against half the box - a chain
              past half its box can meet its own periodic image
  kk          per-pair, per-frame K-K contact probability at 12 A (the lysine
              sweep's cutoff), all pairs and non-local pairs (>= 25 residues
              apart), the latter being the ones an adjacent-pentad block cannot
              satisfy trivially
  sasa        SASA summed over all copies of the RGD tripeptide and, where the
              sequence has it, of the 23-mer cell-binding domain, every 10th
              frame - the convention of sasa-24-arrangements.csv

Per sequence: d_Rg = Rg(370) - Rg(285), errors in quadrature, LCST if it is
more than 2 s.e. below zero.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from pathlib import Path

import numpy as np

from tools.tea_analysis import block_average, collapse_fraction, delta_rg, rg_theta

LOW, HIGH = 285.0, 370.0
DISCARD_NS = 50.0
FRAMES_EXPECTED = 5714
KK_CUTOFF_NM = 1.2
NONLOCAL = 25
SASA_STRIDE = 10
DOMAIN = "LDASTVYAVTGRGDSPASSAASA"


def _whole(pos: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Unwrap one chain by walking its bonds under minimum image."""
    bond = np.diff(pos, axis=0)
    bond -= np.round(bond / box) * box
    return np.concatenate([pos[:1], pos[:1] + np.cumsum(bond, axis=0)], axis=0)


def analyse_run(name: str) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads

    try:
        top, dcd = trajectory_paths(name)
    except FileNotFoundError:
        return None                          # not run yet: no topology or trajectory
    if not Path(dcd).is_file():
        return None
    traj = md.load(str(dcd), top=str(top))
    if traj.n_frames < FRAMES_EXPECTED:
        return None                          # still running, or cut short

    runtime = Path("simulations") / name / "runtime"
    meta = read_metadata(runtime)
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    ps = float(meta["time_per_frame_ps"])
    t = np.arange(traj.n_frames) * ps / 1000.0
    traj = traj[t >= DISCARD_NS]
    box = traj.unitcell_lengths[0]           # nm, cubic and constant

    xyz = np.stack([_whole(f, box) for f in traj.xyz])
    dev = xyz - xyz.mean(axis=1, keepdims=True)
    rg = np.sqrt((dev ** 2).sum(axis=2).mean(axis=1))
    span = (xyz.max(axis=1) - xyz.min(axis=1)).max(axis=1)
    rg_mean, rg_sem = block_average(rg)

    lys = [i for i, c in enumerate(seq) if c == "K"]
    pairs = list(itertools.combinations(lys, 2))
    if pairs:
        a, b = np.array(pairs).T
        d = np.linalg.norm(xyz[:, a] - xyz[:, b], axis=2)
        contact = d < KK_CUTOFF_NM
        far = (b - a) >= NONLOCAL
        p_all = float(contact.mean())
        p_far = float(contact[:, far].mean()) if far.any() else float("nan")
    else:
        p_all = p_far = float("nan")

    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    rgd = [m.start() for m in re.finditer("RGD", seq)]
    dom = [m.start() for m in re.finditer(DOMAIN, seq)]
    rgd_idx = [i for s in rgd for i in range(s, s + 3)]
    dom_idx = [i for s in dom for i in range(s, s + len(DOMAIN))]
    want = sorted(set(rgd_idx) | set(dom_idx))
    sasa_rgd, sasa_dom = [], []
    for frame in traj.xyz[::SASA_STRIDE]:
        per = dict(zip(want, sasa_of_beads(frame, radii, which=want, box=box)))
        sasa_rgd.append(sum(per[i] for i in rgd_idx))
        if dom:
            sasa_dom.append(sum(per[i] for i in dom_idx))

    n = len(seq)
    return {
        "run": name,
        "temperature_K": float(meta["temperature_K"]),
        "tea_gamma": meta.get("tea_gamma", ""),
        "n_residues": n,
        "n_lys": len(lys),
        "frames_used": traj.n_frames,
        "rg_nm": rg_mean,
        "rg_sem_nm": rg_sem,
        "rg_theta_nm": rg_theta(n),
        "collapse_fraction": collapse_fraction(rg, n),
        "span_p99_nm": float(np.percentile(span, 99)),
        "half_box_nm": float(box[0] / 2),
        "kk_p_pair": p_all,
        "kk_p_pair_nonlocal": p_far,
        "rgd_sasa_nm2": block_average(sasa_rgd)[0],
        "rgd_sasa_sem_nm2": block_average(sasa_rgd)[1],
        "domain_sasa_nm2": block_average(sasa_dom)[0] if sasa_dom else float("nan"),
        "domain_sasa_sem_nm2": block_average(sasa_dom)[1] if sasa_dom else float("nan"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="*")
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args(argv)

    names = args.runs or sorted(p.name for p in Path("simulations").glob("g6-*")
                                if (p / "runtime").is_dir())
    rows, skipped = [], []
    for name in names:
        row = analyse_run(name)
        (rows if row else skipped).append(row or name)
    if not rows:
        print("no finished g6-* runs")
        return 1

    print(f"{'run':<32}{'T':>5}{'Rg (nm)':>15}{'coll':>7}{'span99':>8}{'K-K p':>8}"
          f"{'K-K far':>9}{'RGD':>7}{'23-mer':>8}")
    for r in rows:
        print(f"{r['run']:<32}{r['temperature_K']:>5.0f}{r['rg_nm']:>9.2f} ±{r['rg_sem_nm']:.2f}"
              f"{r['collapse_fraction']:>7.0%}{r['span_p99_nm']:>8.1f}"
              f"{r['kk_p_pair']:>8.3f}{r['kk_p_pair_nonlocal']:>9.3f}"
              f"{r['rgd_sasa_nm2']:>7.2f}{r['domain_sasa_nm2']:>8.1f}")
    print(f"\n{'sequence':<28}{'Rg 285':>8}{'Rg 370':>8}{'d_Rg (nm)':>16}  verdict")
    for d in delta_rg(rows, LOW, HIGH):
        print(f"{d['sequence']:<28}{d['rg_low_nm']:>8.2f}{d['rg_high_nm']:>8.2f}"
              f"{d['d_rg_nm']:>10.2f} ±{d['d_rg_sem_nm']:.2f}  {d['verdict']}")
    if skipped:
        print(f"\nnot finished ({len(skipped)}): {', '.join(skipped)}")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
