"""Analysis for the CALVADOS-TEA sweep: Rg, d_Rg, collapse fraction, nu.

    .venv/bin/python tools/tea_analysis.py            # every prepared sweep run
    .venv/bin/python tools/tea_analysis.py --csv out.csv

The primary metric is the authors' own LCST/UCST discriminator::

    d_Rg = Rg(350 K) - Rg(300 K)        < 0 = LCST,  > 0 = UCST

Errors are block-averaged, not the naive SD over correlated frames: with a
~ns correlation time and frames every 100 ps the naive interval is several
times too narrow, and the sequences here differ by one residue in five, so the
interval decides whether anything is being measured at all.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

# Theta-state Flory reference (the task's values): Rg = rho0 * N^nu.
RHO0_NM = 0.21
NU_THETA = 0.5
DISCARD_NS = 50.0


def rg_theta(n_residues: int) -> float:
    """Theta-state Rg in nm. 420 residues -> 4.30 nm."""
    return RHO0_NM * n_residues ** NU_THETA


def rg_series(sim_name: str, discard_ns: float = DISCARD_NS) -> tuple[np.ndarray, np.ndarray]:
    """(time in ns, Rg in nm) for one run, production window only.

    Rg is the geometric, unweighted radius of gyration on coordinates rebuilt by
    walking the backbone under the minimum-image convention — the same
    definition the rest of this project uses, so TEA numbers are comparable with
    the existing free-chain results.
    """
    import MDAnalysis as mda

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata

    top, dcd = trajectory_paths(sim_name)
    u = mda.Universe(str(top), str(dcd))
    if len(u.segments) != 1:
        raise ValueError(f"{sim_name}: {len(u.segments)} chains, expected a single free chain")
    box = u.dimensions[:3].copy()
    idx = u.segments[0].atoms.indices

    rg = np.empty(len(u.trajectory))
    for i, _ in enumerate(u.trajectory):
        pos = u.atoms.positions[idx]
        bond = np.diff(pos, axis=0)
        bond -= np.round(bond / box) * box
        whole = np.concatenate([pos[:1], pos[:1] + np.cumsum(bond, axis=0)], axis=0)
        dev = whole - whole.mean(axis=0)
        rg[i] = np.sqrt((dev ** 2).sum(axis=1).mean())

    ps = float(read_metadata(Path("simulations") / sim_name / "runtime")["time_per_frame_ps"])
    t = np.arange(len(rg)) * ps / 1000.0
    keep = t >= discard_ns
    return t[keep], rg[keep] / 10.0          # Angstrom -> nm


def block_average(x: np.ndarray, n_blocks: int = 20) -> tuple[float, float]:
    """(mean, standard error) from `n_blocks` consecutive blocks.

    The blocks are the samples: with enough of them each is long compared with
    the correlation time, so their scatter is the real uncertainty on the mean
    rather than the frame-to-frame noise.
    """
    x = np.asarray(x, dtype=float)
    if len(x) < n_blocks * 2:
        n_blocks = max(2, len(x) // 2)
    size = len(x) // n_blocks
    blocks = x[:size * n_blocks].reshape(n_blocks, size).mean(axis=1)
    return float(x.mean()), float(blocks.std(ddof=1) / math.sqrt(n_blocks))


def block_curve(x: np.ndarray) -> list[tuple[int, float]]:
    """(block size, s.e.) as the blocks grow — the plateau is the converged error."""
    x = np.asarray(x, dtype=float)
    out = []
    size = 1
    while size <= len(x) // 8:
        n = len(x) // size
        blocks = x[:size * n].reshape(n, size).mean(axis=1)
        out.append((size, float(blocks.std(ddof=1) / math.sqrt(n))))
        size *= 2
    return out


def collapse_fraction(rg_nm: np.ndarray, n_residues: int) -> float:
    """Fraction of frames below the theta-state reference."""
    return float(np.mean(np.asarray(rg_nm) < rg_theta(n_residues)))


def internal_scaling(sim_name: str, discard_ns: float = DISCARD_NS,
                     stride: int = 10, fit_range: tuple[int, int] = (5, 100)):
    """nu from <|r_i - r_j|> vs |i - j|, fitted over `fit_range` in residues.

    Returns (separations, mean distance in nm, nu, prefactor). The short-range
    end is excluded because bond geometry, not chain statistics, sets it.
    """
    import MDAnalysis as mda

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata

    top, dcd = trajectory_paths(sim_name)
    u = mda.Universe(str(top), str(dcd))
    box = u.dimensions[:3].copy()
    idx = u.segments[0].atoms.indices
    n = len(idx)
    ps = float(read_metadata(Path("simulations") / sim_name / "runtime")["time_per_frame_ps"])

    total = np.zeros(n)
    counts = np.zeros(n)
    used = 0
    for i, _ in enumerate(u.trajectory):
        if i * ps / 1000.0 < discard_ns or i % stride:
            continue
        pos = u.atoms.positions[idx]
        bond = np.diff(pos, axis=0)
        bond -= np.round(bond / box) * box
        whole = np.concatenate([pos[:1], pos[:1] + np.cumsum(bond, axis=0)], axis=0)
        for sep in range(1, n):
            d = np.linalg.norm(whole[sep:] - whole[:-sep], axis=1)
            total[sep] += d.sum()
            counts[sep] += len(d)
        used += 1
    if not used:
        raise ValueError(f"{sim_name}: no frames after discarding {discard_ns} ns")

    seps = np.arange(1, n)
    mean_d = (total[1:] / counts[1:]) / 10.0      # nm
    lo, hi = fit_range
    window = (seps >= lo) & (seps <= min(hi, n - 1))
    nu, log_rho = np.polyfit(np.log(seps[window]), np.log(mean_d[window]), 1)
    return seps, mean_d, float(nu), float(np.exp(log_rho))


def analyse(sim_names: list[str], discard_ns: float = DISCARD_NS) -> list[dict]:
    from tools.metadata import read_metadata

    rows = []
    for name in sim_names:
        meta = read_metadata(Path("simulations") / name / "runtime")
        t, rg = rg_series(name, discard_ns)
        mean, sem = block_average(rg)
        n_res = int(float(meta.get("n_residues") or 0))
        rows.append({
            "run": name,
            "temperature_K": float(meta.get("temperature_K") or "nan"),
            "tea": meta.get("tea", ""),
            "tea_gamma": meta.get("tea_gamma", ""),
            "n_residues": n_res,
            "frames_used": len(rg),
            "ns_used": float(t[-1] - t[0]) if len(t) else 0.0,
            "rg_nm": mean,
            "rg_sem_nm": sem,
            "rg_theta_nm": rg_theta(n_res) if n_res else float("nan"),
            "collapse_fraction": collapse_fraction(rg, n_res) if n_res else float("nan"),
        })
    return rows


def delta_rg(rows: list[dict], low: float = 300.0, high: float = 350.0) -> list[dict]:
    """d_Rg = Rg(high) - Rg(low) per sequence, with errors added in quadrature.

    Negative = LCST (collapses on heating), positive = UCST.
    """
    by_sequence: dict[str, dict[float, dict]] = {}
    for row in rows:
        stem = row["run"].rsplit("-", 1)[0]
        by_sequence.setdefault(stem, {})[row["temperature_K"]] = row
    out = []
    for stem, at in sorted(by_sequence.items()):
        if low not in at or high not in at:
            continue
        a, b = at[low], at[high]
        d = b["rg_nm"] - a["rg_nm"]
        err = math.hypot(a["rg_sem_nm"], b["rg_sem_nm"])
        out.append({"sequence": stem, "rg_low_nm": a["rg_nm"], "rg_high_nm": b["rg_nm"],
                    "d_rg_nm": d, "d_rg_sem_nm": err,
                    "verdict": ("LCST" if d < -2 * err else
                                "UCST" if d > 2 * err else "no resolved change")})
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="*", help="simulation names (default: every tea-* run)")
    parser.add_argument("--csv", type=Path, help="write the per-run table here")
    parser.add_argument("--discard-ns", type=float, default=DISCARD_NS)
    args = parser.parse_args(argv)

    names = args.runs or sorted(
        p.name for p in Path("simulations").glob("tea-*")
        if (p / "runtime").exists())
    if not names:
        print("no prepared tea-* runs found")
        return 1

    rows = analyse(names, args.discard_ns)
    print(f"{'run':<22}{'T (K)':>7}{'Rg (nm)':>18}{'theta':>8}{'collapsed':>11}{'frames':>9}")
    for r in rows:
        print(f"{r['run']:<22}{r['temperature_K']:>7.0f}"
              f"{r['rg_nm']:>12.3f} ±{r['rg_sem_nm']:.3f}{r['rg_theta_nm']:>8.2f}"
              f"{r['collapse_fraction']:>10.1%}{r['frames_used']:>9,}")

    print(f"\n{'sequence':<22}{'Rg(300)':>10}{'Rg(350)':>10}{'d_Rg (nm)':>18}  verdict")
    for d in delta_rg(rows):
        print(f"{d['sequence']:<22}{d['rg_low_nm']:>10.3f}{d['rg_high_nm']:>10.3f}"
              f"{d['d_rg_nm']:>12.3f} ±{d['d_rg_sem_nm']:.3f}  {d['verdict']}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
