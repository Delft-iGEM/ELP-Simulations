"""Compare the domain exposure of two designs, with a real error bar.

``analysis_batch``'s domain-exposure step pools all 16 chains of a run into one
mean per copy position. That is the right thing for ranking the four copy slots
against each other, but it throws away the only genuine replicate a single run
has: the 16 chains. Ask "is design A better than design B" of those pooled
numbers and there is nothing to put an interval on.

So this keeps the per-chain totals. Each chain contributes one number, the SASA
summed over its four copies of the motif, time-averaged over the production
window. Sixteen chains give sixteen samples, and two designs can then be
compared with a Welch t-test.

What the p-value does and does not establish
--------------------------------------------
Two designs are two *different runs*, with different seeds and different initial
placements. A significant difference therefore says "these two runs differ", not
"the design caused it" — seed and design are confounded and no amount of chains
inside one run separates them. Read it as a magnitude, with the caveat attached.

The causal claim comes from the paired test instead: inside a single run, compare
one copy slot against its own siblings. Same seed, same box, same chains, so the
only thing varying is position along the sequence. That is what established the
C-terminal slot's advantage (8 runs out of 8). Use this module for "how big",
that test for "is it real".

The chains within a run are also not fully independent of each other — they share
a box, a surface and a thermostat — so the sixteen samples are somewhat
correlated and the quoted standard error is, if anything, optimistic.

    .venv/bin/python tools/sasa_two_sample.py rgdarr08-all-cterm rgdarr23-staircase

Nothing else in the pipeline imports this; it is a reporting tool.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from tools.analysis_batch import trajectory_paths
from tools.sasa import PROBE_NM, bead_radii, chain_sequence, find_motif, sasa_of_beads

MOTIF = "LDASTVYAVTGRGDSPASSAASA"
RGD_OFFSET = MOTIF.find("RGD")

# Production window, matching the domain-exposure step: frames from 100 ns on,
# every 8th. The stride is there because SASA is the expensive part and the
# series is autocorrelated over far longer than 8 frames anyway.
FIRST_FRAME = 1428
STRIDE = 8


def per_chain_totals(sim_name: str, first_frame: int = FIRST_FRAME,
                     stride: int = STRIDE) -> dict[str, Any]:
    """Time-averaged SASA of all four motif copies, one number per chain.

    Returns the per-chain totals for the whole 23-mer and for the bare RGD
    tripeptides inside it, plus the free-domain reference each would reach with
    nothing around it.
    """
    import mdtraj as md

    top, dcd = trajectory_paths(sim_name)
    traj = md.load(str(dcd), top=str(top))
    if traj.unitcell_lengths is None:
        raise ValueError(f"{dcd} has no box information; occlusion by periodic "
                         "images cannot be evaluated.")

    sequence = chain_sequence(traj.topology, 0)
    hits = find_motif(sequence, MOTIF)
    if not hits:
        raise ValueError(f"{sim_name}: chain 0 does not contain {MOTIF}")
    n_chains = traj.topology.n_chains
    n_res = traj.topology.chain(0).n_residues
    radii = bead_radii(traj.topology)

    # Bead indices of every copy in every chain, in (chain, copy, position)
    # order so the flat SASA array can be reshaped straight back to that.
    which = np.array([chain * n_res + start + k
                      for chain in range(n_chains)
                      for start, _ in hits
                      for k in range(len(MOTIF))])

    frames = range(first_frame, traj.n_frames, stride)
    tot_motif = np.zeros((len(frames), n_chains))
    tot_rgd = np.zeros((len(frames), n_chains))
    for t, i in enumerate(frames):
        area = sasa_of_beads(
            np.asarray(traj.xyz[i], dtype=float), radii, which=which,
            box=np.asarray(traj.unitcell_lengths[i], dtype=float),
        ).reshape(n_chains, len(hits), len(MOTIF))
        tot_motif[t] = area.sum(axis=(1, 2))
        tot_rgd[t] = area[:, :, RGD_OFFSET:RGD_OFFSET + 3].sum(axis=(1, 2))

    first_copy = which[:len(MOTIF)]
    free = 4.0 * np.pi * (radii[first_copy] + PROBE_NM) ** 2
    return {
        "sim_name": sim_name,
        "n_copies": len(hits),
        "n_chains": n_chains,
        "n_frames_used": len(frames),
        "motif_per_chain": tot_motif.mean(axis=0),
        "rgd_per_chain": tot_rgd.mean(axis=0),
        "motif_free": float(len(hits) * free.sum()),
        "rgd_free": float(len(hits) * free[RGD_OFFSET:RGD_OFFSET + 3].sum()),
    }


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Welch's unequal-variance t-test: (t, degrees of freedom, two-sided p)."""
    from scipy import stats

    t, p = stats.ttest_ind(a, b, equal_var=False)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    return float(t), float(df), float(p)


def _block(label: str, a: np.ndarray, b: np.ndarray, free: float,
           name_a: str, name_b: str) -> str:
    t, df, p = welch(a, b)
    lines = [f"  {label}",
             f"    {name_a:<24} {a.mean():8.3f} +/- {a.std(ddof=1):5.3f} nm^2"
             f"   {100 * a.mean() / free:5.2f}% of free   range {a.min():.2f}-{a.max():.2f}",
             f"    {name_b:<24} {b.mean():8.3f} +/- {b.std(ddof=1):5.3f} nm^2"
             f"   {100 * b.mean() / free:5.2f}% of free   range {b.min():.2f}-{b.max():.2f}",
             f"    difference               {a.mean() - b.mean():+8.3f} nm^2"
             f"        {100 * (a.mean() - b.mean()) / free:+5.2f} points"
             f"   {100 * (a.mean() / b.mean() - 1):+5.2f}% relative",
             f"    Welch t = {t:.2f}, df = {df:.1f}, p = {p:.2g}"]
    return "\n".join(lines)


def report(name_a: str, name_b: str) -> str:
    """Human-readable comparison of two runs' domain exposure."""
    a = per_chain_totals(name_a)
    b = per_chain_totals(name_b)
    if a["n_copies"] != b["n_copies"]:
        raise ValueError(f"{name_a} has {a['n_copies']} copies of the motif but "
                         f"{name_b} has {b['n_copies']}; the totals are not comparable.")
    head = (f"{name_a} vs {name_b}: total SASA of all {a['n_copies']} motif copies\n"
            f"  one sample per chain, {a['n_chains']} and {b['n_chains']} chains, "
            f"time-averaged over {a['n_frames_used']} and {b['n_frames_used']} frames "
            f"from frame {FIRST_FRAME} stride {STRIDE}\n")
    tail = ("\n  Caveat: two runs means two seeds, so this tests whether these runs\n"
            "  differ, not whether the design is the cause. See the module docstring.")
    return "\n\n".join([
        head,
        _block(f"whole 23-mer {MOTIF}", a["motif_per_chain"], b["motif_per_chain"],
               a["motif_free"], name_a, name_b),
        _block("RGD tripeptides only", a["rgd_per_chain"], b["rgd_per_chain"],
               a["rgd_free"], name_a, name_b),
    ]) + tail


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        print(f"usage: {Path(sys.argv[0]).name} <run-a> <run-b>", file=sys.stderr)
        return 2
    print(report(argv[0], argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
