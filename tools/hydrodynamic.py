"""Hydrodynamic radius from a trajectory — the quantity light scattering measures.

Every other analysis in this project reports the radius of gyration, because Rg
is what a simulation naturally gives you: a mass distribution about a centre.
Dynamic light scattering does not measure that. It measures how fast the
molecule diffuses, and reports the radius of the sphere that would diffuse at
the same rate — the *hydrodynamic* radius Rh.

The two are not interchangeable. Rg/Rh depends on chain statistics: about 1.5
for a self-avoiding coil, 0.77 for a hard sphere, and anywhere between for a
real polymer. Converting one to the other needs that ratio assumed in advance,
which would beg the question in a validation against measured Rh. So this
computes Rh directly.

Kirkwood-Riseman
----------------
For a chain of N beads, the pre-averaged hydrodynamic radius is

    1 / Rh = (1 / N^2) * sum_{i != j} < 1 / r_ij >

averaged over frames. The double sum runs over all ordered pairs, so each
unordered pair appears twice; the sum is done here over the upper triangle and
doubled, which is the same thing and half the work. This is the estimator the
CALVADOS papers use, so a comparison against their numbers, or against an
experiment they were fitted to, is like for like.

Known limitations, none of which this code can fix:

* It is *pre-averaged*: the ensemble average of 1/r is taken rather than the
  hydrodynamic interaction being solved per configuration. That overestimates
  Rh slightly for expanded chains.
* It has no bead radius in it, so it is the radius of the chain's shape, not of
  the chain plus its hydration shell. For a coarse-grained chain of ~100+ beads
  the shape term dominates.
* It says nothing about solvent quality that the force field does not already
  say.

Periodic images
---------------
Distances are computed on the *whole* chain, rebuilt bond by bond under the
minimum image convention, not on the wrapped coordinates a DCD stores. A chain
that straddles the box boundary is written in two pieces, and the raw distance
between those pieces is a box length rather than a bond length — which would
make 1/r_ij enormous for exactly the pairs that matter most. Rebuilding first
is not optional.

This module is only used by the free-solution validation runs. It touches
nothing the surface-attached pipeline depends on.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

SIMULATIONS = Path("simulations")


def unwrap_chain(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Rebuild one chain's coordinates so consecutive beads are never a box apart.

    ``positions`` is (n_beads, 3) as stored in the trajectory, ``box`` the three
    side lengths in the same units. Walks the backbone taking the minimum image
    of every bond vector, so the result is a continuous chain wherever the
    trajectory happened to cut it.
    """
    bond = np.diff(positions, axis=0)
    bond -= np.round(bond / box) * box
    out = np.empty_like(positions)
    out[0] = positions[0]
    out[1:] = positions[0] + np.cumsum(bond, axis=0)
    return out


def rh_of_frame(positions: np.ndarray) -> float:
    """Kirkwood-Riseman Rh of one configuration, in the units of `positions`.

    Expects coordinates already made whole (see `unwrap_chain`).
    """
    n = len(positions)
    if n < 2:
        return float("nan")
    diff = positions[:, None, :] - positions[None, :, :]
    r = np.sqrt((diff ** 2).sum(-1))
    iu = np.triu_indices(n, k=1)
    inv = (1.0 / r[iu]).sum() * 2.0          # both orderings of every pair
    return float(n * n / inv)


def rg_of_frame(positions: np.ndarray) -> float:
    """Geometric radius of gyration, for the Rg/Rh ratio only."""
    dev = positions - positions.mean(axis=0)
    return float(np.sqrt((dev ** 2).sum(-1).mean()))


def analyse(sim_name: str, start_fraction: float = 0.2,
            stride: int = 1) -> dict[str, Any]:
    """Rh and Rg over a finished free-solution run.

    ``start_fraction`` of the trajectory is discarded as equilibration (0.2 =
    the first 20%, which is what the validation protocol asks for). ``stride``
    subsamples frames; the pair sum is O(N^2) per frame, so a 603-bead chain
    over 2857 frames is about a billion distances at stride 1.

    Returns the means, the standard error on *independent* samples estimated
    from the autocorrelation of the per-frame series, and the numbers needed to
    say what was measured.
    """
    import mdtraj as md

    runtime = SIMULATIONS / sim_name / "runtime"
    dcds = [p for p in runtime.glob("*.dcd") if not p.name.startswith("backup_")]
    if not dcds:
        raise FileNotFoundError(f"no trajectory in {runtime} — has {sim_name} run?")
    traj = md.load(str(dcds[0]), top=str(runtime / "top.pdb"))
    if traj.unitcell_lengths is None:
        raise ValueError(f"{dcds[0]} has no box information; cannot unwrap the chain.")

    n_frames = traj.n_frames
    first = int(round(start_fraction * n_frames))
    frames = range(first, n_frames, stride)

    rh, rg = [], []
    for i in frames:
        whole = unwrap_chain(np.asarray(traj.xyz[i], dtype=np.float64),
                             np.asarray(traj.unitcell_lengths[i], dtype=np.float64))
        rh.append(rh_of_frame(whole))
        rg.append(rg_of_frame(whole))
    rh = np.array(rh)
    rg = np.array(rg)

    return {
        "sim_name": sim_name,
        "n_beads": int(traj.n_atoms),
        "n_frames_total": n_frames,
        "n_frames_used": len(rh),
        "first_frame": first,
        "stride": stride,
        "rh_nm": float(rh.mean()),
        "rh_sem_nm": _sem(rh),
        "rh_std_nm": float(rh.std(ddof=1)),
        "rg_nm": float(rg.mean()),
        "rg_sem_nm": _sem(rg),
        "rg_over_rh": float(rg.mean() / rh.mean()),
        "n_eff": _n_eff(rh),
    }


def _statistical_inefficiency(series: np.ndarray) -> float:
    """g = 1 + 2 sum_t (1 - t/N) C(t), truncated at the first non-positive C(t).

    The same estimator the notebook's equilibration cell uses, so "independent
    samples" means the same thing in both places.
    """
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 3:
        return 1.0
    dx = x - x.mean()
    var = dx.dot(dx) / n
    if var <= 0:
        return 1.0
    g = 1.0
    for lag in range(1, n - 1):
        c = dx[:-lag].dot(dx[lag:]) / ((n - lag) * var)
        if c <= 0 and lag > 1:
            break
        g += 2.0 * c * (1.0 - lag / n)
    return max(1.0, g)


def _n_eff(series: np.ndarray) -> float:
    return float(len(series) / _statistical_inefficiency(series))


def _sem(series: np.ndarray) -> float:
    """Standard error on the mean of a correlated series."""
    n_eff = _n_eff(series)
    return float(series.std(ddof=1) / math.sqrt(max(n_eff, 1.0)))
