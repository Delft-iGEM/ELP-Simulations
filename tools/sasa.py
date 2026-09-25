"""Solvent-accessible surface area of a sequence motif, for coarse-grained beads.

"Is the cell-binding domain actually on the outside?" is a question about
exposure, and exposure is SASA. The complication is that every off-the-shelf
SASA routine assumes atoms: mdtraj's ``shrake_rupley`` looks the radius up by
chemical element, and every bead in a CALVADOS chain is element carbon. Run it
here and every residue would get carbon's 0.17 nm radius, when a CALVADOS bead
stands for a whole residue and is 0.45 to 0.68 nm across. The numbers would come
out looking fine and mean nothing.

So this implements Shrake-Rupley directly with the force field's own radii:
each bead gets ``sigma / 2`` from ``residues_CALVADOS2.csv``, the same sigma the
Ashbaugh-Hatch term uses, so the surface being measured is the surface the
simulation actually had.

What the number is and is not
-----------------------------
* It is a *coarse-grained* SASA. One bead per residue cannot resolve a side
  chain's own surface, so this measures whether the domain is buried inside the
  chain or the brush, not which of its atoms face the solvent.
* The probe radius is a convention. 0.14 nm (water) is the default and is what
  atomistic work uses, but with 0.3 nm beads it is a small correction. Absolute
  values depend on it; comparisons between runs at the same probe do not.
* There is no solvent in CALVADOS. "Accessible" means geometrically unoccluded
  by other beads, which is the right question for whether a ligand or a cell
  receptor could reach the motif.

Periodic images are respected: a bead can be occluded by a neighbouring chain
that the trajectory happens to have wrapped to the far side of the box.

Used by the notebook's RGD-domain cell. Nothing else in the pipeline depends on
it.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

# Water, the usual convention. See the module docstring on what it does and does
# not mean for beads this size.
PROBE_NM = 0.14

# Points per sphere in the Shrake-Rupley quadrature. 960 puts the discretisation
# error well under a percent, which is far smaller than the coarse-graining
# itself; raise it only if you are chasing a difference of that size.
N_SPHERE_POINTS = 960


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    return Path.cwd()


def sigmas_by_resname(csv_path: str | Path | None = None) -> dict[str, float]:
    """three-letter residue name -> CALVADOS sigma in nm.

    The table keys on the one-letter code and two rows share a three-letter
    name: "Z" is the tagged valine used as the surface anchor and "X" a tagged
    lysine. Both carry the same sigma as the residue they tag, so collapsing on
    the three-letter name is lossless — asserted here rather than assumed.
    """
    path = Path(csv_path) if csv_path else _project_root() / "residues_CALVADOS2.csv"
    out: dict[str, float] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            three, sigma = row["three"], float(row["sigmas"])
            if three in out and abs(out[three] - sigma) > 1e-9:
                raise ValueError(f"{three} has two different sigmas in {path}")
            out[three] = sigma
    return out


def bead_radii(topology, csv_path: str | Path | None = None) -> np.ndarray:
    """Radius of every bead in an mdtraj topology, nm — sigma/2 from the table."""
    sigma = sigmas_by_resname(csv_path)
    missing = sorted({a.residue.name for a in topology.atoms} - sigma.keys())
    if missing:
        raise ValueError(f"no CALVADOS sigma for residue name(s): {', '.join(missing)}")
    return np.array([sigma[a.residue.name] / 2.0 for a in topology.atoms])


def _sphere_points(n: int) -> np.ndarray:
    """`n` roughly equidistant points on the unit sphere (golden spiral)."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


def sasa_of_beads(positions: np.ndarray, radii: np.ndarray,
                  which: Sequence[int] | None = None,
                  box: np.ndarray | None = None,
                  probe: float = PROBE_NM,
                  n_points: int = N_SPHERE_POINTS) -> np.ndarray:
    """Shrake-Rupley SASA (nm^2) of the beads in `which`, or of all of them.

    `positions` (n_beads, 3) and `radii` (n_beads,) in nm. `box` gives the three
    side lengths for the minimum-image convention; pass None for a
    non-periodic system.

    Only the requested beads are sampled, but *every* bead can occlude them, so
    asking for a 23-residue domain out of a 16-chain brush costs 23 spheres, not
    thousands.
    """
    positions = np.asarray(positions, dtype=float)
    radii = np.asarray(radii, dtype=float)
    idx = np.arange(len(positions)) if which is None else np.asarray(which, dtype=int)
    grown = radii + probe
    unit = _sphere_points(n_points)

    out = np.empty(len(idx), dtype=float)
    for k, i in enumerate(idx):
        # only beads whose grown sphere can reach bead i's can occlude it
        d = positions - positions[i]
        if box is not None:
            d -= box * np.round(d / box)
        dist = np.linalg.norm(d, axis=1)
        reach = grown + grown[i]
        near = np.flatnonzero((dist < reach) & (np.arange(len(positions)) != i))
        if near.size == 0:
            out[k] = 4.0 * np.pi * grown[i] ** 2
            continue
        test = positions[i] + unit * grown[i]          # points on bead i's sphere
        # distance from each test point to each nearby bead, minimum-imaged
        rel = test[:, None, :] - (positions[i] + d[near])[None, :, :]
        buried = (np.einsum("pnc,pnc->pn", rel, rel) < grown[near] ** 2).any(axis=1)
        out[k] = 4.0 * np.pi * grown[i] ** 2 * (1.0 - buried.mean())
    return out


def find_motif(sequence: str, motif: str) -> list[tuple[int, int]]:
    """Every occurrence of `motif` in `sequence`, as (start, end) 0-based half-open."""
    hits, start = [], sequence.find(motif)
    while start != -1:
        hits.append((start, start + len(motif)))
        start = sequence.find(motif, start + 1)
    return hits


def chain_sequence(topology, chain_index: int = 0,
                   csv_path: str | Path | None = None) -> str:
    """One-letter sequence of a chain, read back from the topology's residue names."""
    path = Path(csv_path) if csv_path else _project_root() / "residues_CALVADOS2.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    # prefer the untagged one-letter code when two share a three-letter name
    three_to_one: dict[str, str] = {}
    for row in rows:
        three, one = row["three"], row["one"]
        if three not in three_to_one or one in "ACDEFGHIKLMNPQRSTVWY":
            three_to_one[three] = one
    chain = topology.chain(chain_index)
    return "".join(three_to_one[r.name] for r in chain.residues)
