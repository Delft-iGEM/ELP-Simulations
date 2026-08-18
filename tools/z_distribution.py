"""
Plot the distribution of a specified amino acid's position along the z axis,
pooled over every frame (timepoint) of a CALVADOS/OpenMM trajectory.

Usage:
    uv run python -m tools.z_distribution --pdb simulations/<sim>/runtime/top.pdb \\
        --dcd simulations/<sim>/runtime/<sim>.dcd --residue K

    uv run sim distribution --pdb ... --dcd ... --residue LYS
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Annotated, Optional

import matplotlib

matplotlib.use("Agg")  # renders without a display (SSH / headless machines)
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import typer

app = typer.Typer(add_completion=False)


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    return Path.cwd()


def _one_to_three_letter() -> dict[str, str]:
    csv_path = _project_root() / "residues_CALVADOS2.csv"
    with open(csv_path, newline="") as f:
        return {row["one"]: row["three"] for row in csv.DictReader(f)}


def _resolve_resname(residue: str) -> str:
    residue = residue.strip().upper()
    if len(residue) == 1:
        mapping = _one_to_three_letter()
        if residue not in mapping:
            raise typer.BadParameter(f"Unknown one-letter amino-acid code: {residue}")
        return mapping[residue]
    return residue


def collect_z_positions(pdb: Path, dcd: Path, residue: str) -> tuple[np.ndarray, int]:
    """Return (z-coordinates pooled over every frame, number of frames)."""
    resname = _resolve_resname(residue)

    u = mda.Universe(str(pdb), str(dcd))
    atoms = u.select_atoms(f"resname {resname} and name CA")
    if len(atoms) == 0:
        raise typer.BadParameter(
            f"No '{resname}' residues found in {pdb} (checked topology, not trajectory)."
        )

    n_frames = len(u.trajectory)
    z = np.empty((n_frames, len(atoms)), dtype=float)
    for i, _ in enumerate(u.trajectory):
        z[i] = atoms.positions[:, 2]

    return z.ravel(), n_frames


def plot_z_distribution(z: np.ndarray, residue: str, n_frames: int, bins: int = 60) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(z, bins=bins, density=True, color="#4c72b0", alpha=0.75, edgecolor="white")

    try:
        from scipy.stats import gaussian_kde

        kde = gaussian_kde(z)
        xs = np.linspace(z.min(), z.max(), 400)
        ax.plot(xs, kde(xs), color="#c44e52", linewidth=2)
    except ImportError:
        pass

    ax.set_xlabel("z (Å)")
    ax.set_ylabel("density")
    ax.set_title(f"z-distribution of {residue.upper()}\n{len(z)} samples across {n_frames} frames")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


@app.command()
def distribution(
    pdb: Annotated[
        Path,
        typer.Option(exists=True, help="Path to the topology PDB file (e.g. runtime/top.pdb)."),
    ],
    dcd: Annotated[
        Path,
        typer.Option(exists=True, help="Path to the trajectory DCD file."),
    ],
    residue: Annotated[
        str,
        typer.Option(help="Amino acid to plot: one-letter (K) or three-letter (LYS) code."),
    ],
    output: Annotated[
        Optional[Path],
        typer.Option("--output", "-o", help="Where to save the plot. Defaults next to the DCD file."),
    ] = None,
    bins: Annotated[int, typer.Option(help="Number of histogram bins.")] = 60,
    show: Annotated[bool, typer.Option(help="Also display the plot interactively.")] = False,
) -> None:
    """Plot the z-axis distribution of a residue, pooled across all trajectory frames."""
    z, n_frames = collect_z_positions(pdb, dcd, residue)
    fig = plot_z_distribution(z, residue, n_frames, bins=bins)

    out_path = output or dcd.with_name(f"{dcd.stem}_{residue.upper()}_z_distribution.png")
    fig.savefig(out_path, dpi=150)
    typer.echo(f"✓  Saved z-distribution plot to {out_path}")

    if show:
        plt.show()


if __name__ == "__main__":
    app()
