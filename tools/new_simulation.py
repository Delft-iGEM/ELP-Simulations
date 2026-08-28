"""
sim new — scaffold a new surface-attached ELP simulation.

Prompts for an ELP sequence, its surface concentration, how many copies to
simulate, and how many steps to run, then creates
simulations/<name>/prepare.py (in the same style as
simulations/surfaceattached-multiple) and immediately runs it to produce a
ready runtime/ folder.

Usage:
    uv run python -m tools.new_simulation
    uv run sim new

Or non-interactively:
    uv run sim new --sequence VPGIGVPGIGVPGIG... --concentration 0.02 \\
        --nmol 500 --steps 7070000 --name my-elp-surface
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Annotated, NamedTuple

import typer

app = typer.Typer(add_completion=False, no_args_is_help=False)

# Fixed physical constants shared with the existing surface-attached templates.
Z_HEIGHT = 50.0  # nm, box height
Z_WALL = 1.9  # nm, depth of the surface-attachment potential well

# OpenMM requires the nonbonded cutoff to be less than half the box size.
# CALVADOS's default electrostatics cutoff (cutoff_yu) is 4.0 nm, so the box
# footprint must be kept above this floor even if the requested
# concentration/nmol would otherwise call for a tighter box.
MIN_BOX_L = 10.0  # nm




def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    return Path.cwd()


def _valid_amino_acids() -> set[str]:
    csv_path = _project_root() / "residues_CALVADOS2.csv"
    with open(csv_path, newline="") as f:
        letters = {row["one"] for row in csv.DictReader(f)}
    # "Z" is applied automatically to mark the surface-attachment point;
    # users shouldn't type it themselves.
    return letters - {"Z"}


def _slugify(name: str) -> str:
    slug = "".join(c if c.isalnum() else "-" for c in name.strip().lower())
    slug = "-".join(filter(None, slug.split("-")))
    if not slug:
        raise typer.BadParameter("Name must contain at least one letter or digit.")
    return slug


def _validate_sequence(value: str) -> str:
    seq = "".join(value.split()).upper()
    if len(seq) < 2:
        raise typer.BadParameter("Sequence must be at least 2 residues long.")
    valid = _valid_amino_acids()
    invalid = sorted(set(seq) - valid)
    if invalid:
        raise typer.BadParameter(
            f"Sequence contains unknown amino-acid letters: {', '.join(invalid)}"
        )
    return seq


class BoxPlan(NamedTuple):
    l_box_x: float
    l_box_y: float
    nx: int
    ny: int
    spacing: float
    clamped: bool
    achieved_concentration: float
    min_sane_spacing: float
    spacing_too_tight: bool
    aspect_ratio: float


def plan_lattice(nmol: int) -> tuple[int, int]:
    """Most-square rectangular grid (nx, ny) with nx * ny == nmol exactly.

    An exactly-filled rectangle is what lets the periodic images tile the
    surface seamlessly: every chain then sees the same neighbours at the same
    spacing, including the ones sitting at the box edge, so no chain is less
    crowded than any other. The old square spiral did not tile — for nmol=12 it
    fills a 4x3 block of cells while the box was only sqrt(12)=3.46 cells wide,
    so each chain at the edge sat ~0.5 spacings too close to its own periodic
    image.
    """
    divisor = max(k for k in range(1, math.isqrt(nmol) + 1) if nmol % k == 0)
    return nmol // divisor, divisor


def suggest_nmol(nmol: int, max_aspect: float = 2.0, search: int = 24) -> int | None:
    """Nearest molecule count whose most-square lattice is no worse than `max_aspect`.

    Prime-ish counts (13, 37, ...) can only be laid out as a long strip, which
    makes for a very elongated box; this points at a nearby count that doesn't.
    """
    for delta in range(1, search + 1):
        for candidate in (nmol - delta, nmol + delta):
            if candidate < 1:
                continue
            nx, ny = plan_lattice(candidate)
            if max(nx, ny) / min(nx, ny) <= max_aspect:
                return candidate
    return None


def smallest_untilted_nmol(concentration: float, nmol: int, search: int = 400) -> int | None:
    """Smallest count >= nmol whose lattice clears MIN_BOX_L at the full concentration.

    Below that the short side of the box would fall under the cutoff floor and
    plan_box has to spread the chains out, quietly lowering the density; going
    up to this many molecules keeps the density exactly as requested instead.
    """
    spacing = 1.0 / math.sqrt(concentration)
    for candidate in range(nmol, nmol + search + 1):
        nx, ny = plan_lattice(candidate)
        if min(nx, ny) * spacing >= MIN_BOX_L and max(nx, ny) / min(nx, ny) <= 2.0:
            return candidate
    return None


def plan_box(n_residues: int, concentration: float, nmol: int) -> BoxPlan:
    """Periodic box + grafting lattice for `nmol` chains at `concentration`.

    Shared by `sim new` and the analyze.ipynb planning cell so both agree on
    the same box math.

    The chains are laid out on an exactly-filled nx-by-ny lattice and the box
    is exactly nx*spacing by ny*spacing, so the periodic images continue the
    lattice without a seam: the simulation represents an infinite grafted
    surface, and a chain at the edge is as crowded as one in the middle. That
    is also why there is no empty margin around the grid — a margin would make
    this a finite patch whose rim chains are under-crowded.

    Molecules do still cross the periodic boundary and get wrapped in the
    trajectory, which looks like them teleporting across the box; that is a
    representation artefact of correct periodic dynamics, and the analyze
    notebook unwraps it rather than the box trying to avoid it.
    """
    nx, ny = plan_lattice(nmol)
    aspect_ratio = max(nx, ny) / min(nx, ny)

    # The lattice spacing *is* the surface concentration: one chain per
    # spacing^2 of surface, so spacing = 1/sqrt(concentration).
    spacing_raw = 1.0 / math.sqrt(concentration)

    # OpenMM needs the nonbonded cutoff below half the box, so the *short* side
    # of the box has to clear MIN_BOX_L. The only way to widen the box without
    # breaking the tiling (nmol is fixed) is to spread the chains further
    # apart, which lowers the concentration — reported loudly by the caller.
    scale = max(1.0, MIN_BOX_L / (min(nx, ny) * spacing_raw))
    clamped = scale > 1.0
    spacing = round(spacing_raw * scale, 3)

    # Exact multiples of the spacing: this is what makes the tiling seamless.
    l_box_x = round(nx * spacing, 6)
    l_box_y = round(ny * spacing, 6)
    achieved_concentration = nmol / (l_box_x * l_box_y)

    # Rough lower bound on a chain's own footprint (random-walk end-to-end
    # estimate at the ~0.38 nm CALVADOS bond length). If molecules start
    # packed tighter than this, they'll begin heavily overlapped and the
    # simulation is very likely to diverge to NaN once dynamics starts.
    min_sane_spacing = 0.38 * math.sqrt(n_residues)
    spacing_too_tight = spacing < min_sane_spacing

    return BoxPlan(
        l_box_x,
        l_box_y,
        nx,
        ny,
        spacing,
        clamped,
        achieved_concentration,
        min_sane_spacing,
        spacing_too_tight,
        aspect_ratio,
    )


def plan_steps(steps: int, max_n_save: int = 7000) -> tuple[int, int, int]:
    """Save frequency/frame count that yields ~1000 saved frames, capped at max_n_save.

    Returns (n_save, n_frames, actual_steps) — actual_steps is steps rounded
    to an exact multiple of n_save.
    """
    n_save = min(max_n_save, max(1, steps // 1000))
    n_frames = max(1, round(steps / n_save))
    actual_steps = n_save * n_frames
    return n_save, n_frames, actual_steps


def lattice_positions(nmol: int, spacing: float) -> list[tuple[float, float]]:
    """Grafting-point (x, y) of every molecule, in nm, on the tiling lattice.

    Row-major over the nx-by-ny lattice from plan_lattice, each molecule at the
    *centre* of its spacing-by-spacing cell. The generated prepare.py's
    build_sim() computes exactly the same positions, so previews (e.g. the
    analyze.ipynb planning cell) can reproduce the real starting layout.
    """
    nx, _ny = plan_lattice(nmol)
    return [((i % nx + 0.5) * spacing, (i // nx + 0.5) * spacing) for i in range(nmol)]


PREPARE_TEMPLATE = '''import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md

def build_sim(sim: Sim):
     components = sim.components

     # Chains sit on an exactly-filled __NX__ x __NY__ lattice, one per
     # spacing x spacing cell, and the box is exactly that lattice wide. The
     # periodic images therefore continue the lattice without a seam, so this
     # is an infinite grafted surface and every chain is equally crowded.
     # Placing each chain at the *centre* of its cell (the +0.5) keeps it as far
     # from the box edge as the lattice allows.
     nx, ny = __NX__, __NY__
     spacing = __SPACING__

     ibead = 0
     i = 0
     for comp in components:
          # CALVADOS builds each chain's starting conformation with a lateral
          # offset and extent of its own, so xinit's x/y centre is not (0, 0).
          # Recentre it on the lattice point, otherwise the whole grafting
          # pattern sits offset from the lattice the box is built around.
          xy_centre = 0.5 * (comp.xinit[:, :2].min(axis=0) + comp.xinit[:, :2].max(axis=0))
          xinit_centred = comp.xinit - np.array([xy_centre[0], xy_centre[1], 0.0])

          for idx in range(comp.nmol):
               j = ibead + comp.nbeads

               x0 = (i % nx + 0.5) * spacing
               y0 = (i // nx + 0.5) * spacing

               # Bead 0 is the "Z"-tagged bead. CALVADOS gives "Z" a molecular
               # weight of -2 which the +2 N-terminus patch cancels to exactly
               # 0, and OpenMM holds zero-mass particles completely fixed — so
               # this lattice point is where the chain stays grafted for the
               # whole run, in x, y and z.
               pos = xinit_centred + np.array([x0, y0, 2.0])
               sim.pos[ibead:j] = pos

               ibead = j
               i += 1

     md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)


# Job settings for Delft Blue
partition = "__PARTITION__"
runtime = "__WALLTIME__"
cpu_per_task = "__CPU_PER_TASK__"

sim_name = Path(__file__).parent.name

box = [__L_X__, __L_Y__, __Z_HEIGHT__]
N_save = __N_SAVE__
N_frames = __N_FRAMES__

# OpenMM runs in nm. The wall keeps the "Z"-tagged end of each ELP anchored
# near the surface at z ~= 0.
z_wall = __Z_WALL__

sequences: dict[str, str] = {
     "__SEQ_NAME__": "__SEQUENCE__"
}
sequences["__SEQ_NAME__"] = f"Z{sequences['__SEQ_NAME__'][1:]}"


if __name__ == "__main__":
     path = Path(__file__).parent.resolve()
     cwd = Path(os.getcwd())

     runtime_dir = path / "runtime"
     runtime_dir.mkdir(exist_ok=True)

     fasta_file = runtime_dir / "molecules.fasta"

     # Config preparation
     config = Config(
          sysname = sim_name,
          box = box,
          temp = 293.15,
          ionic = 0.19,
          pH = 7.5,
          ext_force = True,
          ext_force_expr = f'step({z_wall}-z)*0.5*({z_wall}-z)^2',
          topol = 'grid',
          wfreq = N_save,
          steps = N_frames*N_save,
          runtime = 0,
          platform = '__PLATFORM__',
          restart = 'checkpoint',
          frestart = 'restart.chk',
          verbose = True
     )

     components = Components(
          # Defaults
          molecule_type = 'protein',
          nmol = __NMOL__, # number of molecules
          restraint = False,
          charge_termini = 'both',
          fresidues = str(cwd / "residues_CALVADOS2.csv"),
          ffasta = str(fasta_file)
     )

     for k in sequences.keys():
          components.add(name=k)

     # Prepare runtime - You probably don't need to edit below this

     # Write config and components
     config.write(str(runtime_dir), "config.yaml")
     components.write(str(runtime_dir), "components.yaml")

     # Create FASTA
     fasta_content = "\\n".join([f">{k}\\n{v}" for k, v in sequences.items()])

     fasta_file.write_text(fasta_content)

     # Move job.sh and run.py
     run_file = (cwd / "template/run.py").read_text()
     (runtime_dir / "run.py").write_text(run_file)

     job_file = (cwd / "template/job.sh").read_text()

     job_file = job_file \\
          .replace("{sim_name}", sim_name) \\
          .replace("{partition}", partition) \\
          .replace("{runtime}", runtime) \\
          .replace("{cpu_per_task}", cpu_per_task)

     (runtime_dir / "job.sh").write_text(job_file)
'''


def create_simulation(
    sequence: str,
    concentration: float,
    nmol: int,
    steps: int,
    name: str,
    platform: str = "CUDA",
    partition: str = "gpu-a100",
    walltime: str = "24:30:00",
    cpu_per_task: str = "18",
    max_n_save: int = 7000,
) -> Path:
    """Scaffold simulations/<name>/prepare.py and run it. Returns the sim folder."""

    if concentration <= 0:
        raise typer.BadParameter("Concentration must be greater than 0.")
    if nmol < 1:
        raise typer.BadParameter("Number of ELPs must be at least 1.")
    if steps < 1:
        raise typer.BadParameter("Number of steps must be at least 1.")
    if platform not in {"CUDA", "CPU"}:
        raise typer.BadParameter("Platform must be 'CUDA' or 'CPU'.")

    seq = _validate_sequence(sequence)
    slug = _slugify(name)

    root = _project_root()
    sim_dir = root / "simulations" / slug
    if sim_dir.exists():
        raise typer.BadParameter(
            f"simulations/{slug}/ already exists — choose a different name."
        )

    box_plan = plan_box(len(seq), concentration, nmol)
    spacing = box_plan.spacing
    min_sane_spacing = box_plan.min_sane_spacing

    n_save, n_frames, actual_steps = plan_steps(steps, max_n_save=max_n_save)

    content = (
        PREPARE_TEMPLATE
        .replace("__SPACING__", str(spacing))
        .replace("__NX__", str(box_plan.nx))
        .replace("__NY__", str(box_plan.ny))
        .replace("__L_X__", str(box_plan.l_box_x))
        .replace("__L_Y__", str(box_plan.l_box_y))
        .replace("__Z_HEIGHT__", str(Z_HEIGHT))
        .replace("__N_SAVE__", str(n_save))
        .replace("__N_FRAMES__", str(n_frames))
        .replace("__Z_WALL__", str(Z_WALL))
        .replace("__SEQ_NAME__", slug)
        .replace("__SEQUENCE__", seq)
        .replace("__NMOL__", str(nmol))
        .replace("__PLATFORM__", platform)
        .replace("__PARTITION__", partition)
        .replace("__WALLTIME__", walltime)
        .replace("__CPU_PER_TASK__", cpu_per_task)
    )

    sim_dir.mkdir(parents=True)
    (sim_dir / "prepare.py").write_text(content)

    typer.echo(f"✓  Created simulations/{slug}/prepare.py")
    typer.echo(f"   box:        [{box_plan.l_box_x}, {box_plan.l_box_y}, {Z_HEIGHT}] nm")
    typer.echo(f"   lattice:    {box_plan.nx} x {box_plan.ny} grafting points, {spacing} nm apart, "
               f"tiling the box exactly (infinite surface via the periodic images)")
    typer.echo(f"   density:    {round(box_plan.achieved_concentration, 5)} chains/nm^2")
    if box_plan.clamped:
        typer.echo(
            f"⚠  the requested {concentration} chains/nm^2 would make the short side of the box "
            f"{round(min(box_plan.nx, box_plan.ny) / math.sqrt(concentration), 2)} nm, below the "
            f"{MIN_BOX_L} nm floor required by CALVADOS's default cutoffs — the chains were spread "
            f"out to {spacing} nm instead ({round(box_plan.achieved_concentration, 5)} chains/nm^2)."
        )
        bigger = smallest_untilted_nmol(concentration, nmol)
        if bigger is not None:
            bx, by = plan_lattice(bigger)
            typer.echo(
                f"   fix:        use nmol={bigger} ({bx} x {by}) to get the full "
                f"{concentration} chains/nm^2 with no clamping."
            )
    if box_plan.aspect_ratio > 2.0:
        better = suggest_nmol(nmol)
        hint = f" — {better} molecules would tile {'x'.join(map(str, plan_lattice(better)))}" if better else ""
        typer.echo(
            f"⚠  {nmol} molecules only tile as {box_plan.nx} x {box_plan.ny}, giving a very "
            f"elongated box{hint}."
        )
    if box_plan.spacing_too_tight:
        typer.echo(
            f"⚠  spacing ({spacing} nm) is tight for a {len(seq)}-residue chain "
            f"(rough own-size estimate: {min_sane_spacing:.2f} nm) — molecules will start "
            f"heavily overlapped and the simulation is likely to diverge to NaN. "
            f"Consider a lower concentration and/or fewer molecules."
        )
    typer.echo(f"   molecules:  {nmol}")
    typer.echo(f"   platform:   {platform}")
    if actual_steps != steps:
        typer.echo(f"   steps:      {actual_steps} (rounded from {steps} to a multiple of the {n_save}-step save frequency)")
    else:
        typer.echo(f"   steps:      {actual_steps}")

    typer.echo(f"▶  Preparing: {slug}")
    result = subprocess.run(
        [sys.executable, "-m", f"simulations.{slug}.prepare"], cwd=str(root)
    )
    if result.returncode != 0:
        raise typer.Exit(result.returncode)

    typer.echo(f"✓  Done — runtime/ folder ready at simulations/{slug}/runtime/")
    typer.echo(f"   Next: sim run {slug}   or   sim submit {slug}")

    return sim_dir


@app.command()
def new(
    sequence: Annotated[
        str,
        typer.Option(prompt="ELP sequence (single-letter amino acids)"),
    ],
    concentration: Annotated[
        float,
        typer.Option(prompt="Surface concentration (ELP chains / nm^2)"),
    ],
    nmol: Annotated[
        int,
        typer.Option(prompt="Number of ELP molecules"),
    ],
    steps: Annotated[
        int,
        typer.Option(prompt="Number of simulation steps"),
    ],
    name: Annotated[
        str,
        typer.Option(prompt="Simulation folder name (under simulations/)"),
    ],
) -> None:
    """Scaffold and prepare a new surface-attached ELP simulation."""
    create_simulation(sequence, concentration, nmol, steps, name)


if __name__ == "__main__":
    app()
