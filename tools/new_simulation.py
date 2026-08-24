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
    l_box: float
    clamped: bool
    spacing: float
    min_sane_spacing: float
    spacing_too_tight: bool


def plan_box(n_residues: int, concentration: float, nmol: int) -> BoxPlan:
    """Surface footprint + starting-grid spacing for `nmol` chains at `concentration`.

    Shared by `sim new` and the analyze.ipynb planning cell so both agree on
    the same box math.
    """
    # Surface footprint sized so nmol molecules sit at the requested
    # concentration (chains / nm^2): L^2 = nmol / concentration.
    l_box_raw = math.sqrt(nmol / concentration)
    clamped = l_box_raw < MIN_BOX_L
    l_box = round(max(l_box_raw, MIN_BOX_L), 2)

    # Average spacing between molecules, recomputed from the final (possibly
    # clamped) box so the initial grid always fits inside it.
    spacing = round(l_box / math.sqrt(nmol), 3)

    # Rough lower bound on a chain's own footprint (random-walk end-to-end
    # estimate at the ~0.38 nm CALVADOS bond length). If molecules start
    # packed tighter than this, they'll begin heavily overlapped and the
    # simulation is very likely to diverge to NaN once dynamics starts.
    min_sane_spacing = 0.38 * math.sqrt(n_residues)
    spacing_too_tight = spacing < min_sane_spacing

    return BoxPlan(l_box, clamped, spacing, min_sane_spacing, spacing_too_tight)


def plan_steps(steps: int, max_n_save: int = 7000) -> tuple[int, int, int]:
    """Save frequency/frame count that yields ~1000 saved frames, capped at max_n_save.

    Returns (n_save, n_frames, actual_steps) — actual_steps is steps rounded
    to an exact multiple of n_save.
    """
    n_save = min(max_n_save, max(1, steps // 1000))
    n_frames = max(1, round(steps / n_save))
    actual_steps = n_save * n_frames
    return n_save, n_frames, actual_steps


def get_clockwise_position(i: int) -> tuple[int, int, int]:
    """Grid offset (in grid units, not nm) of molecule `i` on the starting grid.

    Same spiral layout the generated prepare.py's build_sim() uses to place
    molecules, exposed here so previews (e.g. the analyze.ipynb planning
    cell) can reproduce the exact starting layout without duplicating it.
    """
    x, y = 0, 0
    dx, dy = 0, -1

    for _ in range(i):
        if x == y or (x < 0 and x == -y) or (x > 0 and x == 1 - y):
            dx, dy = -dy, dx
        x, y = x + dx, y + dy

    return (x, y, 0)


PREPARE_TEMPLATE = '''import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md

# Parameters
def get_clockwise_position(i):
     x, y = 0, 0
     dx, dy = 0, -1  # Starting direction logic

     for _ in range(i):
          if x == y or (x < 0 and x == -y) or (x > 0 and x == 1 - y):
               # Rotate direction 90 degrees clockwise: (dx, dy) -> (-dy, dx)
               dx, dy = -dy, dx
          x, y = x + dx, y + dy

     return np.array((x, y, 0))

def build_sim(sim: Sim):
     components = sim.components

     ibead = 0
     i = 0
     for comp in components:
          for idx in range(comp.nmol):
               j = ibead + comp.nbeads

               x_c = (sim.box * 0.5)
               x_c[2] = 2

               pos = comp.xinit + x_c + get_clockwise_position(i) * __SPACING__
               sim.pos[ibead:j] = pos
               ibead = j
               i += 1

     md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)


# Job settings for Delft Blue
partition = "__PARTITION__"
runtime = "__WALLTIME__"
cpu_per_task = "__CPU_PER_TASK__"

sim_name = Path(__file__).parent.name

box = [__L__, __L__, __Z_HEIGHT__]
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
    l_box, clamped, spacing, min_sane_spacing, spacing_too_tight = box_plan

    n_save, n_frames, actual_steps = plan_steps(steps, max_n_save=max_n_save)

    content = (
        PREPARE_TEMPLATE
        .replace("__SPACING__", str(spacing))
        .replace("__L__", str(l_box))
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
    typer.echo(f"   box:        [{l_box}, {l_box}, {Z_HEIGHT}] nm")
    if clamped:
        achieved = round(nmol / (l_box * l_box), 5)
        typer.echo(
            f"   note:       requested concentration needed a {round(math.sqrt(nmol / concentration), 2)} nm box, "
            f"below the {MIN_BOX_L} nm floor required by CALVADOS's default cutoffs — "
            f"clamped to {l_box} nm (achieved concentration: {achieved} chains/nm^2)"
        )
    typer.echo(f"   spacing:    {spacing} nm between starting positions")
    if spacing_too_tight:
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
