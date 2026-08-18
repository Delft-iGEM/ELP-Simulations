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
from typing import Annotated

import typer

app = typer.Typer(add_completion=False, no_args_is_help=False)

# Fixed physical constants shared with the existing surface-attached templates.
Z_HEIGHT = 50.0  # nm, box height
Z_WALL = 1.9  # nm, depth of the surface-attachment potential well


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
partition = "gpu-a100"
runtime = "24:30:00"
cpu_per_task = "18"

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
          platform = 'CUDA',
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
) -> Path:
    """Scaffold simulations/<name>/prepare.py and run it. Returns the sim folder."""

    if concentration <= 0:
        raise typer.BadParameter("Concentration must be greater than 0.")
    if nmol < 1:
        raise typer.BadParameter("Number of ELPs must be at least 1.")
    if steps < 1:
        raise typer.BadParameter("Number of steps must be at least 1.")

    seq = _validate_sequence(sequence)
    slug = _slugify(name)

    root = _project_root()
    sim_dir = root / "simulations" / slug
    if sim_dir.exists():
        raise typer.BadParameter(
            f"simulations/{slug}/ already exists — choose a different name."
        )

    # Surface footprint sized so nmol molecules sit at the requested
    # concentration (chains / nm^2): L^2 = nmol / concentration.
    l_box = math.sqrt(nmol / concentration)

    # OpenMM requires the nonbonded cutoff to be less than half the box size.
    # CALVADOS's default electrostatics cutoff (cutoff_yu) is 4.0 nm, so the
    # box footprint must be kept above this floor even if the requested
    # concentration/nmol would otherwise call for a tighter box.
    min_l = 10.0
    clamped = l_box < min_l
    l_box = round(max(l_box, min_l), 2)

    # Average spacing between molecules, recomputed from the final (possibly
    # clamped) box so the initial grid always fits inside it.
    spacing = round(l_box / math.sqrt(nmol), 3)

    # Rough lower bound on a chain's own footprint (random-walk end-to-end
    # estimate at the ~0.38 nm CALVADOS bond length). If molecules start
    # packed tighter than this, they'll begin heavily overlapped and the
    # simulation is very likely to diverge to NaN once dynamics starts.
    min_sane_spacing = 0.38 * math.sqrt(len(seq))
    spacing_too_tight = spacing < min_sane_spacing

    # Pick a save frequency that yields ~1000 saved frames, capped at the
    # 7000 used by the existing templates, and derive the actual step count
    # (kept as an exact multiple of N_save, same convention as prepare.py).
    n_save = min(7000, max(1, steps // 1000))
    n_frames = max(1, round(steps / n_save))
    actual_steps = n_save * n_frames

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
    )

    sim_dir.mkdir(parents=True)
    (sim_dir / "prepare.py").write_text(content)

    typer.echo(f"✓  Created simulations/{slug}/prepare.py")
    typer.echo(f"   box:        [{l_box}, {l_box}, {Z_HEIGHT}] nm")
    if clamped:
        achieved = round(nmol / (l_box * l_box), 5)
        typer.echo(
            f"   note:       requested concentration needed a {round(math.sqrt(nmol / concentration), 2)} nm box, "
            f"below the {min_l} nm floor required by CALVADOS's default cutoffs — "
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
