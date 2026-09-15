"""
sim new — scaffold a new surface-attached ELP simulation.

Prompts for an ELP sequence, how far apart to graft the chains (an absolute
surface concentration, a spacing given as a fraction of the ELP's own length,
or a grafted mass per area), how many copies to simulate, and how many steps to
run, then creates
simulations/<name>/prepare.py (in the same style as
simulations/surfaceattached-multiple) and immediately runs it to produce a
ready runtime/ folder.

Usage:
    uv run python -m tools.new_simulation
    uv run sim new

Or non-interactively:
    uv run sim new --sequence VPGIGVPGIGVPGIG... --concentration 0.02 \\
        --nmol 500 --steps 7070000 --name my-elp-surface

    # spacing that scales with the ELP: half its fully-extended length, so a
    # 300 nm ELP is grafted 150 nm apart and a 600 nm one 300 nm apart
    uv run sim new --sequence VPGIGVPGIGVPGIG... --spacing-fraction 0.5 \\
        --nmol 500 --steps 7070000 --name my-elp-surface

    # same amount of polymer per area whatever the chain length
    uv run sim new --sequence VPGIGVPGIGVPGIG... --mass-concentration 0.2 \\
        --nmol 500 --steps 7070000 --name my-elp-surface

    # opt in to extra empty surface around the lattice, so chains stop wrapping
    # across the periodic boundary (off by default — the lattice tiles exactly)
    uv run sim new ... --margin-fraction 0.5     # half a spacing of rim per side

    # opt in to reactive crosslinking: lysine pairs within 0.8 nm bond during the
    # run (off by default; see tools/crosslink.py for the caveats)
    uv run sim new ... --crosslink-distance 0.8 --crosslink-start-step 1000000
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Annotated, NamedTuple

import typer

from tools.crosslink import CrosslinkSettings

app = typer.Typer(add_completion=False, no_args_is_help=False)

# Fixed physical constants shared with the existing surface-attached templates.
Z_HEIGHT = 50.0  # nm, box height
Z_WALL = 1.9  # nm, depth of the surface-attachment potential well

# OpenMM requires the nonbonded cutoff to be less than half the box size.
# CALVADOS's default electrostatics cutoff (cutoff_yu) is 4.0 nm, so the box
# footprint must be kept above this floor even if the requested
# concentration/nmol would otherwise call for a tighter box.
MIN_BOX_L = 10.0  # nm

# Empty rim left between the outermost grafting points and the box wall, as a
# fraction of the lattice spacing (one margin on *each* side). Off by default:
# with no rim the lattice tiles its own periodic images seamlessly, so every
# chain — rim or middle — sits exactly one spacing from its neighbours and the
# run is a uniform infinite grafted surface. Ask for a rim (0.5 = the edge
# chains a full spacing from the wall) to keep chains away from the periodic
# boundary so they stop wrapping to the far side of the trajectory, at the cost
# of that uniformity. See plan_box for the trade-off.
DEFAULT_MARGIN_FRACTION = 0.0

# CALVADOS's CA-CA bond length. Sets both the fully-extended ("contour") length
# of a chain and the random-walk ("coil") estimate of its footprint.
BOND_L = 0.38  # nm

LENGTH_MEASURES = ("contour", "coil")

# For the mass-per-area spacing mode: one water is released per peptide bond, so
# a chain weighs the sum of its residue masses plus one water.
WATER_DA = 18.015
AVOGADRO = 6.02214076e23




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


def residue_masses(csv_path: Path | str | None = None) -> dict[str, float]:
    """one-letter code -> residue molar mass in Da, with the anchor tag resolved.

    "Z" is CALVADOS's surface-anchor tag and the table gives it MW -2, so that
    the +2 N-terminus patch cancels it to exactly zero and OpenMM pins the bead.
    That is a simulation trick, not chemistry: the row's three-letter code says
    which amino acid the bead really is (VAL), and that is what the polymer
    actually weighs, so "Z" is counted as that residue here.
    """
    path = Path(csv_path) if csv_path else _project_root() / "residues_CALVADOS2.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    by_three = {row["three"]: float(row["MW"]) for row in rows if float(row["MW"]) > 0}
    masses = {}
    for row in rows:
        mw = float(row["MW"])
        masses[row["one"]] = mw if mw > 0 else by_three.get(row["three"], 0.0)
    return masses


def chain_molar_mass(sequence: str, csv_path: Path | str | None = None) -> float:
    """Molar mass of one chain in Da (g/mol) — its residues plus one water."""
    if not sequence:
        return 0.0
    masses = residue_masses(csv_path)
    return sum(masses.get(letter, 0.0) for letter in sequence) + WATER_DA


def mass_per_area(chain_mass_da: float, spacing_nm: float) -> float:
    """Grafted mass per unit area in ug/cm^2 — one chain per spacing^2 of surface.

    This is the quantity a QCM-D or ellipsometry measurement reports, and the
    one to hold fixed when comparing ELPs of different lengths: the same
    ug/cm^2 means the same amount of polymer on the surface, whether that is
    many short chains or few long ones.
    """
    if spacing_nm <= 0 or chain_mass_da <= 0:
        return 0.0
    # g per chain / nm^2 -> ug/cm^2: x 1e14 nm^2/cm^2 x 1e6 ug/g.
    return chain_mass_da * 1e20 / (AVOGADRO * spacing_nm ** 2)


def spacing_for_mass_per_area(chain_mass_da: float, ug_per_cm2: float) -> float:
    """Inverse of mass_per_area: the spacing that puts `ug_per_cm2` on the surface."""
    return math.sqrt(chain_mass_da * 1e20 / (AVOGADRO * ug_per_cm2))


def chain_length_nm(n_residues: int, measure: str = "contour") -> float:
    """The chain's own length scale in nm — what `spacing_fraction` is a fraction of.

    "contour" (default): the fully-extended length of the molecule, (N-1) bonds
    at the CALVADOS bond length. A 790-residue ELP is ~300 nm, a 1580-residue
    one ~600 nm, so spacing_fraction=0.5 puts their anchors 150 nm and 300 nm
    apart respectively — the spacing scales with the ELP instead of being the
    same absolute number for every run.

    "coil": the random-walk end-to-end size, 0.38*sqrt(N) — the same rough
    footprint estimate `min_sane_spacing` uses (~11 nm for that 790-mer).
    """
    if measure == "contour":
        return BOND_L * max(1, n_residues - 1)
    if measure == "coil":
        return BOND_L * math.sqrt(n_residues)
    raise typer.BadParameter(
        f"length_measure must be one of {', '.join(LENGTH_MEASURES)} — got {measure!r}."
    )


class BoxPlan(NamedTuple):
    l_box_x: float
    l_box_y: float
    nx: int
    ny: int
    spacing: float
    clamped: bool
    # Grafting density of the lattice itself, 1/spacing^2 — the density the
    # chains actually feel from their neighbours, and what the concentration /
    # mass knobs are asked in. With a margin the box is bigger than the lattice,
    # so the *box-averaged* density (box_concentration below) is lower.
    achieved_concentration: float
    min_sane_spacing: float
    spacing_too_tight: bool
    aspect_ratio: float
    # The chain's own length scale (chain_length_nm) and the fraction of it the
    # spacing actually came out at — meaningful in both modes, so a
    # concentration-driven plan still reports how the spacing compares to the
    # size of the chain.
    chain_length: float
    achieved_spacing_fraction: float
    length_measure: str
    # Molar mass of one chain and the grafted mass per area it works out to.
    # Like the fraction above, reported in every mode, so a run planned by
    # concentration still says how much polymer per cm^2 that puts down.
    chain_mass_da: float
    achieved_mass_concentration: float
    # Empty rim between the outermost grafting points and the box wall, on each
    # side, and the fraction of the spacing it was asked for as.
    margin: float
    margin_fraction: float
    # Densities averaged over the whole box, margin included — below the lattice
    # values above whenever margin > 0.
    box_concentration: float
    box_mass_concentration: float


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


def smallest_untilted_nmol(
    target_spacing: float,
    nmol: int,
    search: int = 400,
    margin_fraction: float = DEFAULT_MARGIN_FRACTION,
) -> int | None:
    """Smallest count >= nmol whose box clears MIN_BOX_L at the requested spacing.

    Below that the short side of the box would fall under the cutoff floor and
    plan_box has to spread the chains out, quietly changing the spacing (and so
    the density); going up to this many molecules keeps the spacing exactly as
    asked for instead. The margin counts towards the box size, same as in
    plan_box.
    """
    for candidate in range(nmol, nmol + search + 1):
        nx, ny = plan_lattice(candidate)
        short_side = (min(nx, ny) + 2 * margin_fraction) * target_spacing
        if short_side >= MIN_BOX_L and max(nx, ny) / min(nx, ny) <= 2.0:
            return candidate
    return None


def plan_box(
    sequence: str | int,
    concentration: float | None,
    nmol: int,
    spacing_fraction: float | None = None,
    length_measure: str = "contour",
    mass_concentration: float | None = None,
    margin_fraction: float = DEFAULT_MARGIN_FRACTION,
) -> BoxPlan:
    """Periodic box + grafting lattice for `nmol` chains.

    `sequence` is the ELP sequence, or just its residue count when the molar
    mass is not needed (everything except the mass mode works from the count
    alone).

    How far apart the chains are grafted is set in one of three mutually
    exclusive ways — pass exactly one:

    * `concentration` (chains/nm^2): an absolute spacing, 1/sqrt(concentration).
      The same number of chains per area for every run, whatever the ELP.
    * `spacing_fraction`: the spacing as a fraction of the chain's *own* length
      (`chain_length_nm(n_residues, length_measure)`), so it scales with the
      ELP. At spacing_fraction=0.5 a 300 nm ELP gets 150 nm between anchors and
      a 600 nm one gets 300 nm.
    * `mass_concentration` (ug/cm^2): the same *mass* of polymer per unit area
      for every run. A chain twice as long weighs twice as much, so it is
      grafted sqrt(2) times further apart — which is what to hold fixed when
      comparing against an experiment that measures adsorbed mass.

    Shared by `sim new` and the analyze.ipynb planning cell so both agree on
    the same box math.

    The chains are laid out on an exactly-filled nx-by-ny lattice, and
    `margin_fraction` sets how much empty surface is left between the outermost
    grafting points and the box wall — `margin = margin_fraction * spacing` on
    each side, so the box is nx*spacing + 2*margin by ny*spacing + 2*margin.

    The default is no margin, and that is the physically clean choice: the box
    is then exactly nx*spacing by ny*spacing, the periodic images continue the
    lattice without a seam, and a rim chain — half a spacing from the wall, half
    a spacing from the wall on the other side — sits exactly one spacing from
    its neighbour across the boundary, as crowded as one in the middle. The run
    is a uniform infinite grafted surface.

    A margin trades that uniformity for a more readable trajectory. A chain that
    leans past the boundary is written wrapped to the opposite side of the box,
    which looks like it teleporting; the dynamics are correct and the analyze
    notebook unwraps it, but it is much easier to watch if it rarely happens at
    all. margin_fraction=0.5 puts the edge chains a full spacing from the wall,
    which is what stops most of the wrapping. The cost is that neighbours across
    the boundary are then spacing + 2*margin apart instead of spacing, so the
    rim is less crowded than the middle: the run becomes a finite grafted patch
    repeated periodically rather than a strictly uniform surface. Worth it for
    per-chain observables (Rg, height, RMSD); not worth it for anything that
    depends on uniform lateral crowding.

    `spacing` always means the nearest-neighbour distance on the lattice, so the
    concentration / mass knobs keep their meaning whatever the margin; the
    margin dilutes the *box-averaged* density, reported separately as
    `box_concentration` / `box_mass_concentration`.
    """
    given = [concentration, spacing_fraction, mass_concentration]
    if sum(value is not None for value in given) != 1:
        raise typer.BadParameter(
            "Give exactly one of concentration (chains/nm^2), spacing_fraction "
            "(fraction of the chain's own length) or mass_concentration (ug/cm^2)."
        )

    if margin_fraction < 0:
        raise typer.BadParameter("margin_fraction cannot be negative.")

    n_residues = len(sequence) if isinstance(sequence, str) else int(sequence)
    chain_mass = chain_molar_mass(sequence) if isinstance(sequence, str) else 0.0

    nx, ny = plan_lattice(nmol)
    aspect_ratio = max(nx, ny) / min(nx, ny)

    chain_length = chain_length_nm(n_residues, length_measure)

    if spacing_fraction is not None:
        # Spacing relative to the molecule itself, so a longer ELP is grafted
        # proportionally further from its neighbours instead of every run
        # sharing one absolute spacing.
        if spacing_fraction <= 0:
            raise typer.BadParameter("spacing_fraction must be greater than 0.")
        spacing_raw = spacing_fraction * chain_length
    elif mass_concentration is not None:
        # Same mass of polymer per unit area whatever the chain length: one
        # chain of chain_mass per spacing^2 of surface.
        if mass_concentration <= 0:
            raise typer.BadParameter("mass_concentration must be greater than 0.")
        if chain_mass <= 0:
            raise typer.BadParameter(
                "mass_concentration needs the sequence itself (its molar mass), "
                "not just a residue count."
            )
        spacing_raw = spacing_for_mass_per_area(chain_mass, mass_concentration)
    else:
        # The lattice spacing *is* the surface concentration: one chain per
        # spacing^2 of surface, so spacing = 1/sqrt(concentration).
        if concentration <= 0:
            raise typer.BadParameter("Concentration must be greater than 0.")
        spacing_raw = 1.0 / math.sqrt(concentration)

    # OpenMM needs the nonbonded cutoff below half the box, so the *short* side
    # of the box has to clear MIN_BOX_L. With nmol fixed the only way to widen
    # the box is to spread the chains further apart, which lowers the
    # concentration (raises the achieved fraction) — reported loudly by the
    # caller. The margin counts towards the floor, so a generous margin can be
    # what gets a small run over it.
    short_side_cells = min(nx, ny) + 2 * margin_fraction
    scale = max(1.0, MIN_BOX_L / (short_side_cells * spacing_raw))
    clamped = scale > 1.0
    spacing = round(spacing_raw * scale, 3)

    # The empty rim on each side. Rounded like the spacing so the numbers
    # written into prepare.py reproduce the box exactly.
    margin = round(margin_fraction * spacing, 3)

    l_box_x = round(nx * spacing + 2 * margin, 6)
    l_box_y = round(ny * spacing + 2 * margin, 6)
    # The density the chains feel is set by the lattice; the box average is what
    # the margin dilutes it to. They are the same number when margin == 0.
    achieved_concentration = 1.0 / spacing ** 2
    box_concentration = nmol / (l_box_x * l_box_y)

    # Rough lower bound on a chain's own footprint (random-walk end-to-end
    # estimate at the CALVADOS bond length). If molecules start
    # packed tighter than this, they'll begin heavily overlapped and the
    # simulation is very likely to diverge to NaN once dynamics starts.
    min_sane_spacing = BOND_L * math.sqrt(n_residues)
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
        chain_length,
        spacing / chain_length,
        length_measure,
        chain_mass,
        mass_per_area(chain_mass, spacing),
        margin,
        margin_fraction,
        box_concentration,
        mass_per_area(chain_mass, spacing) * box_concentration / achieved_concentration,
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


def lattice_positions(nmol: int, spacing: float, margin: float = 0.0) -> list[tuple[float, float]]:
    """Grafting-point (x, y) of every molecule, in nm, on the lattice.

    Row-major over the nx-by-ny lattice from plan_lattice, each molecule at the
    *centre* of its spacing-by-spacing cell, with the whole lattice shifted in
    by `margin` so the empty rim plan_box left is on all four sides. The
    generated prepare.py's build_sim() computes exactly the same positions, so
    previews (e.g. the analyze.ipynb planning cell) reproduce the real starting
    layout.
    """
    nx, _ny = plan_lattice(nmol)
    return [
        (margin + (i % nx + 0.5) * spacing, margin + (i // nx + 0.5) * spacing)
        for i in range(nmol)
    ]


PREPARE_TEMPLATE = '''import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md
from tools.paths import ensure_runtime_dir
from tools.metadata import write_lattice, write_metadata
from tools.crosslink import CrosslinkSettings, write_settings

# Grafting lattice: __NX__ x __NY__ points, one chain per spacing x spacing
# cell, with an empty `margin` rim between the outermost points and the box
# wall on every side (box = nx*spacing + 2*margin). The rim keeps chains away
# from the periodic boundary, so they rarely lean across it and get wrapped to
# the far side of the trajectory; the cost is that neighbours across the
# boundary sit spacing + 2*margin apart, so rim chains are a little less
# crowded than the ones in the middle. margin = 0 restores the seamless tiling.
nx, ny = __NX__, __NY__
# __SPACING_NOTE__
spacing = __SPACING__
# __MARGIN_NOTE__
margin = __MARGIN__

# Reactive crosslinking. None = off, and then run.py takes the stock CALVADOS
# path unchanged. A dict turns it on — see tools/crosslink.py for what the
# numbers mean and, just as importantly, what they don't.
# __CROSSLINK_NOTE__
crosslink = __CROSSLINK__


def build_sim(sim: Sim):
     components = sim.components

     # Each chain goes at the *centre* of its own cell (the +0.5), offset by the
     # margin — exactly what tools.new_simulation.lattice_positions previews.

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

               x0 = margin + (i % nx + 0.5) * spacing
               y0 = margin + (i // nx + 0.5) * spacing

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

     runtime_dir = ensure_runtime_dir(path)

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

     # The box alone doesn't say where the lattice inside it sits (spacing and
     # margin are two unknowns in one number), so record them next to the
     # trajectory for metadata.csv to pick up.
     write_lattice(runtime_dir, nx=nx, ny=ny, spacing=spacing, margin=margin)

     # run.py reads this file to decide whether the run reacts. Written (or
     # removed) every time, so a regenerated non-reactive run can never inherit
     # a stale crosslink.yaml from a previous attempt.
     write_settings(runtime_dir,
                    CrosslinkSettings.from_dict(crosslink) if crosslink else None)

     # Settings snapshot next to the trajectory-to-be. run.py rewrites it when
     # the run finishes (status/frame count); `sim metadata` refreshes it any time.
     write_metadata(runtime_dir)
'''


def create_simulation(
    sequence: str,
    concentration: float | None,
    nmol: int,
    steps: int,
    name: str,
    platform: str = "CUDA",
    partition: str = "gpu-a100",
    walltime: str = "24:30:00",
    cpu_per_task: str = "18",
    max_n_save: int = 7000,
    spacing_fraction: float | None = None,
    length_measure: str = "contour",
    mass_concentration: float | None = None,
    margin_fraction: float = DEFAULT_MARGIN_FRACTION,
    crosslink_distance: float | None = None,
    crosslink_valence: int = 1,
    crosslink_prob: float = 1.0,
    crosslink_check_every: int = 1000,
    crosslink_k: float = 2000.0,
    crosslink_r0: float = 0.6,
    crosslink_ramp_steps: int = 500,
    crosslink_selection: list[int] | None = None,
    crosslink_start_step: int = 0,
    crosslink_seed: int | None = None,
) -> Path:
    """Scaffold simulations/<name>/prepare.py and run it. Returns the sim folder.

    Pass exactly one of `concentration` (chains/nm^2), `spacing_fraction`
    (fraction of this ELP's own length) or `mass_concentration` (ug/cm^2, the
    same grafted mass per area whatever the chain length) — see plan_box.

    `margin_fraction` is the empty rim left around the lattice, as a fraction of
    the spacing (see plan_box). It defaults to 0 — a seamless, uniformly crowded
    lattice; raise it to keep chains off the periodic boundary so they stop
    wrapping to the far side of the trajectory.

    `crosslink_distance` (nm) turns on reactive crosslinking: lysine pairs that
    come within it bond permanently *during* the run. None — the default — leaves
    the run non-reactive and bit-for-bit what it was before the feature existed.
    The other crosslink_* arguments only matter when it is set; see
    tools/crosslink.py for what they mean and how far the results can be trusted.
    """

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

    box_plan = plan_box(seq, concentration, nmol, spacing_fraction, length_measure,
                        mass_concentration, margin_fraction)
    spacing = box_plan.spacing
    margin = box_plan.margin
    min_sane_spacing = box_plan.min_sane_spacing

    n_save, n_frames, actual_steps = plan_steps(steps, max_n_save=max_n_save)

    # Validated here rather than at run time: a typo in a crosslink setting
    # should stop you now, not eight hours into a GPU job.
    crosslink_settings = None
    if crosslink_distance is not None:
        crosslink_settings = CrosslinkSettings(
            distance=crosslink_distance,
            valence=crosslink_valence,
            prob=crosslink_prob,
            check_every=crosslink_check_every,
            k=crosslink_k,
            r0=crosslink_r0,
            ramp_steps=crosslink_ramp_steps,
            selection=tuple(crosslink_selection) if crosslink_selection else None,
            start_step=crosslink_start_step,
            seed=crosslink_seed,
        )

    if spacing_fraction is not None:
        spacing_note = (f"spacing = {spacing_fraction} x the chain's own {length_measure} length "
                        f"({box_plan.chain_length:.2f} nm for these {len(seq)} residues)")
    elif mass_concentration is not None:
        spacing_note = (f"spacing = one {box_plan.chain_mass_da / 1000:.2f} kDa chain per "
                        f"spacing^2, i.e. {mass_concentration} ug/cm^2 of grafted polymer")
    else:
        spacing_note = f"spacing = 1/sqrt({concentration} chains/nm^2)"
    if box_plan.clamped:
        spacing_note += f", widened to {spacing} nm to clear the {MIN_BOX_L} nm box floor"
    if margin:
        margin_note = (f"margin = {margin_fraction} x the {spacing} nm spacing, empty on every "
                       f"side, so chains start {margin} nm clear of the periodic boundary")
    else:
        margin_note = ("margin = 0: the lattice fills the box exactly, so the periodic images "
                       "continue it seamlessly and every chain is equally crowded")

    if crosslink_settings is None:
        crosslink_note = "crosslinking off (crosslink_distance = None)"
        crosslink_literal = "None"
    else:
        crosslink_note = (
            f"lysine pairs within {crosslink_settings.distance} nm "
            f"({crosslink_settings.distance_angstrom} A) bond during the run, "
            f"valence {crosslink_settings.valence}, p={crosslink_settings.prob}")
        crosslink_literal = repr(crosslink_settings.to_dict())

    content = (
        PREPARE_TEMPLATE
        .replace("__SPACING_NOTE__", spacing_note)
        .replace("__SPACING__", str(spacing))
        .replace("__MARGIN_NOTE__", margin_note)
        .replace("__MARGIN__", str(margin))
        .replace("__CROSSLINK_NOTE__", crosslink_note)
        .replace("__CROSSLINK__", crosslink_literal)
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
    typer.echo(f"   lattice:    {box_plan.nx} x {box_plan.ny} grafting points, {spacing} nm apart")
    if margin:
        typer.echo(f"   margin:     {margin} nm of empty surface on every side "
                   f"({margin_fraction} x the spacing) — chains this far from the periodic "
                   f"boundary, so neighbours across it are {round(spacing + 2 * margin, 3)} nm apart")
    else:
        typer.echo("   margin:     none (default) — the lattice tiles the box exactly, so the "
                   "periodic images continue it seamlessly and every chain is equally "
                   "crowded; chains do reach the boundary and wrap")
    typer.echo(f"   chain:      {len(seq)} residues, {box_plan.chain_length:.2f} nm "
               f"({length_measure} length) — spacing is "
               f"{box_plan.achieved_spacing_fraction:.3f} x that")
    typer.echo(f"   density:    {round(box_plan.achieved_concentration, 5)} chains/nm^2 on the lattice"
               + (f", {round(box_plan.box_concentration, 5)} averaged over the box "
                  f"(the margin dilutes it)" if margin else ""))
    typer.echo(f"   mass:       {box_plan.chain_mass_da / 1000:.2f} kDa per chain, "
               f"{box_plan.achieved_mass_concentration:.4f} ug/cm^2 grafted "
               f"({box_plan.achieved_mass_concentration * 10:.3f} mg/m^2)"
               + (f"; {box_plan.box_mass_concentration:.4f} ug/cm^2 over the whole box"
                  if margin else ""))
    if box_plan.clamped:
        if spacing_fraction is not None:
            requested_spacing = spacing_fraction * box_plan.chain_length
            asked_for = (f"a spacing of {spacing_fraction} x the chain's "
                         f"{box_plan.chain_length:.2f} nm {length_measure} length "
                         f"({requested_spacing:.2f} nm)")
        elif mass_concentration is not None:
            requested_spacing = spacing_for_mass_per_area(box_plan.chain_mass_da, mass_concentration)
            asked_for = f"{mass_concentration} ug/cm^2 ({requested_spacing:.2f} nm apart)"
        else:
            requested_spacing = 1.0 / math.sqrt(concentration)
            asked_for = f"{concentration} chains/nm^2 ({requested_spacing:.2f} nm apart)"
        typer.echo(
            f"⚠  the requested {asked_for} would make the short side of the box "
            f"{round((min(box_plan.nx, box_plan.ny) + 2 * margin_fraction) * requested_spacing, 2)} nm "
            f"(lattice + margin), below the "
            f"{MIN_BOX_L} nm floor required by CALVADOS's default cutoffs — the chains were spread "
            f"out to {spacing} nm instead ({round(box_plan.achieved_concentration, 5)} chains/nm^2, "
            f"{box_plan.achieved_mass_concentration:.4f} ug/cm^2, "
            f"{box_plan.achieved_spacing_fraction:.3f} x the chain length)."
        )
        bigger = smallest_untilted_nmol(requested_spacing, nmol, margin_fraction=margin_fraction)
        if bigger is not None:
            bx, by = plan_lattice(bigger)
            typer.echo(
                f"   fix:        use nmol={bigger} ({bx} x {by}) to get the full "
                f"{requested_spacing:.2f} nm spacing with no clamping."
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
            f"Consider a lower concentration / mass loading, a larger spacing fraction, "
            f"and/or fewer molecules."
        )
    typer.echo(f"   molecules:  {nmol}")
    if crosslink_settings is not None:
        typer.echo(f"   crosslink:  ON — lysine pairs within "
                   f"{crosslink_settings.distance} nm ({crosslink_settings.distance_angstrom} A) "
                   f"bond during the run; valence {crosslink_settings.valence}, "
                   f"p={crosslink_settings.prob}, checked every "
                   f"{crosslink_settings.check_every} steps")
        if crosslink_settings.start_step:
            typer.echo(f"               reactions ignored before step "
                       f"{crosslink_settings.start_step} "
                       f"({crosslink_settings.start_step * 0.01 / 1000:.3f} ns)")
        else:
            typer.echo("⚠  crosslink_start_step is 0, so contacts left over from the initial "
                       "placement count as reactions. Set it past equilibration — the analyze "
                       "notebook's Rg/RMSD check says where that is.")
    else:
        typer.echo("   crosslink:  off (default) — nothing reacts, identical to the "
                   "pre-crosslinking pipeline")
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
    concentration: Annotated[
        float | None,
        typer.Option(help="Surface concentration (ELP chains / nm^2) — absolute spacing."),
    ] = None,
    spacing_fraction: Annotated[
        float | None,
        typer.Option(
            help="Spacing between grafting points as a fraction of the chain's own length "
                 "(0.5 => a 300 nm ELP is grafted 150 nm apart, a 600 nm one 300 nm apart).",
        ),
    ] = None,
    mass_concentration: Annotated[
        float | None,
        typer.Option(
            help="Grafted mass per area (ug/cm^2) — holds the amount of polymer on the "
                 "surface fixed, so a chain twice as long is grafted sqrt(2) further apart.",
        ),
    ] = None,
    length_measure: Annotated[
        str,
        typer.Option(help="Which chain length spacing-fraction is a fraction of: "
                          "'contour' (fully extended) or 'coil' (random-walk size)."),
    ] = "contour",
    margin_fraction: Annotated[
        float,
        typer.Option(
            help="Empty surface left between the outermost chains and the box wall, as a "
                 "fraction of the spacing, on each side. Default 0 = the lattice tiles the "
                 "box exactly, so every chain is equally crowded (but chains reach the "
                 "boundary and wrap). 0.5 keeps chains a full spacing clear of the boundary "
                 "so they stop wrapping, at the cost of a less crowded rim.",
        ),
    ] = DEFAULT_MARGIN_FRACTION,
    crosslink_distance: Annotated[
        float | None,
        typer.Option(
            help="Turn on reactive crosslinking: lysine pairs that come within this distance "
                 "IN NANOMETRES bond permanently during the run (0.8 nm = 8 A). Default off, "
                 "which leaves the run bit-for-bit identical to the non-reactive pipeline.",
        ),
    ] = None,
    crosslink_valence: Annotated[
        int,
        typer.Option(help="Max bonds per lysine. 1 = bifunctional crosslinker at 1:1 "
                          "(glutaraldehyde-like); 2 = trifunctional junctions (THPP-like)."),
    ] = 1,
    crosslink_prob: Annotated[
        float,
        typer.Option(help="Probability a within-cutoff pair reacts at a given check. "
                          "1.0 = diffusion-limited; below that the reaction rate is "
                          "decoupled from the diffusion rate (sweep for Damkohler)."),
    ] = 1.0,
    crosslink_check_every: Annotated[
        int,
        typer.Option(help="MD steps between reaction checks."),
    ] = 1000,
    crosslink_k: Annotated[
        float,
        typer.Option(help="Final crosslink stiffness, kJ/mol/nm^2."),
    ] = 2000.0,
    crosslink_r0: Annotated[
        float,
        typer.Option(help="Crosslink equilibrium length, nm."),
    ] = 0.6,
    crosslink_ramp_steps: Annotated[
        int,
        typer.Option(help="Steps over which a new bond's stiffness rises from 0 to the full "
                          "value. Ramping avoids the energy spike of switching a stiff bond "
                          "on instantaneously."),
    ] = 500,
    crosslink_start_step: Annotated[
        int,
        typer.Option(help="Ignore reactions before this step. Set it past equilibration: "
                          "contacts in the early frames are artefacts of the initial "
                          "placement, not encounters the dynamics produced."),
    ] = 0,
    crosslink_selection: Annotated[
        list[int] | None,
        typer.Option(
            help="Explicit reactive bead indices (0-based), repeat the flag per bead. "
                 "Default: auto-detect lysines. Beads at residue 0 of a chain are dropped "
                 "either way — those are the fixed surface anchors.",
        ),
    ] = None,
    crosslink_seed: Annotated[
        int | None,
        typer.Option(help="Seed for the reaction RNG (only used when --crosslink-prob < 1). "
                          "Defaults to the run's OpenMM seed, else a recorded random draw."),
    ] = None,
) -> None:
    """Scaffold and prepare a new surface-attached ELP simulation.

    Spacing comes from exactly one of --concentration, --spacing-fraction or
    --mass-concentration. Give none of them interactively and you're asked
    which one you want.
    """
    # Typer prompts at parse time, which can't ask for one option *or* another,
    # so the choice between the three is prompted for here instead.
    modes = [concentration, spacing_fraction, mass_concentration]
    if sum(value is not None for value in modes) > 1:
        raise typer.BadParameter(
            "Give only one of --concentration, --spacing-fraction or --mass-concentration."
        )
    if not any(value is not None for value in modes):
        typer.echo("How should the chains be spaced?")
        typer.echo("  1  surface concentration  (chains/nm^2 — same chain count per area)")
        typer.echo("  2  fraction of the chain's own length  (scales with the ELP)")
        typer.echo("  3  mass per area  (ug/cm^2 — same amount of polymer per area)")
        choice = typer.prompt("Choose", type=int, default=1)
        if choice == 2:
            length_measure = typer.prompt(
                f"Which chain length to measure against ({'/'.join(LENGTH_MEASURES)})",
                default=length_measure,
            )
            spacing_fraction = typer.prompt("Spacing as a fraction of that length", type=float)
        elif choice == 3:
            mass_concentration = typer.prompt("Grafted mass per area (ug/cm^2)", type=float)
        else:
            concentration = typer.prompt("Surface concentration (ELP chains / nm^2)", type=float)

    create_simulation(
        sequence,
        concentration,
        nmol,
        steps,
        name,
        spacing_fraction=spacing_fraction,
        length_measure=length_measure,
        mass_concentration=mass_concentration,
        margin_fraction=margin_fraction,
        crosslink_distance=crosslink_distance,
        crosslink_valence=crosslink_valence,
        crosslink_prob=crosslink_prob,
        crosslink_check_every=crosslink_check_every,
        crosslink_k=crosslink_k,
        crosslink_r0=crosslink_r0,
        crosslink_ramp_steps=crosslink_ramp_steps,
        crosslink_start_step=crosslink_start_step,
        crosslink_selection=crosslink_selection,
        crosslink_seed=crosslink_seed,
    )


if __name__ == "__main__":
    app()
