"""Per-run metadata — ``runtime/metadata.csv``, next to the trajectory.

Every finished run should be readable months later without having to
reverse-engineer its ``prepare.py``, so the settings that define the run are
written as a small CSV in the *same folder the .dcd lands in*
(``simulations/<sim>/runtime/``, i.e. on ``/scratch`` wherever the data root
is configured — see tools.paths).

Everything here is recovered from what CALVADOS itself writes into that folder
(``config.yaml``, ``components.yaml``, ``molecules.fasta``) plus the DCD header,
rather than from the generating ``prepare.py``. That means one code path serves
both new runs and every simulation that already exists: `sim metadata --all`
back-fills the lot.

The file is written three times over a run's life, each time with better
information:

1. by ``prepare.py``, as soon as ``runtime/`` exists (settings, status "not
   started"),
2. by ``template/run.py`` when ``sim.simulate()`` returns (status "completed",
   real frame count),
3. whenever you ask, via ``sim metadata <sim>``.

Layout is one row per setting — ``key,value,unit,description`` — so a single
run is readable by eye and a whole batch is one ``pd.concat`` away from a
comparison table:

    import pandas as pd
    from pathlib import Path
    runs = {p.parent.parent.name: pd.read_csv(p, index_col="key")["value"]
            for p in Path("simulations").glob("*/runtime/metadata.csv")}
    pd.DataFrame(runs).T
"""

from __future__ import annotations

import csv
import math
import re
import struct
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NamedTuple

from yaml import safe_load

from tools.crosslink import (
    CROSSLINK_FILENAME,
    EVENTS_FILENAME,
    load_settings as load_crosslink_settings,
)
from tools.new_simulation import (
    BOND_L,
    Z_HEIGHT,
    chain_length_nm,
    chain_molar_mass,
    mass_per_area,
    plan_lattice,
)
from tools.results import publish

METADATA_FILENAME = "metadata.csv"

# The grafting lattice, written by prepare.py. The box alone can't be decoded
# into a spacing once there is a margin around the lattice (two unknowns, one
# number), so prepare.py states both. Runs made before margins existed have no
# such file and are read back as a seamless, margin-free lattice — which is
# exactly what they are.
LATTICE_FILENAME = "lattice.yaml"

# CALVADOS's Langevin integrator uses a fixed 0.01 ps step (calvados/sim.py);
# it is not part of config.yaml, so it has to be recorded here.
DT_PS = 0.01


class Field(NamedTuple):
    key: str
    value: Any
    unit: str
    description: str


def _read_yaml(path: Path) -> dict:
    if not path.is_file():
        return {}
    with open(path) as f:
        return safe_load(f) or {}


def _read_fasta(path: Path) -> dict[str, str]:
    """{name: sequence} — the sequences as actually simulated (anchor "Z" included)."""
    if not path.is_file():
        return {}
    sequences: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            sequences[name] = ""
        elif name is not None:
            sequences[name] += line
    return sequences


def dcd_n_frames(path: Path) -> int | None:
    """Frames written so far, from the DCD header — no MDAnalysis/mdtraj needed.

    A DCD starts with a 4-byte Fortran record length (84), the tag "CORD", and
    then the control array whose first entry is NSET, the frame count. OpenMM's
    reporter rewrites NSET as it appends, so this is current even mid-run and
    costs one 12-byte read instead of parsing a multi-GB trajectory.
    """
    if not path.is_file():
        return None
    with open(path, "rb") as f:
        head = f.read(12)
    if len(head) < 12 or head[4:8] != b"CORD":
        return None
    endian = "<" if struct.unpack("<i", head[:4])[0] == 84 else ">"
    return struct.unpack(endian + "i", head[8:12])[0]


def _z_wall_from_expr(expr: str | None) -> float | str:
    """The wall position out of CALVADOS's ext_force_expr, e.g. step(1.9-z)*... -> 1.9."""
    if not expr:
        return ""
    match = re.search(r"step\(\s*([0-9.eE+-]+)\s*-\s*z\s*\)", expr)
    return float(match.group(1)) if match else ""


def _package_version(name: str) -> str:
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(name)
    except PackageNotFoundError:
        return ""


def _round(value: float, digits: int = 6) -> float:
    return round(value, digits)


def _count_rows(path: Path) -> int | None:
    """Data rows in a CSV, or None if it isn't there (i.e. the run didn't react)."""
    if not path.is_file():
        return None
    with open(path, newline="") as f:
        return sum(1 for _ in csv.DictReader(f))


def collect_metadata(runtime_dir: Path) -> list[Field]:
    """Every recorded setting for the run in `runtime_dir`, in file order.

    Missing inputs degrade to empty values rather than raising — a half-prepared
    or interrupted folder should still produce a readable metadata.csv.
    """
    runtime_dir = Path(runtime_dir)
    config = _read_yaml(runtime_dir / "config.yaml")
    components = _read_yaml(runtime_dir / "components.yaml")

    # CALVADOS's sysname is the simulation name *and* what the trajectory file is
    # named after. Take it from there rather than from the folder: on DelftBlue
    # runtime/ is a symlink onto /scratch, so the resolved path run.py passes in
    # has the data root as its parent, not simulations/<sim>.
    sim_name = str(config.get("sysname") or "").strip()
    if not sim_name:
        sim_name = runtime_dir.parent.name if runtime_dir.name == "runtime" else runtime_dir.name
    sequences = _read_fasta(runtime_dir / "molecules.fasta")

    defaults = components.get("defaults", {}) or {}
    system = components.get("system", {}) or {}

    # nmol lives in defaults unless a molecule overrides it; with several
    # molecule types the chains share the lattice, so the total is what the
    # box/spacing math is about.
    per_type = {
        name: int((entry or {}).get("nmol", defaults.get("nmol", 0)))
        for name, entry in system.items()
    } or {sim_name: int(defaults.get("nmol", 0))}
    nmol = sum(per_type.values())

    seq = "".join(sequences.get(name, "") for name in per_type) if len(per_type) == 1 else ""
    seq_names = "|".join(per_type)
    n_residues = len(seq)

    box = config.get("box") or [0.0, 0.0, Z_HEIGHT]
    l_x, l_y, l_z = (float(v) for v in box[:3])

    # The lattice the chains were grafted on. prepare.py writes it down
    # (lattice.yaml) because a box with a margin around the lattice can't be
    # decoded back into a spacing on its own. Without that file — every run made
    # before margins existed — the lattice filled the box exactly, so the
    # spacing is just the box over the grid.
    lattice = _read_yaml(runtime_dir / LATTICE_FILENAME)
    nx, ny = plan_lattice(nmol) if nmol else (0, 0)
    if lattice:
        nx = int(lattice.get("nx", nx) or nx)
        ny = int(lattice.get("ny", ny) or ny)
        spacing = _round(float(lattice.get("spacing_nm", 0.0) or 0.0), 6)
        margin = _round(float(lattice.get("margin_nm", 0.0) or 0.0), 6)
        # The stated lattice should rebuild the box it was written next to.
        lattice_matches_box = bool(nx) and all(
            abs(side - (n * spacing + 2 * margin)) <= 1e-3 * max(side, 1.0)
            for side, n in ((l_x, nx), (l_y, ny))
        )
    else:
        margin = 0.0
        spacing_x = l_x / nx if nx else 0.0
        spacing_y = l_y / ny if ny else 0.0
        # Both sides of the box are exact multiples of one spacing for anything
        # `sim new` generated back then. A run scaffolded by hand (or by an
        # older layout that didn't tile) can disagree, and then the spacing is a
        # best guess rather than a fact — say so instead of quietly averaging.
        lattice_matches_box = bool(nx) and abs(spacing_x - spacing_y) <= 1e-3 * max(spacing_x, 1.0)
        spacing = _round(0.5 * (spacing_x + spacing_y), 6)
    # Two densities once there is a margin: what the chains feel from their
    # lattice neighbours, and what the box averages out to over the empty rim.
    area = l_x * l_y
    concentration = 1.0 / spacing ** 2 if spacing else 0.0
    box_concentration = nmol / area if area else 0.0

    contour = chain_length_nm(n_residues, "contour") if n_residues else 0.0
    coil = chain_length_nm(n_residues, "coil") if n_residues else 0.0

    # Weigh the chain with the residue table the run itself used (components.yaml
    # records its path), falling back to the project's copy if that has moved.
    fresidues = Path(str(defaults.get("fresidues", "")))
    chain_mass = chain_molar_mass(seq, fresidues if fresidues.is_file() else None) if seq else 0.0
    mass_concentration = mass_per_area(chain_mass, spacing)
    box_mass_concentration = (
        mass_concentration * box_concentration / concentration if concentration else 0.0
    )

    # Reactive crosslinking, if the run had any. Absent crosslink.yaml = the run
    # was not reactive, which is the default and most runs.
    crosslink = load_crosslink_settings(runtime_dir)
    n_crosslinks = _count_rows(runtime_dir / EVENTS_FILENAME)

    steps = int(config.get("steps", 0) or 0)
    wfreq = int(config.get("wfreq", 0) or 0)
    n_frames_planned = steps // wfreq if wfreq else 0

    dcd = runtime_dir / f"{sim_name}.dcd"
    n_frames_saved = dcd_n_frames(dcd)
    if n_frames_saved is None:
        status = "not started"
    elif n_frames_planned and n_frames_saved >= n_frames_planned:
        status = "completed"
    else:
        status = "incomplete"

    return [
        Field("sim_name", sim_name, "", "Simulation folder under simulations/"),
        Field("metadata_written_utc", datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "", "When this file was last (re)generated"),
        Field("status", status, "", "not started / incomplete / completed, from the DCD frame count"),

        # ---- what was simulated ----
        Field("sequence", seq, "", "Sequence as simulated; residue 0 is 'Z', the grafted anchor bead"),
        Field("sequence_names", seq_names, "", "FASTA entry name(s) in molecules.fasta"),
        Field("n_residues", n_residues, "residues", "Chain length in residues"),
        Field("nmol", nmol, "chains", "Number of ELP chains in the box"),
        Field("chain_contour_length_nm", _round(contour, 3), "nm",
              f"Fully-extended chain length, (N-1) x {BOND_L} nm"),
        Field("chain_coil_size_nm", _round(coil, 3), "nm",
              f"Random-walk end-to-end estimate, {BOND_L} x sqrt(N)"),
        Field("chain_molar_mass_da", _round(chain_mass, 2), "Da",
              "Molar mass of one chain (residues + one water). The anchor bead is "
              "counted as the valine it is a tagged copy of, not as its MW = -2 "
              "simulation placeholder"),

        # ---- how they were grafted ----
        Field("spacing_nm", spacing, "nm", "Distance between neighbouring grafting points"),
        Field("spacing_fraction_contour", _round(spacing / contour, 4) if contour else "", "",
              "spacing_nm as a fraction of chain_contour_length_nm (the --spacing-fraction knob)"),
        Field("spacing_fraction_coil", _round(spacing / coil, 4) if coil else "", "",
              "spacing_nm as a fraction of chain_coil_size_nm"),
        Field("margin_nm", _round(margin, 6), "nm",
              "Empty surface between the outermost grafting points and the box wall, each side"),
        Field("margin_fraction", _round(margin / spacing, 4) if spacing else "", "",
              "margin_nm as a fraction of the spacing (the --margin-fraction knob). 0 = the "
              "lattice tiles the box exactly"),
        Field("periodic_neighbour_spacing_nm", _round(spacing + 2 * margin, 6), "nm",
              "Distance between chains on opposite edges, across the periodic boundary — "
              "spacing_nm + 2 x margin_nm, so equal to spacing_nm only when margin is 0"),
        Field("concentration_chains_per_nm2", _round(concentration, 6), "chains/nm^2",
              "Grafting density on the lattice, 1 / spacing_nm^2 — what the chains feel "
              "from their neighbours, and what --concentration asks for"),
        Field("box_concentration_chains_per_nm2", _round(box_concentration, 6), "chains/nm^2",
              "Density averaged over the whole box, nmol / (box_x * box_y) — below the "
              "lattice value by however much the margin dilutes it"),
        Field("mass_concentration_ug_cm2", _round(mass_concentration, 6), "ug/cm^2",
              "Grafted polymer mass per area (the --mass-concentration knob) — one "
              "chain of chain_molar_mass_da per spacing^2"),
        Field("mass_concentration_mg_m2", _round(mass_concentration * 10, 6), "mg/m^2",
              "The same quantity in the other common unit (1 ug/cm^2 = 10 mg/m^2)"),
        Field("box_mass_concentration_ug_cm2", _round(box_mass_concentration, 6), "ug/cm^2",
              "Grafted mass per area averaged over the whole box, margin included"),
        Field("lattice_matches_box", lattice_matches_box, "",
              "True = nx x ny cells plus the margins rebuild the box exactly, so spacing_nm "
              "is exact; False = spacing_nm is inferred and approximate"),
        Field("lattice_nx", nx, "", "Grafting points along x"),
        Field("lattice_ny", ny, "", "Grafting points along y"),
        Field("box_x_nm", l_x, "nm", "Periodic box, x"),
        Field("box_y_nm", l_y, "nm", "Periodic box, y"),
        Field("box_z_nm", l_z, "nm", "Periodic box, z (height above the surface)"),

        # ---- did it react? ----
        Field("crosslinking", bool(crosslink), "",
              f"True = lysines bonded during the run (settings in {CROSSLINK_FILENAME}); "
              "False = nothing reacted, the default"),
        Field("crosslink_distance_nm", crosslink.distance if crosslink else "", "nm",
              "Reaction distance between two lysine beads"),
        Field("crosslink_valence", crosslink.valence if crosslink else "", "bonds",
              "Max bonds per lysine — 1 = bifunctional crosslinker, 2 = trifunctional"),
        Field("crosslink_prob", crosslink.prob if crosslink else "", "",
              "P(react | within the cutoff at a check). 1.0 = diffusion-limited"),
        Field("crosslink_check_every", crosslink.check_every if crosslink else "", "steps",
              "MD steps between reaction checks"),
        Field("crosslink_start_step", crosslink.start_step if crosslink else "", "steps",
              "Reactions before this step were ignored (initial-placement contacts)"),
        Field("crosslink_k", crosslink.k if crosslink else "", "kJ/mol/nm^2",
              "Stiffness of a formed crosslink"),
        Field("crosslink_r0_nm", crosslink.r0 if crosslink else "", "nm",
              "Equilibrium length of a formed crosslink"),
        Field("crosslinks_formed", "" if n_crosslinks is None else n_crosslinks, "bonds",
              f"Rows in {EVENTS_FILENAME} — see crosslink_summary.txt for the intra/inter split"),

        # ---- how long it ran ----
        Field("total_steps", steps, "steps", "Total MD steps requested"),
        Field("steps_per_frame", wfreq, "steps", "Save interval (config wfreq) — steps between saved frames"),
        Field("n_frames_planned", n_frames_planned, "frames", "total_steps / steps_per_frame"),
        Field("n_frames_saved", "" if n_frames_saved is None else n_frames_saved, "frames",
              "Frames actually in the DCD (header NSET)"),
        Field("timestep_ps", DT_PS, "ps", "Langevin timestep — fixed by CALVADOS, not configurable"),
        Field("time_per_frame_ps", _round(wfreq * DT_PS, 6), "ps", "Simulated time between saved frames"),
        Field("total_time_ns", _round(steps * DT_PS / 1000, 6), "ns", "Total simulated time"),

        # ---- physics ----
        Field("temperature_K", config.get("temp", ""), "K", "Thermostat temperature"),
        Field("ionic_strength_M", config.get("ionic", ""), "M", "Ionic strength (Debye screening)"),
        Field("pH", config.get("pH", ""), "", "pH used for residue charges"),
        Field("z_wall_nm", _z_wall_from_expr(config.get("ext_force_expr")), "nm",
              "Depth of the surface-attachment potential well"),
        Field("ext_force_expr", config.get("ext_force_expr", ""), "", "External (surface) potential"),
        Field("cutoff_lj_nm", config.get("cutoff_lj", ""), "nm", "Ashbaugh-Hatch cutoff"),
        Field("cutoff_yu_nm", config.get("cutoff_yu", ""), "nm", "Yukawa (electrostatics) cutoff"),
        Field("friction_coeff_per_ps", config.get("friction_coeff", ""), "1/ps", "Langevin friction"),
        Field("force_field", Path(str(defaults.get("fresidues", ""))).name, "",
              "Residue parameter table used"),
        Field("topol", config.get("topol", ""), "", "CALVADOS initial-topology mode"),

        # ---- environment ----
        Field("platform", config.get("platform", ""), "", "OpenMM platform requested (CUDA/CPU)"),
        Field("random_number_seed", config.get("random_number_seed", "") or "", "",
              "Integrator seed (empty = OpenMM picked one, so the run is not reproducible bit-for-bit)"),
        Field("calvados_version", _package_version("calvados"), "", "calvados package version"),
        Field("openmm_version", _package_version("openmm"), "", "openmm package version"),
    ]


def write_lattice(runtime_dir: Path, *, nx: int, ny: int, spacing: float, margin: float) -> Path:
    """Record the grafting lattice in runtime/lattice.yaml. Returns the path written."""
    runtime_dir = Path(runtime_dir)
    path = runtime_dir / LATTICE_FILENAME
    path.write_text(
        "# Grafting lattice as prepare.py laid it out (nm).\n"
        f"nx: {nx}\n"
        f"ny: {ny}\n"
        f"spacing_nm: {spacing}\n"
        f"margin_nm: {margin}\n"
    )
    return path


def write_metadata(runtime_dir: Path) -> Path:
    """(Re)write runtime/metadata.csv, and copy it into results/. Returns the runtime path.

    The copy is what makes a committed result interpretable: a plot in
    ``results/`` is worth little without the spacing, box, temperature and
    crosslinking settings that produced it sitting next to it.
    """
    runtime_dir = Path(runtime_dir)
    path = runtime_dir / METADATA_FILENAME
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["key", "value", "unit", "description"])
        for field in collect_metadata(runtime_dir):
            writer.writerow([field.key, field.value, field.unit, field.description])
    publish(runtime_dir, METADATA_FILENAME)
    return path


def read_metadata(runtime_dir: Path) -> dict[str, str]:
    """{key: value} from runtime/metadata.csv, or {} if it hasn't been written yet.

    Values come back as strings exactly as stored; callers that want numbers
    convert what they need (see `as_float`).
    """
    path = Path(runtime_dir) / METADATA_FILENAME
    if not path.is_file():
        return {}
    with open(path, newline="") as f:
        return {row["key"]: row["value"] for row in csv.DictReader(f)}


def as_float(meta: dict[str, str], key: str, default: float = math.nan) -> float:
    """Numeric lookup in a read_metadata() dict, tolerant of missing/blank values."""
    try:
        return float(meta[key])
    except (KeyError, TypeError, ValueError):
        return default
