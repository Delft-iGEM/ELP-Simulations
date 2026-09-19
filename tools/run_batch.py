"""Many runs from one CSV — read it, check it, preview it, remember it.

One job per simulation stopped making sense once a job is a *set* of sequences
to compare. This module is the batch counterpart of ``tools.new_simulation``:
the planning notebook reads a CSV where every row is one run, previews all of
them at once, and then submits the lot — instead of editing one cell per
sequence and hoping the name hasn't been used before.

The CSV (see ``Example-csv.csv``) has one row per run and holds only what
actually differs between runs::

    Name ;Sequence;concentration;nmolecules;steps;crosslink distance;crosslink prob;mode;walltime

``mode`` picks how the chains attach to the surface: ``brush`` (blank — residue
0 is the fixed anchor, the pipeline as it always was), ``free`` or
``preattached`` (chains hang from surface-bonded lysines; see tools/surface.py).

Everything else — the margin, the length measure, the crosslink bond
parameters, the surface-attachment settings, the SLURM partition — is the same
for every run in the batch and is set once, in the notebook's preparation cell,
as a `Shared` object. A column
that isn't in the CSV, or a cell left blank, falls back to that same object, so
``concentration`` can live in the CSV for one batch and in the cell for the
next without changing anything else.

Three helpers matter to the notebook:

* `read_runs` — CSV -> list of `RunSpec`, with every problem in the file
  reported at once rather than one exception per rerun;
* `plan_run` / `print_summary` / `plot_lattice` — the per-run stats and
  space-filling preview, so previewing N runs is a loop and not N copies of a
  cell;
* `check_names` / `record_runs` — names have to be unique *before* anything is
  queued, and what did get queued is appended to ``used-runs.csv`` in the
  project root so the next batch knows those names and sequences are taken.

Concentration here is always chains/nm^2. `tools.new_simulation` also supports
spacing-as-a-fraction-of-chain-length and grafted mass per area; a batch mixing
those modes would not be comparable row to row, so the CSV doesn't offer them.
"""

from __future__ import annotations

import csv
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NamedTuple

import typer

from tools.elibpy import build_sequence_with_features, humanize_seq
from tools.new_simulation import (
    MIN_BOX_L,
    BoxPlan,
    _project_root,
    _slugify,
    _validate_sequence,
    lattice_positions,
    plan_box,
    plan_steps,
)
from tools.surface import MODES, SurfaceSettings, check_capacity

# Every name and sequence ever submitted, in the project root next to the CSV
# you plan from. Committed (unlike runtime/), so a name burned on DelftBlue is
# also burned in a fresh checkout somewhere else.
REGISTRY_FILENAME = "used-runs.csv"

REGISTRY_COLUMNS = [
    "name",
    "sequence",
    "n_residues",
    "concentration_chains_per_nm2",
    "nmol",
    "steps",
    "crosslink_distance_nm",
    "crosslink_prob",
    "mode",
    "walltime",
    "submitted_utc",
    "job_id",
    "csv",
]

# CALVADOS's fixed Langevin timestep (calvados/sim.py), for steps -> ns.
DT_PS = 0.01


class Shared(NamedTuple):
    """Settings that are the same for every run in a batch.

    Set once in the notebook's preparation cell. The first block is used as the
    fallback for CSV columns that are missing or left blank; the rest never
    appear in the CSV at all.
    """

    # ---- fallbacks for CSV columns ----
    concentration: float | None = None      # chains/nm^2
    nmol: int | None = None
    steps: int | None = None
    walltime: str | None = None             # SLURM HH:MM:SS
    crosslink_distance: float | None = None  # nm; None = no crosslinking
    crosslink_prob: float = 1.0
    mode: str = "brush"                      # brush | free | preattached
    surface_preattached_fraction: float | None = None   # preattached: fraction bound at t=0

    # ---- never in the CSV: identical for the whole batch ----
    length_measure: str = "contour"
    margin_fraction: float = 0.0
    box_height: float | None = None         # nm; None = 0.3 x contour, >= 50 nm (plan_box_height)
    max_n_save: int = 7000
    crosslink_valence: int = 1
    crosslink_check_every: int = 1000
    crosslink_k: float = 2000.0
    crosslink_r0: float = 0.6
    crosslink_ramp_steps: int = 500
    crosslink_selection: list[int] | None = None
    crosslink_start_step: int = 0
    # Surface attachment through lysines (modes free / preattached only; ignored
    # for brush). Siblings of the crosslink_* settings — see tools/surface.py.
    surface_distance: float = 0.8           # nm — a lysine this close to the tether plane binds
    surface_prob: float = 1.0               # P(bind | within that distance at a check)
    surface_max_fraction: float = 0.5       # global cap, fraction of ALL lysines
    surface_dynamic: bool | None = None     # binding during the run: None = free yes, preattached no
    surface_tether: float = 0.6             # nm above the wall onset a bound lysine sits
    surface_k: float = 2000.0               # kJ/mol/nm^2
    surface_ramp_steps: int = 500
    surface_check_every: int = 1000
    surface_start_step: int = 0
    surface_attraction: float = 0.0         # kJ/mol; 0 = wall stays purely repulsive
    surface_attraction_width: float = 0.5   # nm
    surface_seed: int | None = None         # None = drawn and recorded per run

    # ---- SLURM ----
    partition: str = "gpu-a100"
    cpu_per_task: str = "18"


class RunSpec(NamedTuple):
    """One CSV row, parsed and defaulted but not yet planned."""

    name: str                       # slugified — the simulations/<name>/ folder
    sequence: str                   # validated single-letter AAs
    sequence_input: str             # the cell as typed (shorthand or raw)
    sequence_label: str             # human-readable form of a shorthand spec
    concentration: float
    nmol: int
    steps: int
    walltime: str
    crosslink_distance: float | None
    crosslink_prob: float
    row: int                        # line number in the CSV, for error messages
    mode: str = "brush"             # brush | free | preattached
    surface_preattached_fraction: float | None = None


class RunPlan(NamedTuple):
    """A `RunSpec` with its box, lattice and step plan worked out."""

    spec: RunSpec
    plan: BoxPlan
    n_save: int
    n_frames: int
    actual_steps: int
    shared: Shared

    @property
    def name(self) -> str:
        return self.spec.name

    @property
    def total_time_ns(self) -> float:
        return self.actual_steps * DT_PS / 1000


# ---------------------------------------------------------------------------
# Reading the CSV
# ---------------------------------------------------------------------------

# Headers are matched loosely — lowercased, with anything that isn't a letter
# or digit treated as a separator — so "Name ", "nmolecules", "crosslink
# distance" and "Crosslink_Distance" all land in the right place.
COLUMN_ALIASES: dict[str, str] = {
    "name": "name", "sim": "name", "sim name": "name", "run": "name", "run name": "name",
    "sequence": "sequence", "seq": "sequence", "components": "sequence", "elp": "sequence",
    "concentration": "concentration", "conc": "concentration",
    "chains per nm2": "concentration", "concentration chains per nm2": "concentration",
    "nmolecules": "nmol", "nmol": "nmol", "n molecules": "nmol", "molecules": "nmol",
    "chains": "nmol", "n chains": "nmol",
    "steps": "steps", "n steps": "steps", "total steps": "steps",
    "crosslink distance": "crosslink_distance", "crosslink dist": "crosslink_distance",
    "xl distance": "crosslink_distance", "crosslink distance nm": "crosslink_distance",
    "crosslink prob": "crosslink_prob", "crosslink probability": "crosslink_prob",
    "xl prob": "crosslink_prob", "crosslink p": "crosslink_prob",
    "walltime": "walltime", "wall time": "walltime", "time": "walltime",
    "runtime": "walltime", "job time": "walltime",
    "mode": "mode", "attachment": "mode", "attachment mode": "mode", "surface mode": "mode",
    "surface": "mode",
    "preattached fraction": "surface_preattached_fraction",
    "preattached": "surface_preattached_fraction",
    "surface preattached fraction": "surface_preattached_fraction",
    "surface fraction": "surface_preattached_fraction",
}

# Cells meaning "not set" — fall back to Shared, exactly as a blank cell does.
BLANK_VALUES = {"", "-", "none", "n/a", "na", "default", "off"}

# SLURM's own grammar (sbatch(1), --time): "minutes", "minutes:seconds",
# "hours:minutes:seconds", "days-hours", "days-hours:minutes" and
# "days-hours:minutes:seconds". Note that a bare "H:MM" is *not* hours and
# minutes to SLURM — two colon-separated fields are minutes:seconds — so
# "2:00" buys two minutes, not two hours. Parsed with exactly SLURM's meaning
# so the notebook's walltime totals are what the scheduler will enforce.
WALLTIME_RE = re.compile(
    r"^(?:(?P<days>\d+)-(?P<dh>\d+)(?::(?P<dm>[0-5]?\d))?(?::(?P<ds>[0-5]?\d))?"
    r"|(?P<a>\d+)(?::(?P<b>[0-5]?\d))?(?::(?P<c>[0-5]?\d))?)$"
)


def _norm_header(header: str) -> str:
    """'Crosslink_distance (nm) ' -> 'crosslink distance' for alias matching."""
    cleaned = re.sub(r"\(.*?\)", " ", header or "")
    cleaned = re.sub(r"[^0-9a-zA-Z]+", " ", cleaned).strip().lower()
    return cleaned


def _is_blank(value: str | None) -> bool:
    return value is None or value.strip().lower() in BLANK_VALUES


def _number(value: str, field: str, row: int) -> float:
    """Parse a numeric cell, tolerating 1e6, 7_070_000 and a decimal comma.

    A semicolon-separated CSV is what Excel writes in a locale where the
    decimal separator *is* the comma, so "0,02" in the concentration column is
    a likely spelling of 0.02 rather than a typo — and since the comma can't be
    a field separator here, reading it as a decimal point is unambiguous.
    """
    text = value.strip().replace(" ", "").replace("_", "")
    if text.count(",") == 1 and "." not in text:
        text = text.replace(",", ".")
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"row {row}: {field} = {value!r} is not a number") from None


def parse_sequence(text: str) -> tuple[str, str]:
    """A Sequence cell -> (validated sequence, human-readable label).

    Two spellings, told apart by whether the cell contains a digit or a
    bracket:

    * a raw single-letter sequence, e.g. ``VPGVGVPGVGVPGVG``;
    * the ELP shorthand `tools.elibpy` uses, e.g. ``V5 K1 (VPGFG) V4`` —
      ``<letter><count>`` is that many ``VPG<letter>G`` repeats and
      ``(FRAGMENT)`` is inserted literally. Tokens may be separated by spaces,
      commas or plus signs.
    """
    raw = (text or "").strip()
    if not raw:
        raise ValueError("empty sequence")

    if re.search(r"[0-9()]", raw):
        tokens = re.findall(r"\([^)]*\)|[A-Za-z]+\d*", raw)
        if not tokens:
            raise ValueError(f"{raw!r} has no recognizable ELP tokens")
        # A bare letter with no count ("K V3") is inserted by elibpy as ONE
        # literal residue, not as a VPGKG pentapeptide — almost never what a
        # shorthand cell means, and invisible in the residue count. Refuse it.
        bare = [t for t in tokens if re.fullmatch(r"[A-Za-z]", t)]
        if bare:
            raise ValueError(
                f"{raw!r}: token(s) {', '.join(repr(t) for t in bare)} have no repeat count "
                f"— write e.g. '{bare[0].upper()}1' for one VPG{bare[0].upper()}G repeat, or "
                f"'({bare[0].upper()})' to insert the single residue literally"
            )
        components = [(token, token) for token in tokens]
        sequence = build_sequence_with_features(components)["seq"]
        return _validate_sequence(sequence), humanize_seq(components)

    sequence = _validate_sequence(raw)
    return sequence, f"{len(sequence)} aa"


def parse_walltime(value: str) -> int:
    """SLURM walltime -> seconds, with SLURM's meaning for every form it accepts.

    MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM and D-HH:MM:SS (see WALLTIME_RE). A
    two-field "H:MM" is minutes:seconds to SLURM, so it is here too — write
    "2:00:00" for two hours.
    """
    text = (value or "").strip()
    match = WALLTIME_RE.match(text)
    if not match:
        raise ValueError(f"{value!r} is not a SLURM walltime (use HH:MM:SS)")
    g = match.groupdict()
    if g["days"] is not None:                          # D-HH[:MM[:SS]]
        days, hours = int(g["days"]), int(g["dh"])
        minutes, seconds = int(g["dm"] or 0), int(g["ds"] or 0)
    elif g["c"] is not None:                           # HH:MM:SS
        days, hours, minutes, seconds = 0, int(g["a"]), int(g["b"]), int(g["c"])
    elif g["b"] is not None:                           # MM:SS
        days, hours, minutes, seconds = 0, 0, int(g["a"]), int(g["b"])
    else:                                              # MM
        days, hours, minutes, seconds = 0, 0, int(g["a"]), 0
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def format_duration(seconds: float) -> str:
    """Seconds -> H:MM:SS, for printing sums of walltimes."""
    seconds = int(round(seconds))
    return f"{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"


def _sniff_delimiter(first_line: str) -> str:
    """';' for the Excel-in-Europe dialect, ',' otherwise."""
    return ";" if first_line.count(";") >= first_line.count(",") else ","


def read_runs(csv_path: str | Path, shared: Shared) -> list[RunSpec]:
    """Parse `csv_path` into one `RunSpec` per row, filling blanks from `shared`.

    Every problem in the file is collected and raised together: fixing a batch
    of twenty runs one exception at a time is twenty reruns of the cell.
    """
    path = Path(csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")

    text = path.read_text(encoding="utf-8-sig")
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"{path} is empty")

    reader = csv.reader(lines, delimiter=_sniff_delimiter(lines[0]))
    rows = list(reader)
    header = [_norm_header(column) for column in rows[0]]

    unknown = [
        raw for raw, norm in zip(rows[0], header)
        if norm and norm not in COLUMN_ALIASES
    ]
    fields = [COLUMN_ALIASES.get(norm) for norm in header]
    if "name" not in fields or "sequence" not in fields:
        raise ValueError(
            f"{path} must have a Name and a Sequence column "
            f"(found: {', '.join(c.strip() for c in rows[0] if c.strip())})"
        )

    problems: list[str] = []
    if unknown:
        problems.append(
            f"unrecognized column(s): {', '.join(repr(c.strip()) for c in unknown)} — "
            f"known columns are {', '.join(sorted(set(COLUMN_ALIASES.values())))}"
        )

    specs: list[RunSpec] = []
    for line_no, row in enumerate(rows[1:], start=2):
        if not any(cell.strip() for cell in row):
            continue
        if row[0].lstrip().startswith("#"):     # a commented-out run
            continue

        cells: dict[str, str] = {}
        for field, cell in zip(fields, row):
            if field:
                cells[field] = cell

        def value(field: str, fallback: Any) -> Any:
            """The cell, or `fallback` from Shared when it's missing/blank."""
            cell = cells.get(field)
            return fallback if _is_blank(cell) else cell.strip()

        raw_name = (cells.get("name") or "").strip()
        if not raw_name:
            problems.append(f"row {line_no}: no name")
            continue

        try:
            name = _slugify(raw_name)
            sequence, label = parse_sequence(cells.get("sequence", ""))

            concentration = value("concentration", shared.concentration)
            nmol = value("nmol", shared.nmol)
            steps = value("steps", shared.steps)
            walltime = value("walltime", shared.walltime)
            for field, got in (("concentration", concentration), ("nmolecules", nmol),
                               ("steps", steps), ("walltime", walltime)):
                if got is None:
                    raise ValueError(
                        f"row {line_no}: {field} is blank (or has no column) in the CSV and "
                        f"has no fallback in the preparation cell — set one of the two"
                    )

            concentration = (concentration if isinstance(concentration, (int, float))
                             else _number(concentration, "concentration", line_no))
            nmol = nmol if isinstance(nmol, int) else int(_number(nmol, "nmolecules", line_no))
            steps = steps if isinstance(steps, int) else int(_number(steps, "steps", line_no))
            walltime = str(walltime)
            parse_walltime(walltime)     # validated here so sbatch can't reject it later

            distance = value("crosslink_distance", shared.crosslink_distance)
            if distance is not None and not isinstance(distance, (int, float)):
                distance = _number(distance, "crosslink distance", line_no)
            prob = value("crosslink_prob", shared.crosslink_prob)
            if not isinstance(prob, (int, float)):
                prob = _number(prob, "crosslink prob", line_no)

            mode = str(value("mode", shared.mode)).strip().lower()
            if mode not in MODES:
                raise ValueError(f"row {line_no}: mode = {mode!r} is not one of "
                                 f"{', '.join(MODES)} (blank = brush)")
            pre = value("surface_preattached_fraction", shared.surface_preattached_fraction)
            if pre is not None and not isinstance(pre, (int, float)):
                pre = _number(pre, "preattached fraction", line_no)
            if mode == "preattached" and pre is None:
                raise ValueError(f"row {line_no}: mode preattached needs a preattached fraction "
                                 "(column, or surface_preattached_fraction in the preparation cell)")
            if mode != "brush":
                if "K" not in sequence:
                    raise ValueError(f"row {line_no}: mode {mode} attaches chains through "
                                     "lysines, but the sequence has no K")

            if concentration <= 0:
                raise ValueError(f"row {line_no}: concentration must be > 0")
            if nmol < 1:
                raise ValueError(f"row {line_no}: nmolecules must be at least 1")
            if steps < 1:
                raise ValueError(f"row {line_no}: steps must be at least 1")
            if not 0.0 <= prob <= 1.0:
                raise ValueError(f"row {line_no}: crosslink prob must be between 0 and 1")
        except (ValueError, typer.BadParameter) as error:
            # _validate_sequence and _slugify come from the CLI and raise click's
            # BadParameter, whose useful text is on .message, not str().
            message = str(getattr(error, "message", None) or error)
            problems.append(message if message.startswith("row ") else f"row {line_no}: {message}")
            continue

        specs.append(RunSpec(
            name=name,
            sequence=sequence,
            sequence_input=(cells.get("sequence") or "").strip(),
            sequence_label=label,
            concentration=float(concentration),
            nmol=int(nmol),
            steps=int(steps),
            walltime=walltime,
            crosslink_distance=None if distance is None else float(distance),
            crosslink_prob=float(prob),
            row=line_no,
            mode=mode,
            surface_preattached_fraction=None if pre is None else float(pre),
        ))

    if not specs and not problems:
        problems.append(f"{path} has a header but no runs")

    if problems:
        raise ValueError(f"{path} could not be read:\n  - " + "\n  - ".join(problems))

    return specs


# ---------------------------------------------------------------------------
# Planning and previewing
# ---------------------------------------------------------------------------

def surface_settings_for(spec: RunSpec, shared: Shared) -> SurfaceSettings | None:
    """The row's surface settings (None for brush), validated against its lysine count."""
    if spec.mode == "brush":
        return None
    from tools.new_simulation import Z_WALL

    settings = SurfaceSettings(
        mode=spec.mode,
        distance=shared.surface_distance,
        prob=shared.surface_prob,
        max_fraction=shared.surface_max_fraction,
        preattached_fraction=spec.surface_preattached_fraction,
        dynamic=shared.surface_dynamic,
        k=shared.surface_k,
        tether=shared.surface_tether,
        ramp_steps=shared.surface_ramp_steps,
        check_every=shared.surface_check_every,
        start_step=shared.surface_start_step,
        attraction=shared.surface_attraction,
        attraction_width=shared.surface_attraction_width,
        seed=shared.surface_seed,
        z_wall=Z_WALL,
    )
    check_capacity(settings, spec.nmol, spec.nmol * spec.sequence.count("K"))
    return settings


def surface_kwargs(spec: RunSpec, shared: Shared) -> dict[str, Any]:
    """The mode and surface_* arguments to hand create_simulation for this row."""
    return dict(
        mode=spec.mode,
        surface_distance=shared.surface_distance,
        surface_prob=shared.surface_prob,
        surface_max_fraction=shared.surface_max_fraction,
        surface_preattached_fraction=spec.surface_preattached_fraction,
        surface_dynamic=shared.surface_dynamic,
        surface_k=shared.surface_k,
        surface_tether=shared.surface_tether,
        surface_ramp_steps=shared.surface_ramp_steps,
        surface_check_every=shared.surface_check_every,
        surface_start_step=shared.surface_start_step,
        surface_attraction=shared.surface_attraction,
        surface_attraction_width=shared.surface_attraction_width,
        surface_seed=shared.surface_seed,
    )


def plan_run(spec: RunSpec, shared: Shared) -> RunPlan:
    """Box, lattice and step plan for one row — the same maths create_simulation does."""
    plan = plan_box(
        spec.sequence, spec.concentration, spec.nmol,
        spacing_fraction=None,
        length_measure=shared.length_measure,
        mass_concentration=None,
        margin_fraction=shared.margin_fraction,
        box_height=shared.box_height,
    )
    n_save, n_frames, actual_steps = plan_steps(spec.steps, max_n_save=shared.max_n_save)
    # Fails here, at planning, if the cap cannot give every chain its bond.
    try:
        surface_settings_for(spec, shared)
    except ValueError as error:
        raise ValueError(f"{spec.name} (CSV row {spec.row}): {error}") from None
    return RunPlan(spec, plan, n_save, n_frames, actual_steps, shared)


def plan_runs(specs: list[RunSpec], shared: Shared) -> list[RunPlan]:
    return [plan_run(spec, shared) for spec in specs]


def print_summary(run: RunPlan, anchor_sigma: float | None = None) -> None:
    """The per-run stats block: chain, box, lattice, density, steps, warnings.

    `anchor_sigma` is the CALVADOS bead diameter of the anchor residue; pass it
    to get the "chains are closer than a single bead is wide" check as well.
    """
    spec, plan, shared = run.spec, run.plan, run.shared

    print(f"simulation:      {spec.name}   (CSV row {spec.row})")
    if spec.sequence_label and spec.sequence_label != f"{len(spec.sequence)} aa":
        # The label collapses every VPGxG pentapeptide to its guest letter, so
        # "G1 A1 G1 A1" reads as "GAGA" — which is NOT a literal GAGA block.
        # Say so, and show the real start of the chain, so the two can't be
        # confused when a literal (GAGAGA) fragment was what was meant.
        print(f"built sequence:  {spec.sequence_label}")
        print(f"                 (one letter per VPGxG pentapeptide — 'GA' above is "
              f"VPGGGVPGAG, not a literal GA; only (…) fragments are literal)")
        head = spec.sequence[:60]
        print(f"sequence starts: {head}{'…' if len(spec.sequence) > 60 else ''}")
    print(f"sequence:        {len(spec.sequence)} residues, {plan.chain_length:.2f} nm "
          f"({shared.length_measure} length), {plan.chain_mass_da / 1000:.2f} kDa")
    print(f"box:             [{plan.l_box_x}, {plan.l_box_y}, {plan.l_box_z}] nm"
          + ("" if shared.box_height is None else "   (height set by hand)"))
    print(f"box height:      {plan.l_box_z} nm — the box is periodic in z too, so the brush "
          f"plus the 4 nm cutoff must stay below it; {plan.l_box_z / plan.chain_length:.2f} x "
          f"the contour (finished runs reached 0.16-0.19 x; see z_max_nm in metadata.csv)")
    print(f"lattice:         {plan.nx} x {plan.ny} grafting points, {plan.spacing} nm apart")
    if plan.margin:
        print(f"margin:          {plan.margin} nm of empty surface on every side "
              f"({plan.margin_fraction} x the spacing) — chains start that far from the")
        print(f"                 periodic boundary, and neighbours across it are "
              f"{plan.spacing + 2 * plan.margin:.3f} nm apart instead of {plan.spacing} nm, so the")
        print(f"                 rim is a little less crowded than the middle.")
    else:
        print(f"margin:          none — the lattice tiles the box exactly, so the periodic "
              f"images continue it seamlessly")
        print(f"                 (infinite surface, every chain equally crowded), but chains "
              f"reach the boundary and wrap.")
    print(f"spacing:         {plan.spacing} nm = {plan.achieved_spacing_fraction:.3f} x the "
          f"chain's {plan.chain_length:.2f} nm length")
    print(f"density:         {plan.achieved_concentration:.5f} chains/nm² on the lattice "
          f"(asked for {spec.concentration})"
          + (f", {plan.box_concentration:.5f} averaged over the box (the margin dilutes it)"
             if plan.margin else ""))
    print(f"surface loading: {plan.achieved_mass_concentration:.4f} µg/cm² "
          f"({plan.achieved_mass_concentration * 10:.3f} mg/m²)"
          + (f"; {plan.box_mass_concentration:.4f} µg/cm² over the whole box"
             if plan.margin else ""))

    if plan.clamped:
        requested_spacing = 1 / spec.concentration ** 0.5
        short_side_cells = min(plan.nx, plan.ny) + 2 * shared.margin_fraction
        print(f"⚠  the requested {spec.concentration} chains/nm² would make the short side of "
              f"the box {short_side_cells} x {requested_spacing:.3f} = "
              f"{short_side_cells * requested_spacing:.2f} nm (lattice + margin), under the "
              f"{MIN_BOX_L} nm floor CALVADOS's cutoffs need — chains were spread to "
              f"{plan.spacing} nm instead ({plan.achieved_concentration:.5f} chains/nm², "
              f"{plan.achieved_mass_concentration:.4f} µg/cm²). Raise nmolecules to keep the "
              f"concentration you asked for.")
    if plan.aspect_ratio > 2.0:
        print(f"⚠  {spec.nmol} molecules only tile as {plan.nx} x {plan.ny} — a very elongated "
              f"box. Pick an nmolecules with a squarer factorisation.")
    if plan.spacing_too_tight:
        print(f"⚠  spacing ({plan.spacing} nm) is tighter than a rough own-size estimate for a "
              f"{len(spec.sequence)}-residue chain ({plan.min_sane_spacing:.2f} nm) — chains "
              f"will start heavily overlapped and the simulation is likely to diverge to NaN. "
              f"Ask for a lower concentration and/or fewer chains.")
    if anchor_sigma is not None and plan.spacing < anchor_sigma:
        print(f"⚠  spacing ({plan.spacing} nm) is smaller than a single bead's own diameter "
              f"({anchor_sigma:.3f} nm) — even the anchor beads themselves will overlap.")

    surface = surface_settings_for(spec, shared)
    if surface is None:
        print("mode:            brush — residue 0 of every chain is the fixed surface anchor")
    else:
        n_lys_chain = spec.sequence.count("K")
        n_lys = n_lys_chain * spec.nmol
        print(f"mode:            {spec.mode} — chains hang from surface-bonded lysines "
              f"({n_lys_chain} per chain, {n_lys} in all); residue 0 is an ordinary residue")
        print(f"                 a lysine within {surface.distance} nm ({surface.distance_angstrom:.1f} Å) "
              f"of the tether plane (z = {surface.z_pin:.3f} nm, {surface.tether} nm above the wall) "
              f"binds, p={surface.prob};")
        print(f"                 cap {surface.max_fraction} of all lysines "
              f"({int(surface.max_fraction * n_lys)}), binding during the run "
              f"{'on' if surface.is_dynamic else 'off'}, k={surface.k} kJ/mol/nm² over "
              f"{surface.ramp_steps} steps, checked every {surface.check_every} steps")
        if spec.mode == "preattached":
            print(f"                 {surface.preattached_fraction} of all lysines "
                  f"(~{round(surface.preattached_fraction * n_lys)}) bonded at t=0, at least one "
                  f"per chain")
        if surface.is_dynamic and not surface.start_step:
            print("⚠  surface_start_step is 0, so lysines the initial construction leaves near the "
                  "surface bind at the first check. Set it past equilibration.")

    if spec.crosslink_distance is not None:
        print(f"crosslinking:    ON — lysine pairs within {spec.crosslink_distance} nm "
              f"({spec.crosslink_distance * 10:.1f} Å) bond during the run; valence "
              f"{shared.crosslink_valence}, p={spec.crosslink_prob},")
        print(f"                 checked every {shared.crosslink_check_every} steps "
              f"({shared.crosslink_check_every * DT_PS:.1f} ps), bond k={shared.crosslink_k} "
              f"kJ/mol/nm² at r0={shared.crosslink_r0} nm over {shared.crosslink_ramp_steps} steps")
        if shared.crosslink_start_step:
            print(f"                 reactions ignored before step "
                  f"{shared.crosslink_start_step:,} "
                  f"({shared.crosslink_start_step * DT_PS / 1000:.3f} ns)")
        else:
            print("⚠  crosslink_start_step is 0, so contacts left over from the initial "
                  "placement count as reactions. Set it past equilibration — the Rg/RMSD "
                  "check in the Analyze section says where that is.")
    else:
        print("crosslinking:    off — nothing reacts")

    rounded = "" if run.actual_steps == spec.steps else f"  (rounded from {spec.steps:,})"
    print(f"steps:           {run.actual_steps:,}{rounded} = {run.total_time_ns:.3f} ns, "
          f"saving every {run.n_save:,} steps ({run.n_save * DT_PS:.3f} ps)")
    print(f"                 -> {run.n_frames:,} saved frames (len(u.trajectory) in the "
          f"Analyze section)")
    print(f"walltime:        {spec.walltime} requested from SLURM "
          f"({parse_walltime(spec.walltime) / 3600:.2f} h)")


def plot_lattice(run: RunPlan, anchor_sigma: float, ax=None):
    """The to-scale space-filling preview of one run's grafting lattice.

    Same positions the generated prepare.py places molecules at: row-major over
    the nx x ny lattice, each chain at the centre of its own spacing x spacing
    cell, the whole lattice shifted in by the margin. The shaded rim is that
    margin — the empty surface a chain has to cross before it reaches the
    periodic boundary and starts wrapping to the far side of the box; the blue
    discs are a rough own-size footprint of a chain, so overlapping discs mean
    chains that start on top of each other.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    spec, plan = run.spec, run.plan
    centers = lattice_positions(spec.nmol, plan.spacing, plan.margin)

    if ax is None:
        _fig, ax = plt.subplots(figsize=(6, 6 * plan.l_box_y / plan.l_box_x))

    if plan.margin:
        ax.add_patch(plt.Rectangle((0, 0), plan.l_box_x, plan.l_box_y, facecolor="#dd8452",
                                   edgecolor="none", alpha=0.10))
        ax.add_patch(plt.Rectangle((plan.margin, plan.margin),
                                   plan.l_box_x - 2 * plan.margin,
                                   plan.l_box_y - 2 * plan.margin,
                                   facecolor="white", edgecolor="#dd8452", linestyle=":",
                                   linewidth=0.8))
    for x, y in centers:
        ax.add_patch(plt.Rectangle((x - plan.spacing / 2, y - plan.spacing / 2),
                                   plan.spacing, plan.spacing, fill=False,
                                   edgecolor="#bbbbbb", linestyle="--", linewidth=0.6))
        ax.add_patch(Circle((x, y), radius=plan.min_sane_spacing, facecolor="#4c72b0",
                            edgecolor="none", alpha=0.12))
    for x, y in centers:
        ax.add_patch(Circle((x, y), radius=anchor_sigma / 2, facecolor="#4c72b0",
                            edgecolor="none", alpha=0.8))

    ax.set_xlim(0, plan.l_box_x)
    ax.set_ylim(0, plan.l_box_y)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title(f"{spec.name}: {plan.nx}x{plan.ny} grafting lattice, spacing={plan.spacing} nm\n"
                 f"= {plan.achieved_spacing_fraction:.3f} x chain length "
                 f"({plan.chain_length:.1f} nm), "
                 f"{plan.achieved_mass_concentration:.3f} µg/cm², "
                 f"box {plan.l_box_x} x {plan.l_box_y} nm\n"
                 + (f"orange rim = {plan.margin} nm margin ({plan.margin_fraction} x spacing); "
                    if plan.margin else "no margin — tiles seamlessly; ")
                 + "blue shading = rough chain footprint", fontsize=10)
    return ax


def print_batch_table(runs: list[RunPlan]) -> None:
    """One line per run plus the walltime totals.

    The totals are the point: a single job has to ask for the *sum* of the
    walltimes (the runs go one after another on its one GPU), while separate
    jobs each ask for their own and queue independently — so the sum is the
    walltime to put on one big job, and the maximum is the longest any single
    job would need.
    """
    header = (f"{'name':<24} {'res':>5} {'nmol':>5} {'chains/nm²':>11} {'spacing':>8} "
              f"{'box (nm)':>21} {'steps':>12} {'ns':>9} {'frames':>7} {'XL':>5} "
              f"{'mode':>11} {'walltime':>9}")
    print(header)
    print("-" * len(header))
    total_seconds = 0
    for run in runs:
        spec, plan = run.spec, run.plan
        seconds = parse_walltime(spec.walltime)
        total_seconds += seconds
        crosslink = "off" if spec.crosslink_distance is None else f"{spec.crosslink_distance:g}"
        mode = spec.mode
        if spec.mode == "preattached" and spec.surface_preattached_fraction is not None:
            mode = f"pre {spec.surface_preattached_fraction:g}"
        print(f"{spec.name:<24} {len(spec.sequence):>5} {spec.nmol:>5} "
              f"{plan.achieved_concentration:>11.5f} {plan.spacing:>8.3f} "
              f"{f'{plan.l_box_x:g} x {plan.l_box_y:g} x {plan.l_box_z:g}':>21} {run.actual_steps:>12,} "
              f"{run.total_time_ns:>9.2f} {run.n_frames:>7,} {crosslink:>5} {mode:>11} "
              f"{spec.walltime:>9}")
    print("-" * len(header))
    longest = max(parse_walltime(r.spec.walltime) for r in runs)
    print(f"{len(runs)} runs — one job back-to-back needs {format_duration(total_seconds)} of "
          f"walltime; as separate jobs the longest is {format_duration(longest)}")
    print(f"(measured walltimes of finished runs are in each run's "
          f"runtime/metadata.csv: run_wall_seconds / ns_per_hour)")


# ---------------------------------------------------------------------------
# Name registry — what has already been used
# ---------------------------------------------------------------------------

def anchor_tagged(sequence: str, mode: str = "brush") -> str:
    """The sequence as CALVADOS actually simulates it: in brush mode residue 0 becomes "Z".

    prepare.py rewrites the first residue of every chain as the "Z" anchor tag
    (a zero-mass bead OpenMM pins to the surface), and that is the spelling
    metadata.csv reports. The registry stores the same form so that a sequence
    read back from a finished run and one typed into the CSV compare equal. In
    the free/preattached modes residue 0 is left alone.
    """
    if mode != "brush":
        return sequence
    return f"Z{sequence[1:]}" if sequence else sequence


def _same_sequence(a: str, b: str) -> bool:
    """Equal up to the anchor tag, so a brush run and a free run of one ELP match."""
    return a == b or (len(a) == len(b) and len(a) > 1 and a[1:] == b[1:]
                      and "Z" in (a[0], b[0]))


def _cell(value: Any) -> Any:
    """Blank instead of 'nan' for a number that was never recorded."""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return "" if value is None else value


def registry_path(root: Path | None = None) -> Path:
    """``used-runs.csv`` in the project root."""
    return (Path(root) if root else _project_root()) / REGISTRY_FILENAME


def read_registry(path: Path | None = None) -> list[dict[str, str]]:
    """Every previously submitted run, oldest first. Missing file = no runs yet."""
    path = Path(path) if path else registry_path()
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8-sig") as f:
        return [row for row in csv.DictReader(f) if (row.get("name") or "").strip()]


def check_names(
    specs: list[RunSpec],
    *,
    registry: Path | None = None,
    simulations_dir: Path | None = None,
) -> list[str]:
    """Every reason these runs can't be submitted as named, or [] if they can.

    Three ways a name can be taken, all of them checked *before* anything is
    queued, because a half-submitted batch is worse than a rejected one:
    duplicated inside this CSV, already a folder under simulations/, or already
    in used-runs.csv from an earlier submission. Reusing a name would make
    create_simulation refuse mid-loop (after earlier rows were already queued),
    and — since job.sh sets --dependency=singleton — a duplicate job name also
    silently serializes the new job behind the old one.
    """
    problems: list[str] = []

    seen: dict[str, RunSpec] = {}
    for spec in specs:
        if spec.name in seen:
            problems.append(
                f"'{spec.name}' appears twice in the CSV (rows {seen[spec.name].row} "
                f"and {spec.row}) — names become simulations/<name>/ and must be unique"
            )
        else:
            seen[spec.name] = spec

    sims = Path(simulations_dir) if simulations_dir else _project_root() / "simulations"
    for spec in specs:
        if (sims / spec.name).exists():
            problems.append(
                f"'{spec.name}' (row {spec.row}) already exists as simulations/{spec.name}/ — "
                f"rename the run, or delete that folder if it was a mistake"
            )

    used = {row["name"]: row for row in read_registry(registry)}
    for spec in specs:
        row = used.get(spec.name)
        if row:
            when = row.get("submitted_utc", "?")
            problems.append(
                f"'{spec.name}' (row {spec.row}) was already submitted on {when} — see "
                f"{REGISTRY_FILENAME}"
            )

    return problems


def sequence_notes(specs: list[RunSpec], registry: Path | None = None) -> list[str]:
    """Sequences reused under another name — worth knowing, not worth blocking on.

    Running the same sequence twice is perfectly legitimate (a different
    concentration, a longer run, a repeat with another seed), so this only
    reports; `check_names` is what actually stops a submission.
    """
    notes: list[str] = []
    previous: dict[str, str] = {}
    for row in read_registry(registry):
        previous.setdefault(row.get("sequence", ""), row.get("name", "?"))

    within: dict[str, RunSpec] = {}
    for spec in specs:
        earlier = next((name for seq, name in previous.items()
                        if _same_sequence(seq, anchor_tagged(spec.sequence, spec.mode))), None)
        if earlier and earlier != spec.name:
            notes.append(f"'{spec.name}' has the same sequence as '{earlier}', "
                         f"which was submitted earlier")
        if spec.sequence in within:
            notes.append(f"'{spec.name}' and '{within[spec.sequence].name}' are the same "
                         f"sequence in this CSV")
        else:
            within[spec.sequence] = spec
    return notes


def record_runs(
    runs: list[tuple[RunSpec, str]],
    *,
    csv_name: str = "",
    path: Path | None = None,
) -> Path:
    """Append submitted runs to ``used-runs.csv``. `runs` is [(spec, job_id), ...].

    Called after `sbatch`, never before: a name is burned when a job carrying it
    is actually queued. Appends rather than rewrites, so the file is a
    chronological log and an interrupted batch keeps whatever did go out.
    """
    path = Path(path) if path else registry_path()
    exists = path.is_file()
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")

    if exists:
        _upgrade_registry_columns(path)

    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=REGISTRY_COLUMNS)
        if not exists:
            writer.writeheader()
        for spec, job_id in runs:
            writer.writerow({
                "name": spec.name,
                "sequence": anchor_tagged(spec.sequence, spec.mode),
                "mode": spec.mode,
                "n_residues": len(spec.sequence),
                "concentration_chains_per_nm2": _cell(spec.concentration),
                "nmol": _cell(spec.nmol),
                "steps": _cell(spec.steps),
                "crosslink_distance_nm": _cell(spec.crosslink_distance),
                "crosslink_prob": _cell(spec.crosslink_prob),
                "walltime": spec.walltime,
                "submitted_utc": now,
                "job_id": job_id,
                "csv": csv_name,
            })
    return path


def _upgrade_registry_columns(path: Path) -> None:
    """Rewrite a registry whose header predates a column, keeping every row.

    Older rows get a blank in the new column (blank mode = brush, which is what
    they were). Done in place, once, the first time a newer writer appends.
    """
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        if set(REGISTRY_COLUMNS) <= set(fieldnames):
            return
        rows = list(reader)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=REGISTRY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in REGISTRY_COLUMNS})


def spec_from_simulation(sim_dir: Path) -> RunSpec | None:
    """Rebuild a `RunSpec` from a simulation that already exists on disk.

    Prefers the run's ``runtime/metadata.csv`` (the settings as actually
    simulated) and falls back to parsing the sequence out of ``prepare.py`` when
    the runtime folder isn't on this machine. None when neither yields a
    sequence — there is nothing to record about a folder that has no run in it.
    """
    from tools.metadata import read_metadata

    sim_dir = Path(sim_dir)
    meta = read_metadata(sim_dir / "runtime")
    sequence = meta.get("sequence", "")
    if not sequence:
        prepare = sim_dir / "prepare.py"
        match = re.search(r'"[^"]+"\s*:\s*"([A-Z]+)"', prepare.read_text()) \
            if prepare.is_file() else None
        sequence = match.group(1) if match else ""
    if not sequence:
        return None

    def number(key: str) -> float:
        try:
            return float(meta.get(key, "") or "nan")
        except ValueError:
            return math.nan

    distance = meta.get("crosslink_distance_nm", "")
    mode = (meta.get("mode") or "brush").strip() or "brush"
    pre = meta.get("surface_preattached_fraction", "")
    return RunSpec(
        name=sim_dir.name,
        sequence=sequence,
        sequence_input=sequence,
        sequence_label=f"{len(sequence)} aa",
        concentration=number("concentration_chains_per_nm2"),
        nmol=int(number("nmol")) if not math.isnan(number("nmol")) else 0,
        steps=int(number("total_steps")) if not math.isnan(number("total_steps")) else 0,
        walltime=meta.get("walltime_requested", ""),
        crosslink_distance=float(distance) if distance else None,
        crosslink_prob=number("crosslink_prob"),
        row=0,
        mode=mode,
        surface_preattached_fraction=float(pre) if pre else None,
    )


def record_simulation(sim_dir: Path, *, job_id: str = "", csv_name: str = "",
                      path: Path | None = None) -> bool:
    """Append an already-prepared simulation to the registry. True if it was added.

    For runs that didn't come from a CSV batch — `sim submit` on a single
    simulation — so that the list of used names and sequences stays complete
    whichever route a job took. A name already in the registry is left alone.
    """
    sim_dir = Path(sim_dir)
    path = Path(path) if path else registry_path()
    if any(row["name"] == sim_dir.name for row in read_registry(path)):
        return False
    spec = spec_from_simulation(sim_dir)
    if spec is None:
        return False
    record_runs([(spec, job_id)], csv_name=csv_name, path=path)
    return True


def backfill_registry(path: Path | None = None, root: Path | None = None) -> list[str]:
    """Add every simulation already on disk to ``used-runs.csv``. Returns the names added.

    The registry only knows what it was told, so a project that predates it
    would start empty and claim that every existing sequence is new. This walks
    simulations/ and records each run it can read. Runs already listed are left
    alone, so it is safe to re-run after adding simulations by hand.
    """
    root = Path(root) if root else _project_root()
    path = Path(path) if path else registry_path(root)
    known = {row["name"] for row in read_registry(path)}

    rows: list[tuple[RunSpec, str]] = []
    for sim_dir in sorted((root / "simulations").iterdir()):
        if not sim_dir.is_dir() or sim_dir.name in known or sim_dir.name.startswith("__"):
            continue
        spec = spec_from_simulation(sim_dir)
        if spec is not None:
            rows.append((spec, ""))

    if rows:
        record_runs(rows, csv_name="(backfilled from simulations/)", path=path)
    return [spec.name for spec, _ in rows]


def parse_job_id(sbatch_output: str) -> str:
    """The job id out of sbatch's "Submitted batch job 1234567", or ""."""
    match = re.search(r"(\d+)", sbatch_output or "")
    return match.group(1) if match else ""


def anchor_bead_diameter(residues_csv: Path | str | None = None) -> float:
    """CALVADOS sigma (nm) of the "Z" anchor bead — the overlap check's yardstick."""
    path = Path(residues_csv) if residues_csv else _project_root() / "residues_CALVADOS2.csv"
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["one"] == "Z":
                return float(row["sigmas"])
    return math.nan
