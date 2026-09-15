"""
sim — ELP Simulations CLI

Commands:
  sim new                    Scaffold a new surface-attached simulation and prepare it
  sim prepare <simulation>   Generate the runtime/ folder (config, FASTA, job.sh, run.py)
  sim run     <simulation>   Run the simulation locally
  sim submit  <simulation>   Submit via sbatch on DelftBlue
  sim clean   <simulation>   Delete the runtime/ folder for a simulation
  sim list                   List all available simulations
  sim migrate <simulation>   Move an in-repo runtime/ folder onto the data root
                             (--all for every simulation)
  sim metadata <simulation>  Write runtime/metadata.csv (settings snapshot next to
                             the .dcd; --all for every simulation)
  sim distribution           Plot a residue's z-axis distribution across all frames
  sim crosslink <simulation> Detect crosslinks post-hoc in a finished trajectory
                             (the comparison against a reactive run)

Flags:
  sim run --clean <simulation>   Delete runtime/, prepare, then run

Autocomplete:
  sim --install-completion   Install shell completion (bash/zsh/fish/PowerShell)
  sim --show-completion      Print the completion script
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Annotated

import typer

from tools.paths import data_root, runtime_target
from tools.new_simulation import new as _new_simulation
from tools.z_distribution import distribution as _distribution

app = typer.Typer(
    name="sim",
    help="ELP Simulations CLI — prepare, run, and submit CALVADOS/OpenMM simulations.",
    no_args_is_help=True,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    """Return the project root (directory that contains pyproject.toml).

    We walk up from this file's location so the command works regardless of
    the current working directory.
    """
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    # Fallback: assume cwd is the project root
    return Path.cwd()


def _simulations_dir() -> Path:
    return _project_root() / "simulations"


def _available_simulations() -> list[str]:
    sims_dir = _simulations_dir()
    if not sims_dir.is_dir():
        return []
    return sorted(d.name for d in sims_dir.iterdir() if d.is_dir())


def _sim_completer() -> list[str]:
    """Typer shell-completion callback for simulation names."""
    return _available_simulations()


def _validate_sim(name: str) -> Path:
    """Check the simulation exists and return its path."""
    sim_path = _simulations_dir() / name
    if not sim_path.is_dir():
        available = _available_simulations()
        hint = "\n  ".join(available) if available else "(none found)"
        typer.echo(f"Error: simulation '{name}' not found in simulations/.", err=True)
        typer.echo(f"Available:\n  {hint}", err=True)
        raise typer.Exit(1)
    return sim_path


def _run(cmd: list[str], cwd: Path) -> None:
    """Run a command, streaming output.  Raises SystemExit on failure."""
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        raise typer.Exit(result.returncode)


def _clean_runtime(simulation: str, sim_path: Path) -> None:
    """Delete the runtime/ folder for a simulation if it exists.

    Where a data root is configured (see tools.paths) runtime/ is a symlink and
    the data itself lives on /scratch, so delete what the link points at before
    the link — shutil.rmtree refuses a symlink anyway, and dropping just the
    link would strand the trajectories on /scratch with nothing referring to
    them.
    """
    runtime_dir = sim_path / "runtime"
    if runtime_dir.is_symlink():
        target = runtime_dir.resolve()
        if target.is_dir():
            shutil.rmtree(target)
        runtime_dir.unlink()
        typer.echo(f"✓  Deleted runtime/ for {simulation} (data at {target})")
    elif runtime_dir.exists():
        shutil.rmtree(runtime_dir)
        typer.echo(f"✓  Deleted runtime/ for {simulation}")


def _migrate_runtime(simulation: str, sim_path: Path) -> bool:
    """Move a real runtime/ folder onto the data root, leaving a symlink behind.

    Returns True if anything moved. Prepare does this for new simulations by
    itself; this is for the ones that already have data sitting in the repo.
    """
    runtime_dir = sim_path / "runtime"
    target = runtime_target(sim_path)

    if target is None:
        typer.echo(
            "Error: no data root — set ELP_DATA_DIR, or run this on a machine "
            "with /scratch/$USER (DelftBlue).",
            err=True,
        )
        raise typer.Exit(1)

    if runtime_dir.is_symlink():
        typer.echo(f"·  {simulation}: already on the data root ({runtime_dir.resolve()})")
        return False
    if not runtime_dir.exists():
        typer.echo(f"·  {simulation}: no runtime/ folder — nothing to move")
        return False
    if target.exists():
        typer.echo(
            f"Error: {target} already exists — move or delete it first, "
            f"refusing to merge it with simulations/{simulation}/runtime/.",
            err=True,
        )
        raise typer.Exit(1)

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(runtime_dir), str(target))
    runtime_dir.symlink_to(target, target_is_directory=True)
    typer.echo(f"✓  {simulation}: moved to {target}")
    return True


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

SimName = Annotated[
    str,
    typer.Argument(
        help="Name of the simulation (folder under simulations/).",
        autocompletion=_sim_completer,
    ),
]


app.command(name="new")(_new_simulation)
app.command(name="distribution")(_distribution)


@app.command()
def prepare(
    simulation: SimName,
) -> None:
    """Generate the runtime/ folder for a simulation.

    Writes config.yaml, components.yaml, molecules.fasta, run.py, and job.sh
    into simulations/<simulation>/runtime/.
    """
    _validate_sim(simulation)
    typer.echo(f"▶  Preparing: {simulation}")
    root = _project_root()
    _run(
        [sys.executable, "-m", f"simulations.{simulation}.prepare"],
        cwd=root,
    )
    typer.echo(f"✓  Done — runtime/ folder ready at simulations/{simulation}/runtime/")


@app.command()
def run(
    simulation: SimName,
    clean: Annotated[bool, typer.Option("--clean", help="Delete runtime/, prepare, then run.")] = False,
) -> None:
    """Run a simulation locally.

    Requires that 'sim prepare <simulation>' has been run first.
    Use --clean to delete the runtime/ folder, re-prepare, and then run.
    """
    sim_path = _validate_sim(simulation)

    root = _project_root()

    if clean:
        _clean_runtime(simulation, sim_path)
        typer.echo(f"▶  Preparing: {simulation}")
        _run(
            [sys.executable, "-m", f"simulations.{simulation}.prepare"],
            cwd=root,
        )
        typer.echo(f"✓  Done — runtime/ folder ready at simulations/{simulation}/runtime/")

    run_script = sim_path / "runtime" / "run.py"
    if not run_script.exists():
        typer.echo(
            f"Error: runtime/run.py not found — run 'sim prepare {simulation}' first.",
            err=True,
        )
        raise typer.Exit(1)

    typer.echo(f"▶  Running: {simulation}")
    _run(
        [sys.executable, "-m", f"simulations.{simulation}.runtime.run"],
        cwd=root,
    )


@app.command()
def submit(
    simulation: SimName,
) -> None:
    """Submit a simulation to SLURM on DelftBlue (sbatch).

    Requires that 'sim prepare <simulation>' has been run first.
    """
    sim_path = _validate_sim(simulation)
    job_script = sim_path / "runtime" / "job.sh"
    if not job_script.exists():
        typer.echo(
            f"Error: runtime/job.sh not found — run 'sim prepare {simulation}' first.",
            err=True,
        )
        raise typer.Exit(1)

    typer.echo(f"▶  Submitting to SLURM: {simulation}")
    root = _project_root()
    _run(["sbatch", str(job_script)], cwd=root)


@app.command()
def clean(
    simulation: SimName,
) -> None:
    """Delete the runtime/ folder for a simulation.

    Has no effect if the runtime/ folder does not exist.
    """
    sim_path = _validate_sim(simulation)
    _clean_runtime(simulation, sim_path)


@app.command()
def migrate(
    simulation: Annotated[
        str | None,
        typer.Argument(help="Simulation to migrate.", autocompletion=_sim_completer),
    ] = None,
    all_sims: Annotated[bool, typer.Option("--all", help="Migrate every simulation.")] = False,
) -> None:
    """Move existing runtime/ data onto the data root and symlink it back.

    Simulations prepared from now on land there by themselves; this is a
    one-off for folders that were created before the data root was configured.
    """
    root = data_root()
    if root is None:
        typer.echo(
            "Error: no data root — set ELP_DATA_DIR, or run this on a machine "
            "with /scratch/$USER (DelftBlue).",
            err=True,
        )
        raise typer.Exit(1)

    if all_sims:
        targets = [(name, _simulations_dir() / name) for name in _available_simulations()]
    elif simulation is not None:
        targets = [(simulation, _validate_sim(simulation))]
    else:
        typer.echo("Error: give a simulation name, or --all.", err=True)
        raise typer.Exit(1)

    typer.echo(f"▶  Data root: {root}")
    moved = sum(_migrate_runtime(name, path) for name, path in targets)
    typer.echo(f"✓  Migrated {moved} simulation(s)")


@app.command()
def metadata(
    simulation: Annotated[
        str | None,
        typer.Argument(help="Simulation to (re)write metadata for.", autocompletion=_sim_completer),
    ] = None,
    all_sims: Annotated[bool, typer.Option("--all", help="Every prepared simulation.")] = False,
    show: Annotated[bool, typer.Option("--show", help="Print the metadata as well as writing it.")] = False,
) -> None:
    """(Re)write runtime/metadata.csv — the run's settings, beside its .dcd.

    New runs get this automatically (prepare writes it, run.py refreshes it when
    the run finishes). Use this to back-fill simulations that predate it, or to
    refresh the status/frame count of a run that has moved on since.
    """
    from tools.metadata import METADATA_FILENAME, read_metadata, write_metadata

    if all_sims:
        targets = [(name, _simulations_dir() / name) for name in _available_simulations()]
    elif simulation is not None:
        targets = [(simulation, _validate_sim(simulation))]
    else:
        typer.echo("Error: give a simulation name, or --all.", err=True)
        raise typer.Exit(1)

    written = 0
    for name, sim_path in targets:
        runtime_dir = sim_path / "runtime"
        if not (runtime_dir / "config.yaml").is_file():
            if not all_sims:
                typer.echo(
                    f"Error: {name} has no runtime/config.yaml — run 'sim prepare {name}' first.",
                    err=True,
                )
                raise typer.Exit(1)
            typer.echo(f"–  {name}: not prepared, skipped")
            continue

        write_metadata(runtime_dir)
        meta = read_metadata(runtime_dir)
        written += 1
        typer.echo(
            f"✓  {name}: {METADATA_FILENAME} written "
            f"({meta.get('status', '?')}, {meta.get('n_frames_saved') or 0}/"
            f"{meta.get('n_frames_planned', '?')} frames, "
            f"{meta.get('n_residues', '?')} residues x {meta.get('nmol', '?')} chains, "
            f"spacing {meta.get('spacing_nm', '?')} nm)"
        )
        if show:
            for key, value in meta.items():
                typer.echo(f"     {key:<30} {value}")

    if all_sims:
        typer.echo(f"✓  {written} simulation(s) documented")


@app.command()
def crosslink(
    simulation: SimName,
    distance: Annotated[
        float,
        typer.Option(help="Reaction distance between two lysine beads, IN NANOMETRES."),
    ],
    valence: Annotated[int, typer.Option(help="Max bonds per lysine.")] = 1,
    start_step: Annotated[
        int,
        typer.Option(help="Ignore contacts before this MD step (initial-placement artefacts)."),
    ] = 0,
    prob: Annotated[
        float,
        typer.Option(help="P(react | within the cutoff at a frame)."),
    ] = 1.0,
    seed: Annotated[int, typer.Option(help="Seed for the reaction RNG (only used if --prob < 1).")] = 0,
) -> None:
    """Detect crosslinks post-hoc in a finished trajectory — the first-passage rule.

    Walks the saved frames and bonds any lysine pair that has come within
    `--distance`, exactly as the reactive mode does, except that here the bond
    changes nothing: the chains carry on as if unbonded.

    This is the comparison a reactive run is for. A post-hoc pass over a
    *non-reactive* trajectory should count more crosslinks than the reactive run
    of the same system, because it counts encounters a real network would have
    prevented by tethering the chains at the first bond. A post-hoc conversion
    near 100% mostly says the cutoff was generous and nothing was ever restrained.

    Writes crosslink_events_posthoc.csv / crosslink_summary_posthoc.txt beside the
    trajectory, leaving any reactive run's own output untouched.
    """
    from tools.crosslink import CrosslinkSettings, post_hoc_events
    from tools.metadata import as_float, read_metadata

    sim_path = _validate_sim(simulation)
    runtime_dir = sim_path / "runtime"
    meta = read_metadata(runtime_dir)
    sim_name = meta.get("sim_name") or simulation
    traj = runtime_dir / f"{sim_name}.dcd"
    top = runtime_dir / "top.pdb"
    for path in (traj, top):
        if not path.is_file():
            typer.echo(f"Error: {path} not found — has {simulation} been run?", err=True)
            raise typer.Exit(1)

    n_save = int(as_float(meta, "steps_per_frame", 0) or 0)
    if not n_save:
        typer.echo("Error: steps_per_frame missing from metadata.csv — run "
                   f"'sim metadata {simulation}' first.", err=True)
        raise typer.Exit(1)

    settings = CrosslinkSettings(distance=distance, valence=valence, prob=prob,
                                 check_every=n_save, start_step=start_step, seed=seed)
    typer.echo(f"▶  Post-hoc crosslink detection on {simulation}: {settings.banner()}")
    xl = post_hoc_events(traj, top, settings, n_save=n_save, seed=seed)

    events = runtime_dir / "crosslink_events_posthoc.csv"
    summary = runtime_dir / "crosslink_summary_posthoc.txt"
    rows = xl.event_rows(n_save)
    import csv as _csv

    from tools.crosslink import CSV_COLUMNS

    with open(events, "w", newline="") as f:
        writer = _csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    steps = int(as_float(meta, "total_steps", 0) or 0)
    summary.write_text("\n".join(
        ["POST-HOC detection on a finished trajectory — the bonds changed nothing.",
         "Compare against a reactive run of the same system; see `sim crosslink --help`.",
         ""] + xl.summary_lines(n_save, steps)) + "\n")

    intra = sum(1 for r in rows if r["kind"] == "intra")
    typer.echo(f"✓  {len(rows)} crosslinks ({intra} intra, {len(rows) - intra} inter), "
               f"conversion {2 * len(rows) / (len(xl.sites) * valence):.3f}")
    typer.echo(f"   {events.name} and {summary.name} written")


@app.command(name="list")
def list_simulations() -> None:
    """List all available simulations."""
    sims = _available_simulations()
    if not sims:
        typer.echo("No simulations found in simulations/.")
        return
    typer.echo("Available simulations:")
    for s in sims:
        typer.echo(f"  {s}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app()
