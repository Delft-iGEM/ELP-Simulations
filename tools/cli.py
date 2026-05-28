"""
sim — ELP Simulations CLI

Commands:
  sim prepare <simulation>   Generate the runtime/ folder (config, FASTA, job.sh, run.py)
  sim run     <simulation>   Run the simulation locally
  sim submit  <simulation>   Submit via sbatch on DelftBlue
  sim clean   <simulation>   Delete the runtime/ folder for a simulation
  sim list                   List all available simulations

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
    """Delete the runtime/ folder for a simulation if it exists."""
    runtime_dir = sim_path / "runtime"
    if runtime_dir.exists():
        shutil.rmtree(runtime_dir)
        typer.echo(f"✓  Deleted runtime/ for {simulation}")


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
