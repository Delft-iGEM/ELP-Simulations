"""Analysis output that outlives the run — kept next to prepare.py, so it reaches git.

Everything a run produces lands in ``simulations/<sim>/runtime/``, which is a
symlink onto ``/scratch`` and is gitignored (see ``tools.paths``): trajectories
run to hundreds of MB and have no business in git. ``sim clean`` deletes that
folder outright.

That leaves the *results* of an analysis — the numbers and plots you would put
in a report — with nowhere durable to live. They sit beside a ``.dcd`` on
scratch until someone cleans the run, and they never reach GitHub, so anyone
cloning the repo gets every simulation's recipe and none of its findings.

This module gives each simulation a ``results/`` folder *beside* ``prepare.py``:

    simulations/<sim>/
        prepare.py      the recipe     (tracked)
        results/        what came out  (tracked)   <- this module
        runtime/        the raw data   (gitignored symlink to /scratch)

``results/`` is small, tracked, and survives ``sim clean``. Text results are
written as ``.txt``, tabular ones as ``.csv`` and plots as ``.png``, all of
which GitHub renders in the browser without anyone cloning or downloading.

Only small, final artefacts belong here. The trajectory, the topology PDB and
the crosslink checkpoint stay in ``runtime/`` — they are raw data or restart
state rather than results, and ``MAX_PUBLISH_BYTES`` refuses anything large
enough to suggest a mistake, because a file committed by accident stays in the
history even after it is deleted.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from tools.paths import project_root

RESULTS_DIRNAME = "results"

# A result is a summary, a small table or a plot; nothing of that sort is
# megabytes. The cap is a guard against a .dcd being published by a typo, not a
# considered opinion about how big a legitimate plot may be.
MAX_PUBLISH_BYTES = 5 * 1024 * 1024


def sim_dir(sim: str | Path) -> Path:
    """The in-repo ``simulations/<sim>`` folder, from a name or any path to it."""
    if isinstance(sim, str) and "/" not in sim:
        return project_root() / "simulations" / sim
    path = Path(sim)
    return path if path.is_absolute() else project_root() / path


def sim_dir_for_runtime(runtime_dir: str | Path) -> Path | None:
    """The simulation folder owning ``runtime_dir``, or None if it can't be found.

    Two shapes have to work. On a laptop ``runtime/`` is a real folder inside
    the simulation, so the parent is the answer. On DelftBlue it is a symlink
    onto the data root — and ``template/run.py`` resolves it, so by the time a
    reactive run calls in here the path is ``/scratch/.../<sim>`` whose parent
    holds no ``prepare.py`` at all. The data root stores each run under its own
    simulation name (``tools.paths.runtime_target``), so the last component is
    the name to look up.

    Returns None rather than guessing if neither shape matches: publishing is a
    convenience, and a run must never die because it could not file a copy.
    """
    runtime_dir = Path(runtime_dir)

    parent = runtime_dir.parent
    if (parent / "prepare.py").is_file():
        return parent

    candidate = project_root() / "simulations" / runtime_dir.name
    if (candidate / "prepare.py").is_file():
        return candidate

    return None


def results_dir(sim: str | Path, create: bool = True) -> Path:
    """``simulations/<sim>/results``, created on demand.

    Refuses a simulation that doesn't exist rather than creating the folder for
    it. A typo'd or stale ``sim_name`` in the notebook would otherwise quietly
    scatter a stray ``simulations/<typo>/results/`` next to the real ones, and
    the first sign of it would be the mystery folder in a commit.
    """
    base = sim_dir(sim)
    if not base.is_dir():
        raise FileNotFoundError(
            f"No simulation folder at {base} — check the name. "
            f"Results are filed beside that simulation's prepare.py."
        )
    path = base / RESULTS_DIRNAME
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def save_text(sim: str | Path, name: str, text: str) -> Path:
    """Write a text result, making sure it ends with exactly one newline."""
    path = results_dir(sim) / name
    path.write_text(text.rstrip("\n") + "\n")
    return path


def save_figure(sim: str | Path, name: str, fig, dpi: int = 150) -> Path:
    """Save a matplotlib figure into the simulation's results folder.

    ``bbox_inches="tight"`` because these are read at thumbnail size on GitHub,
    where a wide band of whitespace costs more than it does on screen.
    """
    path = results_dir(sim) / name
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def publish(runtime_dir: str | Path, *names: str, quiet: bool = False) -> list[Path]:
    """Copy finished result files out of ``runtime/`` into the tracked ``results/``.

    Named files that don't exist are skipped without complaint — a non-reactive
    run has no crosslink summary to publish, and that is not an error. Anything
    over ``MAX_PUBLISH_BYTES`` is skipped loudly instead of silently, since the
    only way to get there is by naming the wrong file.

    Returns the paths actually written, so a caller can report them.
    """
    runtime_dir = Path(runtime_dir)
    target_sim = sim_dir_for_runtime(runtime_dir)
    if target_sim is None:
        if not quiet:
            print(f"!  can't tell which simulation {runtime_dir} belongs to — "
                  "results not published")
        return []

    written: list[Path] = []
    for name in names:
        source = runtime_dir / name
        if not source.is_file():
            continue
        if source.stat().st_size > MAX_PUBLISH_BYTES:
            print(f"!  {name} is larger than {MAX_PUBLISH_BYTES // (1024 * 1024)} MB — "
                  "left in runtime/ rather than committed")
            continue
        destination = results_dir(target_sim) / name
        shutil.copyfile(source, destination)
        written.append(destination)

    return written
