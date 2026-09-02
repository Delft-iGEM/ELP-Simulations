"""Where simulation output lives.

Everything a run produces — trajectories, checkpoints, the CG PDB — is written
into ``simulations/<sim_name>/runtime/``, and that adds up fast (a few hundred
MB per simulation). On DelftBlue ``/home`` is capped at 30 GB while ``/scratch``
gives you 5 TB, so the data belongs on ``/scratch``.

It can't simply *be* an absolute path there, though: ``template/run.py`` does
``from ..prepare import build_sim`` and is launched as
``python -m simulations.<sim_name>.runtime.run``, so ``runtime/`` has to stay
inside the ``simulations.<sim_name>`` package tree. So the folder stays where it
is and becomes a *symlink* into the data root — Python's importer, CALVADOS,
mdtraj and MDAnalysis all follow it, so nothing else in the project has to know.

The data root is, in order:

1. ``$ELP_DATA_DIR`` if set,
2. ``/scratch/$USER/elp-data/oefeningen`` if ``/scratch/$USER`` exists (DelftBlue),
3. nothing — on a laptop there's no ``/scratch``, so ``runtime/`` stays a plain
   folder inside the repo and everything works exactly as it did before.

``runtime/`` is gitignored, so the symlinks are never committed and each machine
resolves this for itself.
"""

from __future__ import annotations

import getpass
import os
from pathlib import Path

DEFAULT_SCRATCH_SUBDIR = "elp-data/oefeningen"


def data_root() -> Path | None:
    """Directory that holds every simulation's runtime data, or None if in-repo."""
    env = os.environ.get("ELP_DATA_DIR")
    if env:
        return Path(env).expanduser()

    scratch = Path("/scratch") / getpass.getuser()
    if scratch.is_dir():
        return scratch / DEFAULT_SCRATCH_SUBDIR

    return None


def runtime_target(sim_path: Path) -> Path | None:
    """Where ``sim_path``'s runtime data is really stored, or None if in-repo."""
    root = data_root()
    return None if root is None else root / sim_path.name


def ensure_runtime_dir(sim_path: Path) -> Path:
    """Return ``sim_path/runtime``, backed by the data root when there is one.

    Always returns the in-repo path so callers keep passing the same
    ``simulations/<sim_name>/runtime`` around; it just happens to be a symlink
    to ``/scratch`` on DelftBlue. Idempotent, like the rest of prepare.
    """
    link = sim_path / "runtime"
    target = runtime_target(sim_path)

    if target is None:
        link.mkdir(exist_ok=True)
        return link

    target.mkdir(parents=True, exist_ok=True)

    if link.is_symlink():
        if link.resolve() != target.resolve():
            link.unlink()
            link.symlink_to(target, target_is_directory=True)
    elif link.exists():
        # A real folder from before the data root was configured. Moving it here
        # would silently relocate an in-progress run's checkpoint, so leave it
        # and let the user migrate it deliberately with `sim migrate`.
        print(
            f"!  {link} is a real folder, not a link to {target} — "
            f"run 'sim migrate {sim_path.name}' to move it onto the data root."
        )
    else:
        link.symlink_to(target, target_is_directory=True)

    return link
