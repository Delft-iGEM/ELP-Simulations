"""Where analysis results live.

``analyze.ipynb`` reads a trajectory out of ``runtime/`` (on ``/scratch``, not
committed — see ``tools.paths``) but the *results* of looking at it — the
equilibration verdict, the crosslinking counts, the plots — are small (a few
KB of text, a few hundred KB of PNG) and are exactly what you'd want to carry
around without re-running the notebook: written into
``simulations/<sim_name>/analysis/``, right next to that simulation's
``prepare.py``, so they travel with the code. `git add`/`push` them from
wherever you ran the notebook and a DelftBlue checkout sees the same numbers
and figures without needing MDAnalysis, a trajectory, or a kernel.

``analysis/`` is deliberately *not* in ``.gitignore`` (unlike ``runtime/``).
"""

from __future__ import annotations

import io
import sys
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt


def analysis_dir(sim_name: str) -> Path:
    """``simulations/<sim_name>/analysis``, created on first use."""
    path = Path("simulations") / sim_name / "analysis"
    path.mkdir(parents=True, exist_ok=True)
    return path


class _Tee:
    """Writes to several streams at once, so capturing output doesn't hide it."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for stream in self._streams:
            stream.write(data)

    def flush(self):
        for stream in self._streams:
            stream.flush()


@contextmanager
def save_text(sim_name: str, name: str):
    """Mirror everything printed in the ``with`` block to ``analysis/<name>.txt``.

    Nothing about the printing itself has to change — this just tees stdout
    to a buffer for the duration of the block and writes it out on exit, so a
    cell keeps showing its normal output in the notebook and gets a copy on
    disk for free.
    """
    path = analysis_dir(sim_name) / f"{name}.txt"
    buffer = io.StringIO()
    real_stdout = sys.stdout
    sys.stdout = _Tee(real_stdout, buffer)
    try:
        yield path
    finally:
        sys.stdout = real_stdout
        path.write_text(buffer.getvalue())


def save_figure(sim_name: str, name: str, fig=None) -> Path:
    """Save ``fig`` (default: the current figure) to ``analysis/<name>.png``."""
    fig = fig if fig is not None else plt.gcf()
    path = analysis_dir(sim_name) / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    return path
