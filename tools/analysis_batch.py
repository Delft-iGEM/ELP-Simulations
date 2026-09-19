"""Pointing the Analyze section at a whole CSV of runs instead of one name.

The Analyze section of ``analyze.ipynb`` is written around a single
``sim_name``. A batch, though, is planned as a CSV (``block-runs.csv``,
``used-runs.csv``, ``Example-csv.csv``, ...), and after it comes back from
DelftBlue you want the same cells run over every row of that file — without
babysitting it. Two things make that awkward, and this module handles both:

*Names.* A CSV row's ``Name`` is not always the folder name. ``create_simulation``
slugifies it (``My_ELP 2`` -> ``my-elp-2``), a row may have been renamed or
never submitted, and a run that is still queued has a folder but no trajectory
yet. Any of those would stop a notebook dead halfway through a batch. So
`resolve_targets` checks every name *up front*, prints one line per run it
cannot use — with the closest existing folder name, which is usually the typo —
and hands back only the runs that are actually there. Nothing is skipped
silently and one bad row never costs you the other nineteen.

*Repetition.* `run_for_each` re-executes the notebook's own analysis cells —
the ones tagged ``analysis`` in their cell metadata — once per simulation, in
order, in the notebook's namespace. There is no second copy of the analysis to
keep in sync: it is literally the cells you just ran interactively, with
``sim_name`` rebound each time. Each simulation's results land in its own
``simulations/<name>/analysis/`` as usual (see `tools.analysis_io`), and a
simulation whose cells raise is reported and abandoned — the batch moves on to
the next one rather than analyzing the next name with the previous run's
trajectory still loaded.
"""

from __future__ import annotations

import csv
import traceback
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Iterable

from tools.new_simulation import _slugify
from tools.run_batch import _norm_header, _sniff_delimiter

SIMULATIONS = Path("simulations")
DEFAULT_NOTEBOOK = "analyze.ipynb"
ANALYSIS_TAG = "analysis"


class SkipRest(Exception):
    """"There is nothing here to analyze" — not a failure.

    Raised by an analysis cell when this particular run has nothing for the
    rest of the section to work on (a sequence with no lysines has no
    crosslink cells to run, for example). In a batch it ends that simulation
    with a one-line note instead of a traceback, because nothing is wrong: the
    run simply isn't that kind of run.
    """

# Everything derived from a loaded trajectory. Dropped between simulations so a
# cell that fails can't leave the previous run's data behind for the next name
# to be analyzed against (silently wrong numbers are worse than a gap).
_CARRIED_OVER = ("u", "chains", "n_chains", "box", "mode", "surface_rows")


def available_simulations() -> list[str]:
    """Every folder under ``simulations/``, sorted."""
    if not SIMULATIONS.is_dir():
        return []
    return sorted(p.name for p in SIMULATIONS.iterdir() if p.is_dir())


def names_from_csv(csv_path: str | Path) -> list[str]:
    """The Name column of a runs CSV, in file order.

    Deliberately more forgiving than `tools.run_batch.read_runs`: analysis only
    needs the names, so a file without a Sequence column (or with columns this
    project has never heard of) is still usable here. Falls back to the first
    column when there is no recognisable name header.
    """
    path = Path(csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")

    lines = [line for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if not lines:
        return []

    rows = list(csv.reader(lines, delimiter=_sniff_delimiter(lines[0])))
    header = [_norm_header(column) for column in rows[0]]
    try:
        name_col = header.index("name")
    except ValueError:
        name_col = 0  # headerless or oddly named: assume the name comes first

    names = []
    for row in rows[1:]:
        if not row or not any(cell.strip() for cell in row):
            continue
        if row[0].lstrip().startswith("#"):     # a commented-out run
            continue
        if name_col < len(row) and row[name_col].strip():
            names.append(row[name_col].strip())
    return names


def expand_target(target: Any) -> list[str]:
    """Whatever the notebook was pointed at -> a list of raw run names.

    Accepts a single name, a CSV path, a shell-style pattern (``triblock-*``),
    or a list mixing all three; duplicates are dropped, order is kept.
    """
    if target is None:
        return []
    if isinstance(target, (str, Path)):
        items: Iterable[Any] = [target]
    else:
        items = target

    names: list[str] = []
    for item in items:
        text = str(item).strip()
        if not text:
            continue
        if text.lower().endswith(".csv") or Path(text).is_file():
            names.extend(names_from_csv(text))
        elif any(ch in text for ch in "*?["):
            matches = sorted(p.name for p in SIMULATIONS.glob(text) if p.is_dir())
            if matches:
                names.extend(matches)
            else:
                names.append(text)     # keep it, so it is reported as not found
        else:
            names.append(text)

    seen, unique = set(), []
    for name in names:
        if name not in seen:
            seen.add(name)
            unique.append(name)
    return unique


def trajectory_paths(sim_name: str) -> tuple[Path, Path]:
    """``(top.pdb, .dcd)`` for a simulation, raising a readable error if either is missing.

    The trajectory is normally ``runtime/<sim_name>.dcd``, but a folder that was
    renamed after the run keeps the old file name inside it — so fall back to
    whatever single ``.dcd`` is there rather than insisting on the folder name.
    """
    runtime = SIMULATIONS / sim_name / "runtime"
    top = runtime / "top.pdb"
    if not top.is_file():
        raise FileNotFoundError(f"{top} is missing — has {sim_name} been prepared and run?")

    named = runtime / f"{sim_name}.dcd"
    if named.is_file():
        return top, named
    dcds = sorted(runtime.glob("*.dcd"))
    if not dcds:
        raise FileNotFoundError(f"no .dcd trajectory in {runtime} — has {sim_name} produced frames yet?")
    if len(dcds) > 1:
        dcds.sort(key=lambda p: p.stat().st_size, reverse=True)
        print(f"!  {sim_name}: several .dcd files in runtime/ "
              f"({', '.join(p.name for p in dcds)}) — using the largest, {dcds[0].name}")
    else:
        print(f"!  {sim_name}: trajectory is {dcds[0].name}, not {sim_name}.dcd (folder renamed?)")
    return top, dcds[0]


def check_target(raw_name: str, existing: list[str] | None = None) -> tuple[str, str | None]:
    """``(folder name, problem)`` for one CSV name — ``problem`` is None when usable."""
    existing = available_simulations() if existing is None else existing
    try:
        slug = _slugify(raw_name)
    except Exception:
        return raw_name, "not a usable folder name"

    sim_dir = SIMULATIONS / slug
    if not sim_dir.is_dir():
        close = get_close_matches(slug, existing, n=1, cutoff=0.6)
        hint = f" — closest folder that does exist: {close[0]}" if close else ""
        return slug, f"no simulations/{slug}/ folder{hint}"

    runtime = sim_dir / "runtime"
    if not runtime.is_dir():
        # A dangling symlink means the data root differs from the machine the run
        # was on (see tools.paths) — worth saying so instead of "not run yet".
        if (sim_dir / "runtime").is_symlink():
            return slug, f"runtime/ points at {Path(sim_dir / 'runtime').readlink()}, which isn't there (wrong ELP_DATA_DIR?)"
        return slug, "prepared but never run (no runtime/)"

    if not (runtime / "top.pdb").is_file():
        return slug, "no runtime/top.pdb — prepared but not run (or run on another machine's data root)"
    if not any(runtime.glob("*.dcd")):
        return slug, "no .dcd in runtime/ yet — still queued or running?"
    return slug, None


def resolve_targets(target: Any, quiet: bool = False) -> list[str]:
    """Names to analyze, with every unusable one reported and dropped.

    This is the whole point of pointing the notebook at a CSV: a row whose
    folder is missing, misnamed, or still waiting on its trajectory costs you
    that row and nothing else.
    """
    raw_names = expand_target(target)
    existing = available_simulations()

    found: list[str] = []
    skipped: list[tuple[str, str, str]] = []
    for raw in raw_names:
        slug, problem = check_target(raw, existing)
        if problem is None:
            found.append(slug)
        else:
            skipped.append((raw, slug, problem))

    if not quiet:
        source = target if isinstance(target, (str, Path)) else f"{len(raw_names)} name(s)"
        print(f"{len(found)} of {len(raw_names)} run(s) from {source} can be analyzed")
        for raw, slug, problem in skipped:
            shown = raw if raw == slug else f"{raw} -> {slug}"
            print(f"   !  skipping {shown}: {problem}")
        for name in found:
            print(f"      {name}")
        if not found and raw_names:
            print("   nothing to analyze — check the names against `ls simulations/`")
    return found


def analysis_cells(notebook: str | Path = DEFAULT_NOTEBOOK,
                   tag: str = ANALYSIS_TAG) -> list[tuple[int, str]]:
    """``(cell number, source)`` for every code cell tagged ``tag`` in the notebook.

    Read off disk, so the batch runs the *saved* version of the cells: save the
    notebook after editing an analysis cell, then run the batch.
    """
    import json

    path = Path(notebook)
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist — pass the notebook's path explicitly")
    cells = json.loads(path.read_text())["cells"]
    return [
        (i, "".join(cell["source"]))
        for i, cell in enumerate(cells)
        if cell["cell_type"] == "code" and tag in (cell.get("metadata", {}).get("tags") or [])
    ]


def run_for_each(sim_names: Iterable[str],
                 namespace: dict[str, Any],
                 notebook: str | Path = DEFAULT_NOTEBOOK,
                 tag: str = ANALYSIS_TAG,
                 close_figures: bool = True) -> dict[str, str]:
    """Run every ``analysis``-tagged cell of the notebook once per simulation.

    ``namespace`` is the notebook's own globals (pass ``globals()``), so the
    cells see and define exactly what they do when you run them by hand. A cell
    that raises ends that simulation — the rest of its cells would be reading
    the previous run's variables — and the batch continues with the next name.
    Returns ``{name: "ok"}`` / ``{name: "failed in cell N: ..."}``.
    """
    import matplotlib.pyplot as plt

    cells = analysis_cells(notebook, tag)
    if not cells:
        raise RuntimeError(
            f"no code cells tagged {tag!r} in {notebook} — tag the analysis cells "
            f"(cell metadata: {{\"tags\": [\"{tag}\"]}}) and save the notebook first"
        )

    names = list(sim_names)
    results: dict[str, str] = {}
    for n, name in enumerate(names, start=1):
        print(f"\n{'=' * 78}\n[{n}/{len(names)}] {name}\n{'=' * 78}")
        for key in _CARRIED_OVER:
            namespace.pop(key, None)
        namespace["sim_name"] = name
        results[name] = "ok"
        for index, source in cells:
            try:
                exec(compile(source, f"<{Path(notebook).name} cell {index}>", "exec"), namespace)
            except SkipRest as reason:
                print(f"\n   {name}: {reason} — nothing further to analyze for this run.")
                results[name] = f"ok (stopped at cell {index}: {reason})"
                break
            except Exception as exc:
                print(f"\n!  {name}: cell {index} raised {type(exc).__name__}: {exc}")
                print("   (rest of this run skipped; continuing with the next one)")
                traceback.print_exc()
                results[name] = f"failed in cell {index}: {type(exc).__name__}: {exc}"
                break
        if close_figures:
            # Twelve figures per run x twenty runs is enough to exhaust the
            # backend's figure limit and the kernel's memory; they are already
            # written to analysis/ by save_figure.
            plt.close("all")

    ok = [name for name, status in results.items() if status.startswith("ok")]
    print(f"\n{'=' * 78}\nanalyzed {len(ok)} of {len(names)} run(s)")
    for name, status in results.items():
        if not status.startswith("ok"):
            print(f"   !  {name}: {status}")
    return results
