"""All 24 RGD arrangements grafted on a surface: 0.03 chains/nm^2, 500 ns, no TEA.

    .venv/bin/python tools/rgdarr_graft.py plan
    .venv/bin/python tools/rgdarr_graft.py init          # writes simulations/<name>/
    .venv/bin/python tools/rgdarr_graft.py submit --yes

24 runs, rgdg<NN>-<label>, one per sequence of the rgdarr sweep (the 352-aa
constructs in simulations/rgdarr<NN>-*/prepare.py, read verbatim):

  mode             preattached: 16 chains, each hangs from surface-bonded
                   lysines; 0.3 of ALL lysines bonded at t=0 (64 lysines ->
                   ~19, at least one per chain), fixed for the run
  concentration    0.03 chains/nm^2 (the rgdarr runs used 0.04)
  crosslinking     off (no lysine-lysine bonds form)
  TEA              off: stock CALVADOS2, 293.15 K, 0.19 M
  length           500 ns = 50,000,000 steps of 10 fs

Built exactly the way analyze.ipynb's submit cell builds a batch -
RunSpec -> tools.new_simulation.create_simulation -> sbatch ->
run_batch.record_runs - so the folders are indistinguishable from a CSV batch.

Walltime: rgdarr01-orig (the same 16 x 352 beads) ran 400 ns in 3650 s on an
A100, so 500 ns is ~76 min; 2:00:00 leaves ~1.6x.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

from tools.new_simulation import _project_root
from tools.run_batch import RunSpec, Shared, check_names, plan_runs

N = 24
SHARED = Shared(
    concentration=0.03,
    nmol=16,
    steps=50_000_000,
    walltime="2:00:00",
    crosslink_distance=None,
    mode="preattached",
    surface_preattached_fraction=0.3,
    tea=False,
)


def sources(root: Path) -> list[tuple[str, str]]:
    """(label, sequence) of rgdarr01..24, read from their prepare.py."""
    out = []
    for i in range(1, N + 1):
        (src,) = sorted((root / "simulations").glob(f"rgdarr{i:02d}-*/prepare.py"))
        seq = re.search(r'"([ACDEFGHIKLMNPQRSTVWY]{50,})"', src.read_text()).group(1)
        out.append((src.parent.name[len("rgdarr"):], seq))
    return out


def specs(root: Path) -> list[RunSpec]:
    s = SHARED
    out = []
    for row, (label, seq) in enumerate(sources(root), 1):
        if len(seq) != 352 or seq.count("RGD") != 4 or seq.count("K") != 4:
            raise ValueError(f"rgdarr{label}: unexpected sequence")
        out.append(RunSpec(
            name=f"rgdg{label}", sequence=seq, sequence_input=seq, sequence_label="",
            concentration=s.concentration, nmol=s.nmol, steps=s.steps, walltime=s.walltime,
            crosslink_distance=s.crosslink_distance, crosslink_prob=s.crosslink_prob,
            row=row, mode=s.mode, surface_preattached_fraction=s.surface_preattached_fraction,
            ionic=s.ionic, temperature=s.temperature, tea=s.tea, tea_gamma=s.tea_gamma))
    return out


def init(root: Path, todo: list[RunSpec]) -> dict[str, Path]:
    from tools.new_simulation import create_simulation
    from tools.run_batch import surface_kwargs

    s = SHARED
    made = {}
    for spec in todo:
        made[spec.name] = create_simulation(
            sequence=spec.sequence, concentration=spec.concentration,
            spacing_fraction=None, mass_concentration=None,
            length_measure=s.length_measure, margin_fraction=s.margin_fraction,
            box_height=s.box_height, crosslink_distance=spec.crosslink_distance,
            crosslink_prob=spec.crosslink_prob, ionic=spec.ionic,
            crosslink_valence=s.crosslink_valence, crosslink_check_every=s.crosslink_check_every,
            crosslink_k=s.crosslink_k, crosslink_r0=s.crosslink_r0,
            crosslink_ramp_steps=s.crosslink_ramp_steps,
            crosslink_start_step=s.crosslink_start_step,
            crosslink_selection=s.crosslink_selection,
            **surface_kwargs(spec, s),
            nmol=spec.nmol, steps=spec.steps, name=spec.name,
            temperature=spec.temperature, tea=spec.tea, tea_gamma=spec.tea_gamma,
            platform="CUDA", partition=s.partition, walltime=spec.walltime,
            cpu_per_task=s.cpu_per_task, max_n_save=s.max_n_save)
    return made


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "init", "submit"])
    parser.add_argument("--yes", action="store_true")
    args = parser.parse_args(argv)
    root = _project_root()
    todo = specs(root)

    if args.command == "plan":
        from tools.run_batch import print_batch_table
        print_batch_table(plan_runs(todo, SHARED))
        problems = check_names(todo)
        print("names free" if not problems else "\n".join(problems))
        return 0

    if args.command == "init":
        problems = check_names(todo)
        if problems:
            print("nothing prepared:\n  " + "\n  ".join(problems))
            return 1
        made = init(root, todo)
        print(f"prepared {len(made)} runs")
        return 0

    from tools.run_batch import parse_job_id, record_runs
    missing = [s.name for s in todo
               if not (root / "simulations" / s.name / "runtime" / "job.sh").is_file()]
    if missing:
        print(f"{len(missing)} run(s) not prepared - run `init` first: {missing[:4]}")
        return 1
    if not args.yes:
        print(f"would sbatch {len(todo)} jobs; re-run with --yes")
        return 0
    submitted = []
    try:
        for spec in todo:
            res = subprocess.run(["sbatch", str(root / "simulations" / spec.name / "runtime"
                                                / "job.sh")],
                                 cwd=root, capture_output=True, text=True)
            job = parse_job_id(res.stdout)
            print(f"  {spec.name:<26}{job or 'FAILED: ' + res.stderr.strip()}")
            if job:
                submitted.append((spec, job))
    finally:
        if submitted:
            record_runs(submitted, csv_name="rgdarr_graft.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
