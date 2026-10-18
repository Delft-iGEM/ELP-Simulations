"""Submodule-4 main batch (Submodule_4_mains_new.txt) on gpu-v100.

    .venv/bin/python tools/sub4_mains.py plan
    .venv/bin/python tools/sub4_mains.py init          # writes simulations/<name>/
    .venv/bin/python tools/sub4_mains.py submit --yes

40 preattached runs, 16 chains, 400 ns, no crosslinking, TEA off, 293.15 K -
the settings of the earlier sub4 batches (batches/submodule4_extra_20.csv), but on the
V100 partition because the A100 queue is busy.

Walltime: the CSV's 1:30/2:00 was sized for A100, where sub4 runs took 56-84
min (270-745 residues). No 16-chain run has timed on V100 yet, so the walltime
is 2:30 below 100 pentads and 3:00 from 100 pentads up. Jobs ask for 8 CPUs x
4 GB, which fits beside other jobs on a V100 node (the template's 8 GB per CPU
does not).
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from tools.new_simulation import _project_root
from tools.run_batch import Shared, check_names, plan_runs, read_runs

CSV = "Submodule_4_mains_new.txt"
PARTITION = "gpu-v100"
CPUS = "8"
MEM_PER_CPU = "4000"
SHARED = Shared(mode="preattached", surface_preattached_fraction=0.3,
                partition=PARTITION, cpu_per_task=CPUS)


def walltime(spec) -> str:
    return "3:00:00" if len(spec.sequence) >= 500 else "2:30:00"


def specs(root: Path):
    return [s._replace(walltime=walltime(s)) for s in read_runs(root / CSV, SHARED)]


def init(root: Path, todo) -> None:
    from tools.new_simulation import create_simulation
    from tools.run_batch import surface_kwargs

    s = SHARED
    for spec in todo:
        sim_dir = create_simulation(
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
            platform="CUDA", partition=PARTITION, walltime=spec.walltime,
            cpu_per_task=CPUS, max_n_save=s.max_n_save)
        # 8 GB per CPU does not fit beside other jobs on a V100 node; patch both
        # the job script and prepare.py, which would rewrite it
        for f in (Path(sim_dir) / "runtime" / "job.sh", Path(sim_dir) / "prepare.py"):
            if f.is_file():
                f.write_text(f.read_text().replace("--mem-per-cpu=8000",
                                                   f"--mem-per-cpu={MEM_PER_CPU}"))
        job = (Path(sim_dir) / "runtime" / "job.sh").read_text()
        assert f"--mem-per-cpu={MEM_PER_CPU}" in job and f"--partition={PARTITION}" in job, sim_dir


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
        init(root, todo)
        print(f"prepared {len(todo)} runs")
        return 0

    from tools.run_batch import parse_job_id, record_runs
    missing = [s.name for s in todo
               if not (root / "simulations" / s.name / "runtime" / "job.sh").is_file()]
    if missing:
        print(f"{len(missing)} run(s) not prepared - run `init` first: {missing[:4]}")
        return 1
    if not args.yes:
        print(f"would sbatch {len(todo)} jobs to {PARTITION}; re-run with --yes")
        return 0
    submitted = []
    try:
        for spec in todo:
            res = subprocess.run(["sbatch", str(root / "simulations" / spec.name / "runtime"
                                                / "job.sh")],
                                 cwd=root, capture_output=True, text=True)
            job = parse_job_id(res.stdout)
            print(f"  {spec.name:<40}{job or 'FAILED: ' + res.stderr.strip()}")
            if job:
                submitted.append((spec, job))
    finally:
        if submitted:
            record_runs(submitted, csv_name=CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
