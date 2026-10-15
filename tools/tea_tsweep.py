"""Temperature sweep of the unchanged 357-aa RGD construct: 275-380 K in 5 K steps.

    .venv/bin/python tools/tea_tsweep.py plan
    .venv/bin/python tools/tea_tsweep.py init       # writes prepare.py and runs it
    .venv/bin/python tools/tea_tsweep.py submit --yes

22 runs, single free chains, 400 ns, CALVADOS-TEA gamma 6 - the conditions of
g6-orig-285k / -370k (tools/tea_g6.py). Every run uses the same seed, so
temperature is the only difference. 290-350 K were run first (2026-10-15);
275-285 and 355-380 K extend them. `init` and `submit` skip runs that are
already prepared / have a trajectory, so extending never touches those.

275 K is below TEA's 280-380 K fit range and is extrapolated (recorded in
runtime/tea.yaml); at gamma 6 Val's lambda there is -0.04, slightly
repulsive. Read that point as an extrapolation, not as data-backed.

Box 80 nm, gpu-v100, 8 CPUs x 4 GB, 0:40:00 - a g6 357-aa chain took
1000-1020 s per 400 ns on V100.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from tools.new_simulation import _project_root
from tools.tea import DEFAULT_GAMMA, T_MAX, T_MIN
from tools.tea_sequences import original
from tools.tea_temps import PREPARE

GAMMA = DEFAULT_GAMMA
TEMPERATURES = [275.0 + 5 * i for i in range(22)]          # 275 .. 380
NS = 400.0
N_SAVE = 7000
BOX_NM = 80.0
WALLTIME = "0:40:00"
PARTITION = "gpu-v100"
CPUS = "8"
MEM_PER_CPU = "4000"
SEED = 20261016


def run_name(t: float) -> str:
    return f"g6-tsw-{int(t)}k"


def cells(root: Path) -> list[dict]:
    seq = original(root)
    return [{"name": run_name(T), "sequence": seq, "temperature": T} for T in TEMPERATURES]


def init(root: Path, made: list[dict]) -> None:
    for c in made:
        if not T_MIN - 5 <= c["temperature"] <= T_MAX:
            raise ValueError(f"{c['temperature']} K is too far outside TEA's "
                             f"{T_MIN:g}-{T_MAX:g} K range")
    n_frames = int(round(NS * 1000.0 / (N_SAVE * 0.01)))
    for c in made:
        sim_dir = root / "simulations" / c["name"]
        sim_dir.mkdir(parents=True, exist_ok=True)
        text = PREPARE.format(
            title=f"Original 357-aa RGD construct, temperature sweep, {c['temperature']:g} K.",
            sequence=c["sequence"], temperature=c["temperature"], gamma=GAMMA, tea=True,
            extrapolate=c["temperature"] < T_MIN, box=BOX_NM, n_save=N_SAVE, n_frames=n_frames, seed=SEED,
            ns=NS, tea_line=f"CALVADOS-TEA on, gamma {GAMMA:g}.", partition=PARTITION,
            walltime=WALLTIME, cpus=CPUS)
        # The template's job.sh asks 8 GB per CPU, more than fits beside other
        # jobs on a V100 node.
        hook = f'.replace("{{cpu_per_task}}", "{CPUS}")'
        assert text.count(hook) == 1
        text = text.replace(hook, hook + f' \\\n        .replace("--mem-per-cpu=8000", '
                                        f'"--mem-per-cpu={MEM_PER_CPU}")')
        (sim_dir / "prepare.py").write_text(text)
        subprocess.run([str(root / ".venv/bin/python"), "-m",
                        f"simulations.{c['name']}.prepare"],
                       cwd=root, check=True, capture_output=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "init", "submit"])
    parser.add_argument("--yes", action="store_true")
    args = parser.parse_args(argv)
    root = _project_root()
    made = cells(root)

    if args.command == "plan":
        print(f"{len(made)} runs, {len(made[0]['sequence'])} aa, single free chains, "
              f"{NS:g} ns, {BOX_NM:g} nm, TEA gamma {GAMMA:g}, {PARTITION} {WALLTIME}")
        print("  " + ", ".join(c["name"] for c in made))
        return 0

    if args.command == "init":
        made = [c for c in made
                if not (root / "simulations" / c["name"] / "runtime" / "job.sh").is_file()]
        init(root, made)
        print(f"prepared {len(made)} runs")
        return 0

    from tools.run_batch import parse_job_id, record_simulation
    missing = [c["name"] for c in made
               if not (root / "simulations" / c["name"] / "runtime" / "job.sh").is_file()]
    if missing:
        print(f"{len(missing)} run(s) not prepared - run `init` first: {missing[:4]}")
        return 1
    if not args.yes:
        new = [c["name"] for c in made
               if not list((root / "simulations" / c["name"] / "runtime").glob("*.dcd"))]
        print(f"would sbatch {len(new)} jobs: {', '.join(new)}; re-run with --yes")
        return 0
    made = [c for c in made
            if not list((root / "simulations" / c["name"] / "runtime").glob("*.dcd"))]
    for c in made:
        sim_dir = root / "simulations" / c["name"]
        res = subprocess.run(["sbatch", str(sim_dir / "runtime" / "job.sh")],
                             cwd=root, capture_output=True, text=True)
        job = parse_job_id(res.stdout)
        print(f"  {c['name']:<16}{job or 'FAILED: ' + res.stderr.strip()}")
        if job:
            record_simulation(sim_dir, job_id=job, csv_name="tea_tsweep.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
