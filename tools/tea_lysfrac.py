"""(VPGIG)96 with the first 0, 10, ..., 100 % of pentads turned into VPGKG.

    .venv/bin/python tools/tea_lysfrac.py plan
    .venv/bin/python tools/tea_lysfrac.py init       # writes prepare.py and runs it
    .venv/bin/python tools/tea_lysfrac.py submit --yes

22 runs: 11 sequences x {285, 370} K, single free chains, 400 ns, CALVADOS-TEA
at gamma 6 - the conditions of tools/tea_g6.py.

Sequences
---------
96 pentads, 480 residues. Sequence k has its first round(k/10 * 96) pentads as
VPGKG and the rest VPGIG: 0, 10, 19, 29, 38, 48, 58, 67, 77, 86, 96 lysines.
The lysines are a block at the N terminus, so this varies the amount of lysine
and the length of one charged block together.

Box
---
200 nm, chosen so the chain cannot meet its own periodic image whatever it
does: the contour is 480 x 0.38 = 182 nm, and 182 nm plus the 4 nm Yukawa
cutoff on each side is 190 nm. (VPGKG)96 carries +96 and is the most swollen
chain this project has run, so a bound from measured extents of neutral ELPs
would not cover it. Box volume is nearly free - CALVADOS uses cutoff
electrostatics, not PME, so cost scales with beads.

Hardware
--------
gpu-v100, 8 CPUs x 4 GB. On the g6 set V100 averaged 19 min per 400 ns for a
357-residue chain against 22 min on A100, and the A100 queue was the slow one.
480 residues is ~1.35x the beads; 1:00:00 is about 2x the expected time.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from tools.new_simulation import _project_root
from tools.tea import DEFAULT_GAMMA, T_MAX, T_MIN
from tools.tea_temps import PREPARE

GAMMA = DEFAULT_GAMMA
TEMPERATURES = [285.0, 370.0]
N_PENTADS = 96
FRACTIONS = [k / 10 for k in range(11)]
NS = 400.0
N_SAVE = 7000
BOX_NM = 200.0
WALLTIME = "1:00:00"
PARTITION = "gpu-v100"
CPUS = "8"
MEM_PER_CPU = "4000"
SEED0 = 20261008


def sequence(fraction: float) -> str:
    n_lys = round(fraction * N_PENTADS)
    return "VPGKG" * n_lys + "VPGIG" * (N_PENTADS - n_lys)


def cells() -> list[dict]:
    out = []
    for i, f in enumerate(FRACTIONS):
        seq = sequence(f)
        for T in TEMPERATURES:
            out.append({
                "name": f"g6-kfrac{int(round(f * 100)):03d}-{int(T)}k",
                "fraction": f, "n_lys": seq.count("K"), "sequence": seq,
                "temperature": T,
                # Same seed at both temperatures of a sequence.
                "seed": SEED0 + i * 10,
            })
    return out


def init(root: Path, made: list[dict]) -> None:
    for T in TEMPERATURES:
        if not T_MIN <= T <= T_MAX:
            raise ValueError(f"{T} K is outside TEA's {T_MIN:g}-{T_MAX:g} K range")
    n_frames = int(round(NS * 1000.0 / (N_SAVE * 0.01)))
    for c in made:
        sim_dir = root / "simulations" / c["name"]
        sim_dir.mkdir(parents=True, exist_ok=True)
        text = PREPARE.format(
            title=(f"(VPGIG)96 with the first {c['n_lys']} pentads as VPGKG "
                   f"({c['fraction']:.0%})."),
            sequence=c["sequence"], temperature=c["temperature"], gamma=GAMMA,
            tea=True, extrapolate=False, box=BOX_NM, n_save=N_SAVE,
            n_frames=n_frames, seed=c["seed"], ns=NS,
            tea_line=f"CALVADOS-TEA on, gamma {GAMMA:g}.", partition=PARTITION,
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
    made = cells()

    if args.command == "plan":
        print(f"{len(made)} runs, single free chains, {NS:g} ns, {BOX_NM:g} nm, "
              f"TEA gamma {GAMMA:g}, {PARTITION}")
        for c in made:
            print(f"  {c['name']:<24}{c['temperature']:>5.0f} K {len(c['sequence']):>4} aa"
                  f"  {c['n_lys']:>3} VPGKG  seed {c['seed']}")
        return 0

    if args.command == "init":
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
        print(f"would sbatch {len(made)} jobs; re-run with --yes")
        return 0
    for c in made:
        sim_dir = root / "simulations" / c["name"]
        res = subprocess.run(["sbatch", str(sim_dir / "runtime" / "job.sh")],
                             cwd=root, capture_output=True, text=True)
        job = parse_job_id(res.stdout)
        print(f"  {c['name']:<24}{job or 'FAILED: ' + res.stderr.strip()}")
        if job:
            record_simulation(sim_dir, job_id=job, csv_name="tea_lysfrac.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
