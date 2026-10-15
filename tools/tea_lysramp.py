"""Lysine ramps on the 357-aa RGD construct: 4 -> 49 VPGKG in 10 steps.

    .venv/bin/python tools/tea_lysramp.py plan
    .venv/bin/python tools/tea_lysramp.py init       # writes prepare.py and runs it
    .venv/bin/python tools/tea_lysramp.py submit --yes

40 runs: 2 series x 10 lysine counts x {290, 370} K, single free chains,
400 ns, CALVADOS-TEA at gamma 6 - the conditions of tools/tea_g6.py.

Sequences
---------
The original is 4 RGD blocks + 51 pentads + "VP" (tools/tea_g6.py). The RGD
blocks, the "VP" tail and the two VPGMG (pentads 19 and 45) are never touched,
so every sequence is 357 aa with the RGD at the same residues. Of the 49 other
pentads, n are VPGKG and the rest VPGIG, n = 4, 9, 14, ..., 49 (49 = every
VPGIG replaced).

  kspc  spaced: the n lysines spread as evenly as possible over the 49
        positions (midpoint rule), so n=4 is not the original arrangement.
  kend  ends: the original's 4 (pentads 40, 41, 43, 46) are kept and each new
        lysine fills the free position nearest a chain end, alternating N and
        C terminus, so the middle stays VPGIG and is encroached from both
        sides. n=4 is the original verbatim.

Box
---
150 nm. The contour is 357 x 0.38 = 136 nm; with the 4 nm cutoff on each side
that is 144 nm, so even fully stretched (the +49 chain) cannot meet its image.

Hardware
--------
gpu-v100, 8 CPUs x 4 GB. 400 ns of a 357-aa g6 chain took 1000-1020 s on V100
(timing.yaml of g6-orig, g6-k8-*), and 480 aa took ~1060 s, so cost is not
bead-bound here. 0:40:00 is ~2.3x that.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

from tools.new_simulation import _project_root
from tools.tea import DEFAULT_GAMMA, T_MAX, T_MIN
from tools.tea_g6 import k8_sequence
from tools.tea_sequences import original
from tools.tea_temps import PREPARE

GAMMA = DEFAULT_GAMMA
TEMPERATURES = [290.0, 370.0]
N_PENTADS = 51
MET = (19, 45)
ORIG_LYS = [40, 41, 43, 46]
FREE = [i for i in range(N_PENTADS) if i not in MET]       # 49 positions
COUNTS = [4 + 5 * k for k in range(10)]                    # 4 .. 49
NS = 400.0
N_SAVE = 7000
BOX_NM = 150.0
WALLTIME = "0:40:00"
PARTITION = "gpu-v100"
CPUS = "8"
MEM_PER_CPU = "4000"
SEED0 = 20261014


def spaced(n: int) -> list[int]:
    return [FREE[int((i + 0.5) * len(FREE) / n)] for i in range(n)]


def ends(n: int) -> list[int]:
    pos = list(ORIG_LYS)
    left = [p for p in FREE if p not in pos]
    from_n = True
    while len(pos) < n:
        pos.append(left.pop(0) if from_n else left.pop())
        from_n = not from_n
    return sorted(pos)


SERIES = {"kspc": spaced, "kend": ends}


def cells(root: Path) -> list[dict]:
    orig = original(root)
    rgd = [m.start() for m in re.finditer("RGD", orig)]
    out = []
    for s, (series, place) in enumerate(SERIES.items()):
        for i, n in enumerate(COUNTS):
            pos = place(n)
            seq = k8_sequence(pos, root)
            if (len(set(pos)) != n or set(pos) & set(MET) or len(seq) != len(orig)
                    or seq.count("VPGKG") != n or seq.count("VPGMG") != 2
                    or [m.start() for m in re.finditer("RGD", seq)] != rgd):
                raise ValueError(f"{series} n={n}: bad sequence")
            for T in TEMPERATURES:
                out.append({
                    "name": f"g6-{series}{n:02d}-{int(T)}k", "series": series,
                    "n_lys": n, "positions": pos, "sequence": seq, "temperature": T,
                    # Same seed at both temperatures of a sequence.
                    "seed": SEED0 + (s * 10 + i) * 10,
                })
    if ends(4) != ORIG_LYS or k8_sequence(ORIG_LYS, root) != orig:
        raise ValueError("kend04 is not the original")
    return out


def init(root: Path, made: list[dict]) -> None:
    for T in TEMPERATURES:
        if not T_MIN <= T <= T_MAX:
            raise ValueError(f"{T} K is outside TEA's {T_MIN:g}-{T_MAX:g} K range")
    n_frames = int(round(NS * 1000.0 / (N_SAVE * 0.01)))
    for c in made:
        sim_dir = root / "simulations" / c["name"]
        sim_dir.mkdir(parents=True, exist_ok=True)
        how = "spaced evenly" if c["series"] == "kspc" else "filled in from both ends"
        text = PREPARE.format(
            title=(f"357-aa RGD construct with {c['n_lys']} VPGKG, {how}, "
                   f"at pentads {c['positions']}."),
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
    made = cells(root)

    if args.command == "plan":
        print(f"{len(made)} runs, single free chains, {NS:g} ns, {BOX_NM:g} nm, "
              f"TEA gamma {GAMMA:g}, {PARTITION} {WALLTIME}")
        for c in made[::2]:
            seq = c["sequence"]
            pent = "".join("K" if i in c["positions"] else ("M" if i in MET else ".")
                           for i in range(N_PENTADS))
            print(f"  {c['name'][:-5]:<14}{len(seq):>4} aa {c['n_lys']:>3} K  {pent}")
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
            record_simulation(sim_dir, job_id=job, csv_name="tea_lysramp.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
