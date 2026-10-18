"""Random arrangements of the 357-aa RGD construct at 285 / 370 K.

    .venv/bin/python tools/tea_random.py plan
    .venv/bin/python tools/tea_random.py init       # writes prepare.py and runs it
    .venv/bin/python tools/tea_random.py submit --yes

40 runs, single free chains, 400 ns, CALVADOS-TEA gamma 6, the conditions of the
g6-k8 set (tools/tea_g6.py) so they compare with it directly:

  k8r  10 random placements of 8 VPGKG, the RGD blocks fixed. Built exactly like
       the g6-k8 designs (tea_g6.k8_sequence): 8 of the 49 pentads that are not
       VPGMG, drawn uniformly. These give the spread that random placement alone
       produces, against which the ten designed k8 arrangements can be judged.
  rk   10 random arrangements of both the 4 RGD blocks and the 4 VPGKG. The
       chain is written as the RGD sweep wrote it (batches/rgd-arrangement-sweep.csv):
       "LQ" + 55 units, each "VP" + something - 4 RGD units (VP + the 23-mer),
       51 pentads (45 VPGIG, 4 VPGKG, 2 VPGMG). The two VPGMG stay at the units
       they occupy in the original (21 and 49); the 4 RGD and 4 VPGKG are drawn
       uniformly from the other 53 units. Every sequence is a permutation of the
       original's 357 residues - checked.

The draws use a fixed seed, so `plan` always shows the same sequences, and no
draw repeats another or a designed g6-k8 sequence.

Box 80 nm, gpu-v100, 8 CPUs x 4 GB, 0:40:00 - a g6 357-aa chain took
1000-1020 s per 400 ns on V100.
"""

from __future__ import annotations

import argparse
import random
import re
import subprocess
from collections import Counter
from pathlib import Path

from tools.new_simulation import _project_root
from tools.tea import DEFAULT_GAMMA, T_MAX, T_MIN
from tools.tea_g6 import K8, k8_sequence
from tools.tea_sequences import original
from tools.tea_temps import PREPARE

GAMMA = DEFAULT_GAMMA
TEMPERATURES = [285.0, 370.0]
N_RANDOM = 10
NS = 400.0
N_SAVE = 7000
BOX_NM = 80.0
WALLTIME = "0:40:00"
PARTITION = "gpu-v100"
CPUS = "8"
MEM_PER_CPU = "4000"
SEED0 = 20261015
DRAW_SEED = 20261015

PENTAD_FREE = [i for i in range(51) if i not in (19, 45)]   # not VPGMG
RGD_UNIT = "VPLDASTVYAVTGRGDSPASSAASA"
N_UNITS = 55
MET_UNITS = (21, 49)


def k8_random(root: Path) -> list[tuple[list[int], str]]:
    rng = random.Random(DRAW_SEED)
    designed = {k8_sequence(pos, root) for pos in K8.values()}
    out, seen = [], set()
    while len(out) < N_RANDOM:
        pos = sorted(rng.sample(PENTAD_FREE, 8))
        seq = k8_sequence(pos, root)
        if seq in seen or seq in designed:
            continue
        seen.add(seq)
        out.append((pos, seq))
    return out


def rk_sequence(rgd: list[int], lys: list[int]) -> str:
    units = []
    for u in range(N_UNITS):
        units.append(RGD_UNIT if u in rgd else "VPGKG" if u in lys
                     else "VPGMG" if u in MET_UNITS else "VPGIG")
    return "LQ" + "".join(units)


def rk_random(root: Path) -> list[tuple[list[int], list[int], str]]:
    rng = random.Random(DRAW_SEED + 1)
    orig = original(root)
    free = [u for u in range(N_UNITS) if u not in MET_UNITS]
    out, seen = [], set()
    while len(out) < N_RANDOM:
        pick = rng.sample(free, 8)
        rgd, lys = sorted(pick[:4]), sorted(pick[4:])
        seq = rk_sequence(rgd, lys)
        if seq in seen:
            continue
        if len(seq) != len(orig) or Counter(seq) != Counter(orig):
            raise ValueError(f"rk {rgd} {lys}: not a permutation of the original")
        if len(re.findall("RGD", seq)) != 4 or seq.count("VPGKG") != 4:
            raise ValueError(f"rk {rgd} {lys}: wrong RGD or lysine count")
        seen.add(seq)
        out.append((rgd, lys, seq))
    return out


def cells(root: Path) -> list[dict]:
    seqs = []
    for i, (pos, seq) in enumerate(k8_random(root), 1):
        seqs.append((f"k8r{i:02d}", f"Original RGD construct with 8 VPGKG at random "
                     f"pentads {pos} (draw {i}).", seq))
    for i, (rgd, lys, seq) in enumerate(rk_random(root), 1):
        seqs.append((f"rk{i:02d}", f"357-aa construct, random arrangement {i}: RGD at "
                     f"units {rgd}, VPGKG at units {lys}, VPGMG at {list(MET_UNITS)}.", seq))
    out = []
    for i, (label, title, seq) in enumerate(seqs):
        for T in TEMPERATURES:
            out.append({"name": f"g6-{label}-{int(T)}k", "label": label, "title": title,
                        "sequence": seq, "temperature": T,
                        # Same seed at both temperatures of a sequence.
                        "seed": SEED0 + i * 10})
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
            title=c["title"], sequence=c["sequence"], temperature=c["temperature"],
            gamma=GAMMA, tea=True, extrapolate=False, box=BOX_NM, n_save=N_SAVE,
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


def _map(seq: str) -> str:
    """One character per unit: R RGD block, K VPGKG, M VPGMG, . VPGIG."""
    out, i = [], 0
    while i < len(seq):
        if seq.startswith("LQLDAS", i):           # the original's first block
            out.append("R"); i += 25
        elif seq.startswith("LQ", i) and i == 0:
            i += 2
        elif seq.startswith(RGD_UNIT, i):
            out.append("R"); i += 25
        elif seq.startswith("VP", i) and len(seq) - i == 2:
            i += 2
        else:
            out.append({"K": "K", "M": "M"}.get(seq[i + 3], ".")); i += 5
    return "".join(out)


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
        print(f"  {'original':<8}{len(original(root)):>4} aa  {_map(original(root))}")
        for c in made[::2]:
            print(f"  {c['label']:<8}{len(c['sequence']):>4} aa  {_map(c['sequence'])}")
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
        print(f"  {c['name']:<20}{job or 'FAILED: ' + res.stderr.strip()}")
        if job:
            record_simulation(sim_dir, job_id=job, csv_name="tea_random.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
