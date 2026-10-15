"""285 / 370 K single-chain set at gamma 6, plus ten 8-lysine RGD variants.

    .venv/bin/python tools/tea_g6.py plan
    .venv/bin/python tools/tea_g6.py init       # writes prepare.py and runs it
    .venv/bin/python tools/tea_g6.py submit --yes

54 runs, all single free chains, 400 ns, CALVADOS-TEA on at gamma 6, at 285 K
and 370 K. Both temperatures are inside TEA's 280-380 K fit range, so nothing
here is extrapolated.

  16  lys      best 4 and worst 4 lysine arrangements of the 298-aa scaffold,
               ranked by inter-chain share of crosslink hits
               (contacts-lysine-arrangement.csv / lysine-arrangement-results.xlsx).
  16  rgd      best 4 and worst 4 RGD arrangements of the 357-aa construct,
               ranked by 23-mer SASA over the 24 free-chain runs
               (figures/free-rgdarr/sasa-24-arrangements.csv).
   2  orig     the 357-aa original RGD construct, verbatim.
  20  k8       ten arrangements of the original with 8 VPGKG instead of 4.

Each sequence uses the same seed at both temperatures, so a pair differs only in
lambda(T) and kT.

The 8-lysine variants
---------------------
The original is 4 RGD blocks + 51 pentads + "VP"; the pentads sit in four gaps
of 14, 11, 15 and 11 after blocks 1-4 (pentad indices 0-13, 14-24, 25-39,
40-50). The variants turn 4 more VPGIG into VPGKG, so all ten have exactly the
same composition (41 VPGIG, 8 VPGKG, 2 VPGMG) and the same RGD blocks at the
same residues. VPGMG stays at pentads 19 and 45 as in the original. Only the
eight lysine positions differ between them - checked by `k8_sequences`.

Box
---
80 nm, as the 270/350 K set. lambda at 285 K, gamma 6 is within 0.03 of lambda
at 270 K, gamma 3 for every ELP residue (V 0.066 vs 0.055, I 0.280 vs 0.255),
and those chains fitted in 80 nm.
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
from tools.tea_sequences import decompose, original
from tools.tea_temps import PREPARE, _sequence_of

GAMMA = DEFAULT_GAMMA
TEMPERATURES = [285.0, 370.0]
NS = 400.0
N_SAVE = 7000                      # 70 ps/frame, as every TEA set so far
BOX_NM = 80.0

# Measured 969-1174 ns/h for a single free chain on an A100: 400 ns is ~25 min
# at the slow end. Roughly 2x headroom.
WALLTIME = "0:50:00"
PARTITION = "gpu-a100"
CPUS = "18"
SEED0 = 20261005

# Ranked by inter_share_of_hits, contacts-lysine-arrangement.csv.
LYS_BEST = ["lysarr12-even-4", "lysarr13-even-4-wide",
            "lysarr19-blockb", "lysarr20-blockc"]
LYS_WORST = ["lysarr02-blk4-mid", "lysarr03-blk4-cterm",
             "lysarr09-blk13-nc", "lysarr01-blk4-nterm"]
# Ranked by motif_23mer_nm2, figures/free-rgdarr/sasa-24-arrangements.csv.
RGD_BEST = ["free-rgdarr19-3n-1c", "free-rgdarr15-pairs-ends",
            "free-rgdarr14-tight-c", "free-rgdarr22-c-heavy"]
RGD_WORST = ["free-rgdarr16-pairs-mid", "free-rgdarr10-clust-q1",
             "free-rgdarr03-even-inset", "free-rgdarr23-staircase"]

# Pentad indices (0-50) of the eight VPGKG in each variant. 19 and 45 are VPGMG.
K8 = {
    # The original's lysine motif (gap 4: offsets 0, 1, 3, 6 around the VPGMG
    # at offset 5) kept, and copied into gap 2, whose VPGMG sits at the same
    # offset. The closest variant to the original.
    "01-orig-dup": [14, 15, 17, 20, 40, 41, 43, 46],
    # Single lysines spread as evenly as the chain allows.
    "02-even": [3, 9, 16, 22, 28, 34, 41, 47],
    # Two per gap, at each gap's thirds.
    "03-gap-thirds": [4, 9, 17, 22, 30, 35, 43, 48],
    # One immediately on each side of every RGD block, plus both chain ends.
    "04-rgd-flank": [0, 13, 14, 24, 25, 39, 40, 50],
    # A KK pair at the centre of each gap (stepping round the VPGMG).
    "05-pairs-centre": [6, 7, 18, 20, 32, 33, 44, 46],
    # One block of eight at the C terminus, the original's region.
    "06-blk8-cterm": [40, 41, 42, 43, 44, 46, 47, 48],
    # One block of eight in the middle of the longest gap.
    "07-blk8-mid": [28, 29, 30, 31, 32, 33, 34, 35],
    # Four at each chain end.
    "08-blk44-ends": [0, 1, 2, 3, 47, 48, 49, 50],
    # Alternating with VPGIG through the middle gap: spaced, but confined.
    "09-alt-mid": [25, 27, 29, 31, 33, 35, 37, 39],
    # Uniform random, fixed seed - the no-design baseline.
    "10-random": sorted(random.Random(SEED0).sample(
        [i for i in range(51) if i not in (19, 45)], 8)),
}


def k8_sequence(positions: list[int], root: Path | None = None) -> str:
    """The original with VPGKG at exactly these pentads and VPGIG elsewhere."""
    seq = original(root)
    blocks, pentads, tail = decompose(seq)
    new = ["VPGKG" if i in positions else ("VPGMG" if p == "VPGMG" else "VPGIG")
           for i, p in enumerate(pentads)]
    # Rebuild in the original's order: gaps of 14, 11, 15, 11 after blocks 1-4.
    out, i = [], 0
    for block, n in zip(blocks, (14, 11, 15, 11)):
        out.append(block)
        out.extend(new[i:i + n])
        i += n
    return "".join(out) + tail


def k8_sequences(root: Path | None = None) -> dict[str, str]:
    orig = original(root)
    rgd = [m.start() for m in re.finditer("RGD", orig)]
    out = {}
    for label, pos in K8.items():
        if len(set(pos)) != 8 or {19, 45} & set(pos) or not all(0 <= p <= 50 for p in pos):
            raise ValueError(f"k8-{label}: bad positions {pos}")
        seq = k8_sequence(pos, root)
        if len(seq) != len(orig) or [m.start() for m in re.finditer("RGD", seq)] != rgd:
            raise ValueError(f"k8-{label}: length or RGD positions moved")
        if seq.count("VPGKG") != 8 or seq.count("VPGMG") != 2:
            raise ValueError(f"k8-{label}: wrong pentad counts")
        out[label] = seq
    if len(set(out.values())) != len(out):
        raise ValueError("two k8 variants are the same sequence")
    if len({frozenset(Counter(s).items()) for s in out.values()}) != 1:
        raise ValueError("k8 variants do not share one composition")
    return out


def cells(root: Path) -> list[dict]:
    seqs: list[tuple[str, str, str, str]] = []      # group, label, title, sequence
    for rank, source in [("best", s) for s in LYS_BEST] + [("worst", s) for s in LYS_WORST]:
        label = source.replace("lysarr", "lys")
        seqs.append(("lys", label, f"Lysine arrangement {source} ({rank} 4).",
                     _sequence_of(source, root)))
    for rank, source in [("best", s) for s in RGD_BEST] + [("worst", s) for s in RGD_WORST]:
        label = source.replace("free-rgdarr", "rgd")
        seqs.append(("rgd", label, f"RGD arrangement {source} ({rank} 4).",
                     _sequence_of(source, root)))
    seqs.append(("orig", "orig", "The 357-aa original RGD construct, verbatim.",
                 original(root)))
    for label, seq in k8_sequences(root).items():
        seqs.append(("k8", f"k8-{label}",
                     f"Original RGD construct with 8 VPGKG, at pentads {K8[label]}.", seq))

    out = []
    for i, (group, label, title, seq) in enumerate(seqs):
        for T in TEMPERATURES:
            out.append({
                "group": group, "name": f"g6-{label}-{int(T)}k", "title": title,
                "sequence": seq, "temperature": T,
                # Same seed at both temperatures of a sequence.
                "seed": SEED0 + i * 10,
            })
    return out


def init(root: Path, made: list[dict]) -> None:
    for T in TEMPERATURES:
        if not T_MIN <= T <= T_MAX:
            raise ValueError(f"{T} K is outside TEA's {T_MIN:g}-{T_MAX:g} K range")
    n_frames = int(round(NS * 1000.0 / (N_SAVE * 0.01)))
    for cell in made:
        sim_dir = root / "simulations" / cell["name"]
        sim_dir.mkdir(parents=True, exist_ok=True)
        (sim_dir / "prepare.py").write_text(PREPARE.format(
            title=cell["title"], sequence=cell["sequence"],
            temperature=cell["temperature"], gamma=GAMMA, tea=True, extrapolate=False,
            box=BOX_NM, n_save=N_SAVE, n_frames=n_frames, seed=cell["seed"], ns=NS,
            tea_line=f"CALVADOS-TEA on, gamma {GAMMA:g}.", partition=PARTITION,
            walltime=WALLTIME, cpus=CPUS))
        subprocess.run([str(root / ".venv/bin/python"), "-m",
                        f"simulations.{cell['name']}.prepare"],
                       cwd=root, check=True, capture_output=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "init", "submit"])
    parser.add_argument("--group", action="append", choices=["lys", "rgd", "orig", "k8"])
    parser.add_argument("--yes", action="store_true")
    args = parser.parse_args(argv)
    root = _project_root()
    made = [c for c in cells(root) if not args.group or c["group"] in args.group]

    if args.command == "plan":
        print(f"{len(made)} runs, single free chains, {NS:g} ns, {BOX_NM:g} nm, "
              f"TEA gamma {GAMMA:g}")
        for c in made:
            print(f"  {c['name']:<34}{c['group']:<6}{c['temperature']:>5.0f} K "
                  f"{len(c['sequence']):>4} aa  seed {c['seed']}  "
                  f"K@{[i for i, ch in enumerate(c['sequence']) if ch == 'K']}")
        print(f"~{len(made) * NS / 969.0:.1f} GPU-hours at the measured slow end")
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
        print(f"  {c['name']:<34}{job or 'FAILED: ' + res.stderr.strip()}")
        if job:
            record_simulation(sim_dir, job_id=job, csv_name="tea_g6.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
