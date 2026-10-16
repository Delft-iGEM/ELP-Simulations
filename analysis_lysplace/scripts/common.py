"""Shared paths, run sheets, sequence parsing and design descriptors for the RGD+ lysine-placement
submodule. Nothing here reads a trajectory or a response. Run every script from the repository
root, e.g. ``python analysis_lysplace/scripts/00_inventory.py``.

Construct (read from each run's simulated molecules.fasta, never from the run name alone):
LQ + 4 x [RGD domain LDASTVYAVTGRGDSPASSAASA + ELP block of (VPGXG)n + VP], blocks of 14, 11,
15 and 11 repeats, X = I except four K (the design) and M at repeat 6 of blocks 2 and 4. The run
name encodes the K repeat positions per block: "13R1G13D4" = block 1 repeat 13, block 2 repeat 1,
block 3 repeat 13, block 4 repeat 4; "0R0G0D3+4+5+7" = four K in block 4; 0 = no K in that block.
"""
from __future__ import annotations

import re
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MOD = ROOT / "analysis_lysplace"
TIMESERIES = MOD / "timeseries"
FIGURES = MOD / "figures"
STATS = MOD / "stats"
SIMS = ROOT / "simulations"

SHEETS = [
    "Chr_RGD+_Lysine_Sweep.csv",              # uniform / focused / clustered / random, 40
    "Chr_RGD+_Lysine_sweep_20.csv",           # uniform / focused / random, 20
    "Chr_RGD+_Lysine_sweep_30.csv",           # clustered / uniform / focused / random, 30
    "Chr_RGD+_lysRandom.csv",                 # random, 20
    "Chr_RGD+_lysine_scramble_sweep_20.csv",  # no class suffix; one K per block -> "scramble"
]
DOMAIN = "LDASTVYAVTGRGDSPASSAASA"
RGD = "RGD"
CLASSES = ["clustered", "focused", "random", "uniform", "scramble"]   # fixed order in tables and figures


def slugify(name: str) -> str:
    """Same rule as tools.new_simulation._slugify, which named the folders."""
    slug = "".join(c if c.isalnum() else "-" for c in name.strip().lower())
    return "-".join(filter(None, slug.split("-")))


def read_sheets() -> pd.DataFrame:
    out = []
    for f in SHEETS:
        d = pd.read_csv(StringIO((ROOT / f).read_text(encoding="utf-8-sig")), sep=";", skipinitialspace=True,
                        dtype=str)
        d.columns = [c.strip() for c in d.columns]
        d["Name"] = d["Name"].str.strip()
        d["sheet"] = f
        out.append(d)
    df = pd.concat(out, ignore_index=True)
    df["slug"] = df.Name.map(slugify)
    assert not df.slug.duplicated().any(), "a run name appears twice"
    m = df.Name.str.extract(r"_(uniform|focused|clustered|random)$", expand=False)
    df["design_class"] = m.fillna("scramble")
    return df


def name_positions(name: str) -> list[list[int]]:
    """K repeat positions per block from the run name: '0R2+3G0D11' -> [[], [2, 3], [], [11]]."""
    core = name.split("_")[0]
    parts = re.split(r"[RGD]", core)
    if len(parts) != 4:
        raise ValueError(f"cannot parse {name!r}")
    return [[] if p == "0" else [int(x) for x in p.split("+")] for p in parts]


def parse_fasta(path: Path) -> dict:
    """Block layout and K / M repeat positions from the simulated sequence."""
    seq = "".join(path.read_text().split("\n", 1)[1].split())
    starts = [m.start() for m in re.finditer(DOMAIN, seq)]
    if len(starts) != 4:
        raise ValueError(f"{path}: expected 4 RGD domains, found {len(starts)}")
    blocks = []
    for b, st in enumerate(starts):
        lo = st + len(DOMAIN)
        hi = starts[b + 1] if b + 1 < 4 else len(seq)
        body = seq[lo:hi]
        if not body.endswith("VP") or (len(body) - 2) % 5:
            raise ValueError(f"{path}: block {b + 1} is not (VPGXG)n + VP")
        reps = [body[i:i + 5] for i in range(0, len(body) - 2, 5)]
        if any(r[:3] != "VPG" or r[4] != "G" for r in reps):
            raise ValueError(f"{path}: block {b + 1} has a non-VPGXG repeat")
        blocks.append({"start": lo, "n_repeats": len(reps),
                       "K": [i + 1 for i, r in enumerate(reps) if r[3] == "K"],
                       "M": [i + 1 for i, r in enumerate(reps) if r[3] == "M"]})
    k_res = [i for i, c in enumerate(seq) if c == "K"]                 # 0-based residue index (the X of VPGKG)
    rgd_res = [i for st in starts for i in range(st + DOMAIN.index(RGD), st + DOMAIN.index(RGD) + 3)]
    return {"seq": seq, "length": len(seq), "blocks": blocks, "k_res": k_res, "rgd_res": rgd_res}


def descriptors(k_res: list[int], blocks: list[dict], length: int) -> dict:
    """Design descriptors of one chain's lysine placement (residue indices along the chain)."""
    k = np.array(sorted(k_res))
    gaps = np.diff(k)
    per_block = [len(b["K"]) for b in blocks]
    nb = int(sum(c > 0 for c in per_block))
    adjacent = nb == 1 and gaps.max() <= 10            # consecutive repeats, allowing one M in between
    structure = ("one_per_block" if per_block == [1, 1, 1, 1] else
                 "single_block_adjacent" if adjacent else "single_block_spread" if nb == 1 else
                 f"{nb}_blocks")
    pairs = np.abs(k[:, None] - k[None, :])[np.triu_indices(len(k), 1)]
    return {"n_lys": len(k), "n_blocks_with_K": nb, "structure": structure,
            "k_per_block": "-".join(map(str, per_block)),
            "min_sep_res": int(gaps.min()), "mean_sep_res": float(gaps.mean()), "max_sep_res": int(gaps.max()),
            "span_res": int(k[-1] - k[0]), "first_k_rel": k[0] / (length - 1), "last_k_rel": k[-1] / (length - 1),
            "centroid_rel": float(k.mean() / (length - 1)),
            "n_pairs_lt_30res": int((pairs < 30).sum()),
            "harmonic_sep": float(1 / np.mean(1 / pairs))}      # pair separations, dominated by the closest pairs
