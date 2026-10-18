"""Shared paths, run-sheet parsing and design-variable derivation for the SubM4 analysis.

Nothing here reads a trajectory or a response. Every script in this folder is run
from the repository root, e.g. ``python analysis/scripts/00_inventory.py``.
"""
from __future__ import annotations

import re
import subprocess
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analysis"
TIMESERIES = ANALYSIS / "timeseries"
FIGURES = ANALYSIS / "figures"
SIMS = ROOT / "simulations"

# Run sheets that were submitted from this checkout (Chris's DelftBlue account).
LOCAL_SHEETS = [
    "Chr_sub4_RUN_1.csv",   # anchors BASE + A01-A09 (= submodule4_anchor_replicates_part1_30)
    "CHR_sub4_RUN_2.csv",   # sub4-extra2 (20) + sub4-extra3 (10)
    "CHR_sub4_RUN_3.csv",   # sub4-extra4 (30)
    "CHR_sub4_RUN_4.csv",   # sub4-extra5 01-14 (14 of the planned 20)
    "rerun_A10-A19_nonreactive_30.csv",   # A10-A19 x3 rerun with crosslinking off (-nr-repN); unfinished ones stay "incomplete"
]
# Run sheets submitted by a teammate (Titus, commit 0d88132, merged in 6725a2d). Their
# trajectories live on /scratch/tnartey, which is not readable from this account; what
# reached the repository is his per-run notebook output in simulations/<name>/analysis/
# (crosslink_contacts_K.txt, load_summary.txt), which gives the team-definition-A
# contact count at 12 A after the same 100 ns burn-in, but not the free-only one.
TEAMMATE_SHEETS = [
    "submodule4_extra_20.csv",                    # sub4-extra-* (20), non-reactive
    "submodule4_anchor_replicates_part2_30.csv",  # anchors A10-A19 x3, reactive (excluded)
]

# Design ranges from the SubM4 protocol (used for coding to [-1, +1]).
RANGES = {
    "concentration": (0.01, 0.06),     # chains / nm^2
    "chain_length": (50, 150),         # VPGXG pentarepeats
    "spacing_pct": (5.0, 25.0),        # lysine-to-lysine distance as % of chain length
    "ionic_strength": (0.05, 0.50),    # M
}
BASE = {"concentration": 0.04, "chain_length": 72, "spacing_pct": 100 * 15 / 72, "ionic_strength": 0.20}


def slugify(name: str) -> str:
    """Same rule as tools.new_simulation._slugify, which named the folders."""
    slug = "".join(c if c.isalnum() else "-" for c in name.strip().lower())
    return "-".join(filter(None, slug.split("-")))


def read_sheet(path: str, ref: str | None = None) -> pd.DataFrame:
    text = (subprocess.run(["git", "show", f"{ref}:{path}"], cwd=ROOT, capture_output=True,
                           text=True, check=True).stdout if ref else (ROOT / path).read_text())
    df = pd.read_csv(StringIO(text.lstrip("﻿")), sep=";", dtype={"crosslink distance": str,
                                                                    "crosslink prob": str})
    df["sheet"] = f"{ref}:{path}" if ref else path
    return df


def parse_sequence(seq: str) -> dict:
    """Design variables from the shorthand, e.g. 'I13 K1 I14 K1 I14 K1 I14 K1 I13'.

    spacing_repeats is the lysine-to-lysine distance in repeats (gap + 1), which is
    what the run names call S (I3 K1 I12 K1 ... -> S13). The protocol requires the
    four lysines to be equally spaced; that is checked, not assumed.
    """
    tokens = [(t[0], int(t[1:])) for t in seq.split()]
    length = sum(n for _, n in tokens)
    pos, k_pos = 0, []
    for letter, n in tokens:
        if letter == "K":
            k_pos.extend(range(pos, pos + n))
        elif letter != "I":
            raise ValueError(f"unexpected token {letter} in {seq!r}")
        pos += n
    gaps = {b - a for a, b in zip(k_pos, k_pos[1:])}
    if len(k_pos) != 4 or len(gaps) != 1:
        raise ValueError(f"expected 4 equally spaced lysines in {seq!r}")
    s = gaps.pop()
    return {"chain_length": length, "n_lys": len(k_pos), "spacing_repeats": s,
            "spacing_pct": 100.0 * s / length, "lead_repeats": k_pos[0],
            "trail_repeats": length - 1 - k_pos[-1]}


def parse_contacts_summary(path: Path) -> dict:
    """The numbers in the notebook's crosslink_contacts_K.txt (teammate runs)."""
    text = path.read_text()
    num = r"([\d,.]+)"

    def grab(pattern, cast=float):
        m = re.search(pattern, text)
        if not m:
            raise ValueError(f"{path}: no match for {pattern!r}")
        return [cast(g.replace(",", "")) for g in m.groups()]

    (cut,) = grab(r"within ([\d.]+) Å")
    done, total, first = grab(r"frames analysed:\s+(\d+)/(\d+) \(from frame (\d+)\)", int)
    n_intra, n_inter, n_both_pinned = grab(r"(\d+) intra-chain \+ (\d+) inter-chain = \d+\s+\((\d+) pairs dropped", int)
    hits_intra, hits_inter = grab(rf"hits \(pair-frames\)\s+{num}\s+{num}")
    ever_intra, ever_inter = grab(r"pairs ever in contact\s+(\d+)/\d+\s+(\d+)/\d+", int)
    return {"cutoff_A": cut, "frames_analysed": done, "frames_total": total, "first_frame": first,
            "pairs_intra": n_intra, "pairs_inter": n_inter, "pairs_dropped_both_pinned": n_both_pinned,
            "hits_intra": hits_intra, "hits_inter": hits_inter,
            "pairs_ever_intra": ever_intra, "pairs_ever_inter": ever_inter}


def parse_load_summary(path: Path) -> dict:
    """Chains, residues, box (nm) and observed headroom from the notebook's load_summary.txt."""
    text = path.read_text()
    m = re.search(r"(\d+) chains x (\d+) residues/chain, box = \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)", text)
    h = re.search(r"highest bead ever at z = ([\d.]+) Å.*headroom ([\d.]+) Å", text)
    if not (m and h):
        raise ValueError(f"{path}: unexpected format")
    return {"n_chains": int(m.group(1)), "n_residues": int(m.group(2)),
            "box_x_nm": float(m.group(3)) / 10, "box_y_nm": float(m.group(4)) / 10,
            "box_z_nm": float(m.group(5)) / 10, "z_max_nm": float(h.group(1)) / 10,
            "box_headroom_nm": float(h.group(2)) / 10}


def anchor_fields(name: str) -> dict:
    m = re.match(r"sub4-anchor-(BASE-CASE|A\d\d)-.*-rep(\d)$", name, re.I)
    if not m:
        return {"anchor_id": None, "replicate": None}
    aid = m.group(1).upper()
    return {"anchor_id": "BASE" if aid.startswith("BASE") else aid, "replicate": int(m.group(2))}


def batch_of(name: str) -> str:
    m = re.match(r"sub4-(anchor|extra\d?|main)", name, re.I)
    return m.group(1).lower() if m else "other"


def coded(df: pd.DataFrame) -> pd.DataFrame:
    """Design variables scaled to [-1, +1] over the protocol ranges (the frozen coding)."""
    out = pd.DataFrame(index=df.index)
    for k, (lo, hi) in RANGES.items():
        out[k] = (df[k] - lo) / (hi - lo) * 2 - 1
    return out


def coded_log(df: pd.DataFrame) -> pd.DataFrame:
    """Natural log of each design variable, scaled to [-1, +1] over the log of its range.

    Not frozen: chosen on 2026-10-16 after the frozen model's lack of fit was seen."""
    out = pd.DataFrame(index=df.index)
    for k, (lo, hi) in RANGES.items():
        out[k] = (np.log(df[k]) - np.log(lo)) / (np.log(hi) - np.log(lo)) * 2 - 1
    return out


CODINGS = {"linear": coded, "log": coded_log}
