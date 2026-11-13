"""Sequence parsing, run-length presets and radius-of-gyration analysis.

Kept free of OpenMM/CALVADOS imports so the UI server can use it and start fast.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np

from tools.elibpy import build_sequence_with_features

RESIDUES_CSV = Path(__file__).resolve().parent.parent / "residues_CALVADOS2.csv"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MIN_RESIDUES = 5

# CALVADOS adds these to the terminal beads with charge_termini="both".
N_TERM_MASS = 2.0
C_TERM_MASS = 16.0

TIMESTEP_PS = 0.01  # fixed by calvados.sim.Sim.simulate
EQUILIBRATION_FRACTION = 0.1
N_BLOCKS = 10

PRESETS = {
    # 2/3 of outreach_sim_1's 21 ns at its 70 ps frame spacing, whatever the chain length.
    "fast": {"label": "Fast", "frames": 200, "wfreq": 7000},
    "quick": {"label": "Quick", "frames": 100},
    "standard": {"label": "Standard", "frames": 300},
    "long": {"label": "Long", "frames": 1000},
}


# Blocks that count towards the block-variation bonus.
# ELP blocks: guest residue + repeat count, e.g. K10 = 10 x VPGKG.
ELP_GUESTS = "IKEFACVQW"
ELP_REPEATS = (5, 10, 20)
# Named blocks. Typing the name inserts `seq`. Adhesion domains carry the tissue they bind,
# matching the Rg regions on the race chart (race.js RG_REGIONS).
NAMED_BLOCKS = {
    "RGD": {"seq": "RGD", "kind": "adhesion", "tissue": "Bone"},
    "IKVAV": {"seq": "IKVAV", "kind": "adhesion", "tissue": "Skin"},
    # Cardin-Weintraub heparin-binding motif (X hydropathic, B basic); AKKARA is one instance.
    "XBBXBX": {"seq": "AKKARA", "kind": "adhesion", "tissue": "Brain"},
    "KR12": {"seq": "KRIVQRIKDFLR", "kind": "amp", "tissue": None},  # LL-37 residues 18-29
    "SAAP": {"seq": "LKRVWKRVFKLLKRYWRQLKKPVR", "kind": "amp", "tissue": None},  # SAAP-148
}
BLOCK_ALIASES = {"SAAP148": "SAAP"}
# Bonus by number of different block types; 7 or more gets the full 50.
VARIATION_POINTS = (0, 0, 5, 10, 20, 35, 45, 50)


class SequenceError(ValueError):
    """The user's input can't be turned into a sequence; the message is user-facing."""


def tokenize(text: str) -> list[str]:
    lines = [line for line in text.splitlines() if not line.lstrip().startswith(">")]
    tokens = [t for t in re.split(r"[\s,\-]+", " ".join(lines).upper()) if t]
    return [BLOCK_ALIASES.get(t, t) for t in tokens]


def block_type(token: str) -> str | None:
    """The block type a token counts as for the variation bonus, or None if it isn't a listed block."""
    if token in NAMED_BLOCKS:
        return token
    match = re.fullmatch(r"([A-Z])(\d+)", token)
    if match and match[1] in ELP_GUESTS and int(match[2]) in ELP_REPEATS:
        return match[1]
    return None


def analyze_blocks(text: str) -> dict:
    """Which listed blocks the input uses, and the variation bonus they earn."""
    types, unscored = [], []
    for token in tokenize(text):
        kind = block_type(token)
        if kind is None:
            if token not in unscored:
                unscored.append(token)
        elif kind not in types:
            types.append(kind)
    tissues = sorted({NAMED_BLOCKS[t]["tissue"] for t in types if t in NAMED_BLOCKS and NAMED_BLOCKS[t]["tissue"]})
    return {
        "types": types,
        "unscored": unscored,
        "adhesion_tissues": tissues,
        "variation_pts": VARIATION_POINTS[min(len(types), len(VARIATION_POINTS) - 1)],
    }


def parse_sequence(text: str) -> str:
    """Turn user input into a one-letter sequence.

    Accepts FASTA, raw sequence, or the elibpy shorthand, separated by spaces,
    commas or hyphens: ``I10 K10 I5 KR12 RGD``. ``K10`` means 10 ``VPGKG`` repeats,
    and named blocks (NAMED_BLOCKS) insert their sequence.
    """
    tokens = tokenize(text)
    if not tokens:
        raise SequenceError("Enter a sequence.")

    for token in tokens:
        if re.fullmatch(r"[A-Z]0+", token):
            raise SequenceError(f"'{token}' has zero repeats.")

    components = [(f"({NAMED_BLOCKS[t]['seq']})" if t in NAMED_BLOCKS else t, t) for t in tokens]
    seq = build_sequence_with_features(components)["seq"]

    invalid = sorted(set(seq) - AMINO_ACIDS)
    if "Z" in invalid:
        raise SequenceError(
            "'Z' is the surface-anchor bead, which pins the chain in place. "
            "This tool simulates free chains; use a prepare.py for anchored ones."
        )
    if invalid:
        raise SequenceError(
            f"Not understood: {', '.join(invalid)}. Use the 20 standard amino acids, "
            "and put spaces between blocks, like 'I10 K10 RGD'."
        )
    if len(seq) < MIN_RESIDUES:
        raise SequenceError(f"Sequence needs at least {MIN_RESIDUES} residues.")
    return seq


def frame_interval(n_residues: int) -> int:
    """Steps between saved frames. Longer chains decorrelate more slowly."""
    if n_residues < 150:
        return 7000
    return math.ceil(3e-4 * n_residues**2) * 1000


def box_length(n_residues: int) -> float:
    """Cubic box edge (nm) large enough that a chain doesn't meet its own image.

    About 6x a rough Rg plus twice the 4 nm Debye-Hueckel cutoff; this gives
    49 nm at 350 residues, close to the 50 nm that outreach_sim_1 uses.
    """
    return float(max(20, math.ceil(6 * 0.2 * n_residues**0.6 + 8)))


def run_plan(n_residues: int, preset: str) -> dict:
    frames = PRESETS[preset]["frames"]
    wfreq = PRESETS[preset].get("wfreq") or frame_interval(n_residues)
    return {
        "frames": frames,
        "wfreq": wfreq,
        "steps": frames * wfreq,
        "box": box_length(n_residues),
        "simulated_ns": frames * wfreq * TIMESTEP_PS / 1000,
    }


def bead_masses(seq: str, residues_csv: Path) -> np.ndarray:
    with open(residues_csv, newline="") as f:
        mw = {row["one"]: float(row["MW"]) for row in csv.DictReader(f)}
    masses = np.array([mw[aa] for aa in seq])
    masses[0] += N_TERM_MASS
    masses[-1] += C_TERM_MASS
    return masses


def unwrap_chain(xyz: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Rebuild one chain across periodic boundaries, bond by bond.

    xyz is (frames, beads, 3) and box is (frames, 3).
    """
    bonds = np.diff(xyz, axis=1)
    bonds -= box[:, None, :] * np.round(bonds / box[:, None, :])
    return np.concatenate([xyz[:, :1], xyz[:, :1] + np.cumsum(bonds, axis=1)], axis=1)


def load_chain(dcd: Path, top: Path) -> np.ndarray:
    """Whole-chain coordinates in nm, shape (frames, beads, 3)."""
    import mdtraj as md  # slow to import, and the UI server rarely needs it

    traj = md.load(str(dcd), top=str(top))
    return unwrap_chain(traj.xyz.astype(float), traj.unitcell_lengths.astype(float))


def radius_of_gyration(xyz: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Mass-weighted Rg per frame, the same definition as calvados.analysis.calc_rg."""
    com = np.average(xyz, axis=1, weights=masses)[:, None, :]
    rg_sq = np.sum(masses[None, :, None] * (xyz - com) ** 2, axis=(1, 2)) / masses.sum()
    return np.sqrt(rg_sq)


def cumulative_mean(rg: np.ndarray, start: int = 0) -> np.ndarray:
    """Mean Rg over frames `start` through each frame; NaN before `start`."""
    out = np.full(len(rg), np.nan)
    out[start:] = np.cumsum(rg[start:]) / np.arange(1, len(rg) - start + 1)
    return out


def summarize(rg: np.ndarray) -> dict:
    """Statistics over the frames left after dropping equilibration.

    The SEM comes from block means because successive frames are correlated.
    `settled` is False when the two halves of the run disagree beyond their noise,
    which means the mean isn't trustworthy yet.
    """
    n_discard = max(1, int(len(rg) * EQUILIBRATION_FRACTION))
    used = rg[n_discard:]
    blocks = np.array_split(used, N_BLOCKS)
    block_means = np.array([b.mean() for b in blocks])
    sem = float(block_means.std(ddof=1) / math.sqrt(N_BLOCKS))
    half = len(used) // 2
    first, second = float(used[:half].mean()), float(used[half:].mean())
    return {
        "rg_mean": float(used.mean()),
        "rg_sem": sem,
        "rg_std": float(used.std(ddof=1)),
        "frames_used": int(len(used)),
        "frames_discarded": n_discard,
        "first_half": first,
        "second_half": second,
        "settled": bool(abs(first - second) <= 4 * sem),
    }
