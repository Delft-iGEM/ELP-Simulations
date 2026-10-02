"""The three designed sequences for the 270/350 K set, built from the original.

The original RGD construct (357 aa, the one the 0r0g* crosslinking sweep used)
decomposes exactly into

    4 x RGD block (25 aa)  +  51 pentads (45 VPGIG, 4 VPGKG, 2 VPGMG)  +  "VP"

The first RGD block carries an LQ cloning scar where the others carry VP; the
trailing "VP" is a partial pentad at the C terminus. Both are kept verbatim so
the designed variants differ from the original only where they are meant to.

`all_terminal` is a pure permutation of the original - same 357 residues, same
composition, checked by `verify()`. `vpgig_only` is not: it is the bare ELP with
every functional group removed, at the same length.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

ORIGINAL_SOURCE = "0r0g0d1-2-4-7-focused"   # one of many runs sharing the sequence
BLOCK_LEN = 25
TAIL = "VP"


def original(root: Path | None = None) -> str:
    """The 357-aa original, read from the run that defines it."""
    base = Path(root) if root else Path(".")
    text = (base / "simulations" / ORIGINAL_SOURCE / "prepare.py").read_text()
    return re.search(r'"([ACDEFGHIKLMNPQRSTVWYZ]{50,})"', text).group(1)


def decompose(seq: str) -> tuple[list[str], list[str], str]:
    """(RGD blocks, pentads, tail) - the inverse of the constructors below."""
    starts = [m.start() - 13 for m in re.finditer("RGD", seq)]
    blocks = [seq[a:a + BLOCK_LEN] for a in starts]
    spans, prev = [], 0
    for a in starts:
        if a > prev:
            spans.append(seq[prev:a])
        prev = a + BLOCK_LEN
    spans.append(seq[prev:])
    pentads, tail = [], ""
    for span in spans:
        for i in range(0, len(span), 5):
            piece = span[i:i + 5]
            if len(piece) == 5:
                pentads.append(piece)
            else:
                tail = piece          # the C-terminal partial pentad
    return blocks, pentads, tail


def all_terminal(root: Path | None = None) -> str:
    """Every RGD block and every lysine pushed to the two chain ends.

    Reading inward from each terminus: two RGD blocks, then two VPGKG. The
    middle is plain VPGIG apart from the two VPGMG the original carries, which
    are left near the centre - methionine is not the variable here, and dropping
    it would change the composition and so the comparison.
    """
    blocks, pentads, tail = decompose(original(root))
    lys = [p for p in pentads if p == "VPGKG"]
    met = [p for p in pentads if p == "VPGMG"]
    ile = [p for p in pentads if p == "VPGIG"]
    if len(blocks) != 4 or len(lys) != 4:
        raise ValueError(f"expected 4 RGD blocks and 4 VPGKG, got {len(blocks)} and {len(lys)}")

    # VPGMG at 1/3 and 2/3 of the plain run, so the middle stays homogeneous.
    third = len(ile) // 3
    middle = ile[:third] + met[:1] + ile[third:2 * third] + met[1:] + ile[2 * third:]
    return "".join(blocks[:2] + lys[:2] + middle + lys[2:] + blocks[2:]) + tail


def vpgig_only(root: Path | None = None) -> str:
    """(VPGIG)n + "VP" at exactly the original's length - the bare ELP baseline."""
    target = len(original(root))
    n, rest = divmod(target - len(TAIL), 5)
    if rest:
        raise ValueError(f"{target} aa is not n*5 + len({TAIL!r})")
    return "VPGIG" * n + TAIL


def verify(root: Path | None = None) -> dict:
    """all_terminal must be a permutation of the original; vpgig_only matches length."""
    orig = original(root)
    term = all_terminal(root)
    bare = vpgig_only(root)
    out = {
        "original_len": len(orig),
        "all_terminal_len": len(term),
        "vpgig_only_len": len(bare),
        "composition_identical": Counter(orig) == Counter(term),
        "sequences_differ": orig != term,
        "original_K": [i for i, c in enumerate(orig) if c == "K"],
        "all_terminal_K": [i for i, c in enumerate(term) if c == "K"],
        "original_RGD": [m.start() for m in re.finditer("RGD", orig)],
        "all_terminal_RGD": [m.start() for m in re.finditer("RGD", term)],
        "vpgig_only_has_no_K_or_RGD": "K" not in bare and "RGD" not in bare,
    }
    if not out["composition_identical"]:
        raise ValueError(f"all_terminal is not a permutation: "
                         f"{Counter(orig) - Counter(term)} / {Counter(term) - Counter(orig)}")
    if len(term) != len(orig) or len(bare) != len(orig):
        raise ValueError("lengths differ")
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(verify(), indent=2))
