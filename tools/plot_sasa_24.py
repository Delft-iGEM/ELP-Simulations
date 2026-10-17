"""Domain exposure across the 24 free-chain RGD arrangements, one figure per motif.

    .venv/bin/python tools/plot_sasa_24.py

Reads figures/free-rgdarr/sasa-24-arrangements.csv (free-rgdarr* runs: one free
chain, 200 ns, 293.15 K, no TEA) and writes, next to it,

  sasa-24-arrangements-rgd.png               the RGD tripeptide, 3 residues x 4 copies
  sasa-24-arrangements-23mer.png             RGD + linker, the 23-mer x 4 copies
  sasa-24-arrangements-{rgd,23mer}-schematic.png   the same, untitled, with a sketch of each chain
                                             (yellow bar, red dot per RGD) beside its row

Each figure is sorted by its own SASA. All are slide versions: transparent,
for a light background, without subtitle, caption or % column. The original, and the winner and loser
of the RGD + linker ranking, are set in bold and ringed in both. Colour and
marker shape: how many of the two chain ends (N and C terminus) carry an RGD block. Dots are the production-window
mean, bars +-1 autocorrelation-corrected SEM, the dashed line the mean over the
24; the percentage is of the same beads' SASA with nothing around them.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

CSV = Path("figures/free-rgdarr/sasa-24-arrangements.csv")
INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
# copies at a chain end - the colours of the original plot
END_COLOUR = {0: "#2a6fd1", 1: "#e8642c", 2: "#1aaa7a"}
END_LABEL = {0: "RGD at neither chain end", 1: "RGD at one chain end",
             2: "RGD at both chain ends"}
# shape as well as colour, so the classes stay apart in greyscale and for CVD
END_MARKER = {0: "o", 1: "s", 2: "^"}
MARKER_SCALE = {"o": 1.0, "s": 0.95, "^": 1.2}     # a triangle reads small at the same ms
ORIGINAL = "01-orig"
# chain schematic: every rgdarr construct is 352 aa
CHAIN_LEN = 352
CHAIN_FILL, CHAIN_EDGE, RGD_DOT = "#f5c518", "#b8920c", "#d7263d"

# the last field marks the slide version: transparent, for a light background,
# no subtitle, caption or % column
PANELS = [
    ("rgd", "rgd_nm2", "rgd_sem_nm2", "rgd_exposed_frac",
     "RGD alone", "3 residues × 4 copies", True),
    ("23mer", "motif_23mer_nm2", "motif_23mer_sem_nm2", "motif_23mer_exposed_frac",
     "RGD + linker", "LDASTVYAVTGRGDSPASSAASA, 23 residues × 4 copies", True),
]
SCHEMATIC = {"rgd", "23mer"}                  # panels that also get a version with chain sketches


def draw_schematic(sx, rows) -> None:
    """Each chain as a yellow bar, N terminus left, with a red dot per RGD."""
    for side in ("top", "right", "left", "bottom"):
        sx.spines[side].set_visible(False)
    sx.tick_params(left=False, bottom=False, labelbottom=False)
    sx.set_facecolor(SURFACE)
    sx.set_xlim(-0.05 * CHAIN_LEN, 1.05 * CHAIN_LEN)
    for i, r in enumerate(rows):
        sx.plot([0, CHAIN_LEN], [i, i], color=CHAIN_EDGE, lw=9.5, solid_capstyle="round",
                zorder=1)
        sx.plot([0, CHAIN_LEN], [i, i], color=CHAIN_FILL, lw=8, solid_capstyle="round",
                zorder=2)
        # copy_starts are 1-based 23-mer starts; R G D are residues 12-14 of the
        # 23-mer, so the G (the dot's centre) is residue start + 12
        rgd = [int(c) + 12 for c in r["copy_starts"].split()]
        sx.plot(rgd, [i] * len(rgd), "o", color=RGD_DOT, ms=6, mec=INK, mew=0.6, zorder=3)
    top = len(rows) - 0.45
    for x, lab in ((0, "N"), (CHAIN_LEN, "C")):
        sx.text(x, top, lab, ha="center", va="bottom", fontsize=9, color=INK,
                fontweight="bold")
    sx.set_xlabel("RGD positions", fontsize=10, color=INK)


def draw(plt, rows, panel, role, schematic: bool, out: Path | None = None,
         xlabel: str = "SASA (nm²)") -> None:
    tag, key, sem, frac, name, sub, slide = panel
    muted = INK if slide else MUTED             # all text black on the slide version
    rows = sorted(rows, key=lambda r: r[key])           # own ranking, bottom to top
    y = np.arange(len(rows))
    v = np.array([r[key] for r in rows])
    e = np.array([r[sem] for r in rows])
    if schematic:
        fig, (ax, sx) = plt.subplots(1, 2, figsize=(9.8, 9.6), sharey=True,
                                     gridspec_kw={"width_ratios": (4.2, 1.3), "wspace": 0.05})
    else:
        fig, ax = plt.subplots(figsize=(8.2, 9.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(INK if slide else GRID)
    ax.tick_params(colors=muted, labelsize=9)
    ax.tick_params(axis="y", length=0)
    ax.grid(color=GRID, linewidth=0.7, axis="x")
    ax.set_axisbelow(True)

    for i, r in enumerate(rows):
        end = int(r["copies_at_chain_end"])
        c = END_COLOUR[end]
        bold = r["run"] in role
        if bold:
            for a in (ax, sx) if schematic else (ax,):
                a.axhspan(i - 0.42, i + 0.42, color=GRID, alpha=0.45, lw=0, zorder=0)
        ax.errorbar(v[i], i, xerr=e[i], fmt=END_MARKER[end], color=c,
                    ms=(8.5 if bold else 6.5) * MARKER_SCALE[END_MARKER[end]],
                    elinewidth=2.2 if bold else 1.5, capsize=3, capthick=2.2 if bold else 1.5,
                    mec=INK if bold else SURFACE, mew=1.6 if bold else 1.0, zorder=3)
    ax.axvline(v.mean(), color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['run']}  ({role[r['run']]})" if r["run"] in role
                        else r["run"] for r in rows])
    for lab, r in zip(ax.get_yticklabels(), rows):
        if r["run"] in role:
            lab.set_fontweight("bold")
            lab.set_color(INK)
        else:
            lab.set_color(muted)
    ax.set_ylim(-0.7, len(rows) - 0.3)

    # exposed fraction, at the right edge
    lo, hi = (v - e).min(), (v + e).max()
    pad = (hi - lo) * 0.08
    ax.set_xlim(lo - pad, hi + pad * (1 if slide else 3.2))
    for i, r in enumerate(rows):
        if slide:
            break
        ax.text(1.0, i, f"{100 * r[frac]:.1f}%", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8.3,
                color=INK if r["run"] in role else MUTED,
                fontweight="bold" if r["run"] in role else "normal")
    ax.set_xlabel(xlabel, fontsize=10, color=INK)
    if schematic:
        draw_schematic(sx, rows)

    handles = [plt.Line2D([], [], marker=END_MARKER[n], ls="", color=c, mec=SURFACE,
                          ms=7.5 * MARKER_SCALE[END_MARKER[n]], label=END_LABEL[n])
               for n, c in END_COLOUR.items()]
    if slide:                                    # no caption to explain the dashed line
        handles.append(plt.Line2D([], [], color=MUTED, lw=1, ls=(0, (4, 3)),
                                  label=f"mean of the {len(rows)}"))
    if schematic:
        handles.append(plt.Line2D([], [], marker="o", ls="", color=RGD_DOT, mec=INK,
                                  mew=0.6, ms=6, label="RGD position (yellow bar = chain)"))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, 1.01),
              ncol=2 if slide else 3, frameon=False, fontsize=8.6, labelcolor=INK,
              columnspacing=1.4)
    if not schematic:                            # the schematic version goes without a title
        fig.suptitle(f"{name}: exposure across all 24 RGD arrangements", x=0.02, y=0.985,
                 ha="left", fontsize=13, color=INK, fontweight="bold")
    spread = v.max() - v.min()
    out = out or CSV.parent / f"sasa-24-arrangements-{tag}{'-schematic' if schematic else ''}.png"
    if slide:
        ax.get_legend().set_bbox_to_anchor((0.42, 1.0))
        if schematic:       # tight_layout cannot handle the shared-y pair; set margins
            fig.subplots_adjust(left=0.2, right=0.985, bottom=0.07, top=0.9)
        else:
            fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(out, dpi=200, transparent=True)
        plt.close(fig)
        print(f"wrote {out}")
        return
    fig.text(0.02, 0.952, sub, fontsize=9.5, color=MUTED, ha="left")
    fig.text(0.02, 0.012,
             "One free chain per run, 200 ns, every 8th frame from 50 ns on. Dot = mean, bar = "
             "±1 SEM corrected for autocorrelation;\ndashed line = mean over the 24; % = share "
             "of the same beads' SASA with nothing around them. Rows sorted by this figure's SASA;\n"
             f"bold: the original and the winner and loser by RGD + linker SASA. Spread {spread:.2f} nm², typical error "
             f"{np.median(e):.2f} nm².", fontsize=8.2, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.075, 1, 0.905))
    fig.subplots_adjust(top=0.875)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [{k: (v if k in ("run", "copy_starts") else float(v)) for k, v in r.items()}
            for r in csv.DictReader(open(CSV))]
    by_23mer = sorted(rows, key=lambda r: r["motif_23mer_nm2"])
    winner, loser = by_23mer[-1]["run"], by_23mer[0]["run"]
    role = {ORIGINAL: "original", winner: "winner", loser: "loser"}
    print(f"winner {winner}, loser {loser}")
    for panel in PANELS:
        draw(plt, rows, panel, role, schematic=False)
        if panel[0] in SCHEMATIC:
            draw(plt, rows, panel, role, schematic=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
