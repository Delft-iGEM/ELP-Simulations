"""Submodule 4: lysine contacts (10 and 12 A), radius of gyration and SASA, every run.

    .venv/bin/python tools/sub4_analysis.py               # analyse every sub4 run found + plot
    .venv/bin/python tools/sub4_analysis.py --plot-only   # re-plot from the CSV
    .venv/bin/python tools/sub4_analysis.py --runs 'sub4-extra2-*'   # a subset (glob)

Every simulations/sub4-* folder with a finished trajectory is analysed; the
others are listed and skipped, so the same command works on any checkout that
holds a different set of runs (results for runs not found here are kept from
the existing CSV, so two people's halves merge into one table).

Window: 100 ns to the end, the same start as the earlier sub4 contact table.

Lysine contacts (the notebook's crosslink-contacts definition, analyze.ipynb)
    beads    every K bead of every chain
    pairs    all K-K pairs; same-chain pairs fewer than 4 residues apart are
             dropped, and so is any pair whose two lysines are BOTH bonded to
             the surface (their separation is fixed by construction)
    hit      minimum-image distance below the cutoff (10 A and 12 A)
    reported per cutoff, intra-chain and inter-chain: hits per frame,
             P(pair in contact) = hits / (pairs x frames), pairs ever in
             contact, closest approach, and the inter-chain share of all hits.
    In crosslinked runs a lysine pair that has bonded stays within ~6 A for
    the rest of the run, so it counts as a contact in every later frame; the
    crosslinks themselves are counted separately (xl_*).

Radius of gyration
    per chain, chain unwrapped bond by bond, every bead weighted equally;
    mean over chains, then over frames.

SASA (tools/sasa.py: coarse-grained Shrake-Rupley, bead radius sigma/2, probe
0.14 nm, minimum image; every bead of all chains can occlude)
    lys_sasa   per lysine, every 10th frame; all lysines and only the ones
               not bonded to the surface
    chain_sasa per chain and per residue, every 80th frame (it is the costly one)

Errors are the SEM from 20 block averages over frames (per run). Crosslinked
runs come in 3 replicates per design; the figures show every replicate and
their mean.

Writes, in figures/sub4/:
    sub4-analysis.csv          one row per run, ';'-separated with decimal
                               commas (opens directly in Dutch Excel)
    sub4-analysis-columns.txt  what every column means
    sub4-{xl,noxl}-*.png       figures; crosslinked and non-crosslinked runs are
                               never in the same figure
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

START_NS = 100.0
CUTOFFS_A = (10.0, 12.0)
MIN_SEQ_SEP = 4
LYS_STRIDE = 10
CHAIN_STRIDE = 80
OUT_DIR = Path("figures/sub4")
OUT_CSV = OUT_DIR / "sub4-analysis.csv"
NAME_RE = re.compile(r"sub4-(?P<family>[a-z0-9]+(?:-a(?P<design>\d+)|-base-case)?)"
                     r"(?:-l(?P<L>\d+)-s(?P<S>\d+)-c(?P<C>[\dp]+)-i(?P<I>[\dp]+))?"
                     r"(?:-(?:rep(?P<rep>\d+)|(?P<idx>\d+)))?$")


def registry() -> dict[str, dict]:
    return {r["name"]: r for r in csv.DictReader(open("used-runs.csv", newline=""))
            if r["name"].startswith("sub4")}


def has_trajectory(name: str) -> bool:
    from tools.analysis_batch import trajectory_paths
    try:
        _, dcd = trajectory_paths(name)
    except (FileNotFoundError, ValueError):
        return False
    return Path(dcd).is_file()


def analyse(name: str) -> dict | None:
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.sasa import bead_radii, sasa_of_beads
    from tools.tea_analysis import block_average
    from tools.tea_g6_analysis import _whole

    reg = registry().get(name, {})
    runtime = Path("simulations") / name / "runtime"
    top, dcd = trajectory_paths(name)
    traj = md.load(str(dcd), top=str(top))
    meta = read_metadata(runtime)
    ns_per_frame = float(meta["time_per_frame_ps"]) / 1000.0
    first = max(0, int(np.ceil(START_NS / ns_per_frame)) - 1)
    n_total = traj.n_frames
    traj = traj[first:]
    seq = "".join((runtime / "molecules.fasta").read_text().splitlines()[1:]).strip()
    n_res = len(seq)
    n_chain = traj.n_atoms // n_res
    box = traj.unitcell_lengths[0]                      # nm

    # ---- which lysines are surface-bonded --------------------------------------
    pinned = set()
    pre_frac = float("nan")
    state = runtime / "surface_state.json"
    if state.is_file():
        st = json.loads(state.read_text())
        pre_frac = float(st.get("settings", {}).get("preattached_fraction") or "nan")
        sites = [int(b) for b in st.get("site_beads", [])]
        pinned = {sites[int(i)] for i in st.get("pinned", []) if int(i) < len(sites)}

    lys_res = [i for i, c in enumerate(seq) if c == "K"]
    lys = np.array([c * n_res + i for c in range(n_chain) for i in lys_res])
    chain_of = lys // n_res
    pos_of = lys % n_res
    is_pin = np.array([b in pinned for b in lys])
    iu, ju = np.triu_indices(len(lys), k=1)
    same = chain_of[iu] == chain_of[ju]
    both_pin = is_pin[iu] & is_pin[ju]
    keep = (~same | (np.abs(pos_of[iu] - pos_of[ju]) >= MIN_SEQ_SEP)) & ~both_pin
    iu, ju, same = iu[keep], ju[keep], same[keep]

    xyz_a = traj.xyz * 10.0                               # A
    box_a = box * 10.0
    hits = {c: np.zeros((traj.n_frames, 2), dtype=int) for c in CUTOFFS_A}
    ever = {c: np.zeros(len(iu), dtype=bool) for c in CUTOFFS_A}
    dmin = np.full(len(iu), np.inf)
    for f in range(traj.n_frames):
        d = xyz_a[f, lys[iu]] - xyz_a[f, lys[ju]]
        d -= box_a * np.round(d / box_a)
        r = np.linalg.norm(d, axis=1)
        np.minimum(dmin, r, out=dmin)
        for c in CUTOFFS_A:
            h = r < c
            hits[c][f] = (np.count_nonzero(h & same), np.count_nonzero(h & ~same))
            ever[c] |= h

    # ---- Rg -------------------------------------------------------------------
    rg = np.empty(traj.n_frames)
    for f in range(traj.n_frames):
        per = []
        for c in range(n_chain):
            w = _whole(traj.xyz[f, c * n_res:(c + 1) * n_res], box)
            per.append(np.sqrt(((w - w.mean(axis=0)) ** 2).sum(axis=1).mean()))
        rg[f] = np.mean(per)

    # ---- SASA ------------------------------------------------------------------
    radii = bead_radii(traj.topology, Path("residues_CALVADOS2.csv"))
    lys_all, lys_free = [], []
    free_mask = ~is_pin
    for f in range(0, traj.n_frames, LYS_STRIDE):
        s = sasa_of_beads(traj.xyz[f], radii, which=lys, box=box)
        lys_all.append(s.mean())
        lys_free.append(s[free_mask].mean() if free_mask.any() else np.nan)
    chain = [sasa_of_beads(traj.xyz[f], radii, box=box).sum() / n_chain
             for f in range(0, traj.n_frames, CHAIN_STRIDE)]

    # ---- crosslinks formed ------------------------------------------------------
    xl_intra = xl_inter = 0
    events = runtime / "crosslink_events.csv"
    if events.is_file():
        for e in csv.DictReader(open(events, newline="")):
            if e.get("kind") == "intra":
                xl_intra += 1
            elif e.get("kind") == "inter":
                xl_inter += 1

    m = NAME_RE.match(name)
    crosslinked = bool(reg.get("crosslink_distance_nm", "").strip())
    row = {
        "run": name,
        "group": "crosslinked" if crosslinked else "not_crosslinked",
        "family": m.group("family") if m else "",
        "design": int(m.group("design")) if m and m.group("design") else "",
        "replicate": int(m.group("rep")) if m and m.group("rep") else "",
        "lysine_spacing_S": int(m.group("S")) if m and m.group("S") else "",
        "pentads_L": int(m.group("L")) if m and m.group("L") else "",
        "residues_per_chain": n_res,
        "lysines_per_chain": len(lys_res),
        "chains": n_chain,
        "concentration_chains_nm2": float(meta.get("concentration_chains_per_nm2") or "nan"),
        "preattached_fraction": pre_frac,
        "ionic_strength_M": float(meta.get("ionic_strength_M") or "nan"),
        "crosslink_distance_nm": float(reg["crosslink_distance_nm"]) if crosslinked else float("nan"),
        "frames_total": n_total,
        "first_frame": first,
        "frames_analysed": traj.n_frames,
        "window_ns": f"{first * ns_per_frame:.0f}-{n_total * ns_per_frame:.0f}",
        "lysines_surface_bonded": int(is_pin.sum()),
        "pairs_intra": int(same.sum()),
        "pairs_inter": int((~same).sum()),
        "pairs_dropped_both_bonded": int(both_pin.sum()),
        "closest_intra_A": float(dmin[same].min()) if same.any() else float("nan"),
        "closest_inter_A": float(dmin[~same].min()) if (~same).any() else float("nan"),
    }
    for c in CUTOFFS_A:
        t = f"{c:g}A"
        hi, he = hits[c][:, 0], hits[c][:, 1]
        n_i, n_e = row["pairs_intra"], row["pairs_inter"]
        tot = hi.sum() + he.sum()
        row |= {
            f"hits_per_frame_intra_{t}": float(hi.mean()),
            f"hits_per_frame_inter_{t}": float(he.mean()),
            f"p_pair_intra_{t}": float(hi.sum() / (n_i * len(hi))) if n_i else float("nan"),
            f"p_pair_inter_{t}": float(he.sum() / (n_e * len(he))) if n_e else float("nan"),
            f"inter_share_of_hits_{t}": float(he.sum() / tot) if tot else float("nan"),
            f"pairs_ever_intra_{t}": int((ever[c] & same).sum()),
            f"pairs_ever_inter_{t}": int((ever[c] & ~same).sum()),
        }
        row[f"inter_share_sem_{t}"] = block_average(
            he / np.maximum(hi + he, 1))[1] if tot else float("nan")
    row["rg_nm"], row["rg_sem_nm"] = block_average(rg)
    row["lys_sasa_nm2"], row["lys_sasa_sem_nm2"] = block_average(np.array(lys_all))
    row["lys_free_sasa_nm2"], row["lys_free_sasa_sem_nm2"] = block_average(np.array(lys_free))
    row["chain_sasa_nm2"], row["chain_sasa_sem_nm2"] = block_average(np.array(chain))
    row["chain_sasa_per_residue_nm2"] = row["chain_sasa_nm2"] / n_res
    row["chain_sasa_per_residue_sem_nm2"] = row["chain_sasa_sem_nm2"] / n_res
    row["xl_formed"] = xl_intra + xl_inter
    row["xl_intra_loops"] = xl_intra
    row["xl_inter_bridges"] = xl_inter
    row["xl_bridge_fraction"] = xl_inter / (xl_intra + xl_inter) if xl_intra + xl_inter else float("nan")
    return row


# ---- CSV in Dutch-Excel format ---------------------------------------------------

def _fmt(v) -> str:
    if isinstance(v, float):
        return "" if np.isnan(v) else f"{v:.6g}".replace(".", ",")
    return str(v)


def write_csv(rows: list[dict], path: Path) -> None:
    fields = list(rows[0])
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(fields)
        for r in rows:
            w.writerow([_fmt(r.get(k, "")) for k in fields])


def read_csv(path: Path) -> list[dict]:
    text = ("run", "group", "family", "window_ns")
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig", newline=""), delimiter=";"):
        out.append({k: (v if k in text else float(v.replace(",", ".")) if v else float("nan"))
                    for k, v in r.items()})
    return out


COLUMNS = """\
sub4-analysis.csv - one row per submodule 4 run. ';'-separated, decimal comma.
Made by tools/sub4_analysis.py. Window: 100 ns to the end of the run.

run, group                 run name; crosslinked / not_crosslinked (physical crosslinking on/off)
family, design, replicate  from the name (anchor a10-a19 x rep1-3, extra, extra2, ...)
lysine_spacing_S           pentads between lysines (the -s in the name)
pentads_L                  the -l in the name
residues_per_chain, lysines_per_chain, chains
concentration_chains_nm2, preattached_fraction, ionic_strength_M, crosslink_distance_nm
frames_total, first_frame, frames_analysed, window_ns
lysines_surface_bonded     lysines pinned to the surface (all chains)
pairs_intra, pairs_inter   K-K pairs monitored within one chain / between chains
                           (same-chain pairs < 4 residues apart and pairs with both
                           lysines surface-bonded are left out)
pairs_dropped_both_bonded
closest_intra_A, closest_inter_A   closest approach over the window, Angstrom

per cutoff (suffix _10A and _12A), K-K bead distance below the cutoff = one hit:
hits_per_frame_intra/inter average number of pairs in contact per frame
p_pair_intra/inter         P(a given pair is in contact) = hits / (pairs x frames)
inter_share_of_hits        share of all hits that are between different chains
inter_share_sem            its 20-block SEM over frames
pairs_ever_intra/inter     pairs that came within the cutoff at least once

rg_nm, rg_sem_nm           radius of gyration per chain (unwrapped, beads unweighted),
                           mean over chains and frames; SEM from 20 blocks
lys_sasa_nm2               SASA per lysine bead (all lysines), every 10th frame
lys_free_sasa_nm2          the same, only lysines not bonded to the surface
chain_sasa_nm2             SASA of one whole chain, every 80th frame
chain_sasa_per_residue_nm2 the same divided by residues_per_chain
(each with a _sem_ column: 20-block SEM over frames)
  SASA: coarse-grained Shrake-Rupley (bead radius sigma/2, probe 0.14 nm); every
  bead of every chain can bury another, the grafting wall does not.

xl_formed, xl_intra_loops, xl_inter_bridges, xl_bridge_fraction
                           crosslinks formed during the run (crosslink_events.csv);
                           0 for runs without crosslinking. In crosslinked runs a
                           bonded lysine pair stays in contact afterwards, so it is
                           also counted in the contact columns.
"""


# ---- figures ---------------------------------------------------------------------

INK, MUTED, GRID, SURFACE = "#1a1a1a", "#6b6f76", "#dcdee2", "#fcfcfb"
GROUP_TITLE = {"not_crosslinked": "no crosslinking", "crosslinked": "physical crosslinking on"}


def plot(rows: list[dict], out_dir: Path) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.ticker import PercentFormatter

    cmap = LinearSegmentedColormap.from_list(
        "density", ["#a9c8f2", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"])
    conc = [r["concentration_chains_nm2"] for r in rows]
    norm = Normalize(min(conc), max(conc))
    written = []

    def figure(group, tag, title, ylabel, x_key, xlabel, series, logy=False, pct=False,
               caption=""):
        sub = [r for r in rows if r["group"] == group]
        if not sub:
            return
        fig, ax = plt.subplots(figsize=(7.6, 5.2))
        fig.patch.set_facecolor(SURFACE)
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        handles = []
        for key, sem, marker, label in series:
            designs = {}
            for r in sub:
                if np.isnan(r[key]):
                    continue
                c = cmap(norm(r["concentration_chains_nm2"]))
                small = group == "crosslinked"
                ax.errorbar(r[x_key], r[key], yerr=None if sem is None else r[sem],
                            fmt=marker, ms=5 if small else 7.5, color=c,
                            mfc=c if marker != "s" else SURFACE, mec=c if marker == "s" else SURFACE,
                            mew=1.5 if marker == "s" else 1, elinewidth=1, capsize=0,
                            alpha=0.75 if small else 1, zorder=3)
                if r["design"] != "" and not np.isnan(r["design"]):
                    designs.setdefault(int(r["design"]), []).append(r)
            for d, reps in designs.items():           # replicate mean, crosslinked runs
                x = reps[0][x_key]
                y = np.mean([q[key] for q in reps])
                c = cmap(norm(reps[0]["concentration_chains_nm2"]))
                ax.plot(x, y, "D", ms=8.5, color=c, mec=INK, mew=1.2, zorder=4)
                ax.annotate(f"a{d}", (x, y), xytext=(7, 0), textcoords="offset points",
                            va="center", fontsize=7.5, color=MUTED)
            if label:
                handles.append(plt.Line2D([], [], ls="", marker=marker, ms=7, color=MUTED,
                                          mfc=MUTED if marker != "s" else SURFACE, mew=1.5,
                                          label=label))
        if group == "crosslinked":
            handles.append(plt.Line2D([], [], ls="", marker="D", ms=8, color=MUTED, mec=INK,
                                      label="mean of 3 replicates"))
        if logy:
            ax.set_yscale("log")
        if pct:
            ax.yaxis.set_major_formatter(PercentFormatter(1.0))
            ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel(xlabel, fontsize=9.5, color=INK)
        ax.set_ylabel(ylabel, fontsize=9.5, color=INK)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=ax, pad=0.015)
        cb.set_label("Grafting density (chains/nm²)", fontsize=8.5, color=INK)
        cb.ax.tick_params(colors=MUTED, labelsize=8)
        cb.outline.set_edgecolor(GRID)
        if handles:
            ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13),
                      ncol=len(handles), frameon=False, fontsize=8.5, labelcolor=INK)
        n = len(sub)
        ax.set_title(f"{title}\nSubmodule 4, {GROUP_TITLE[group]} ({n} runs)", loc="left",
                     fontsize=10.5, color=INK, fontweight="bold")
        fig.text(0.01, 0.01, caption + f"\nWindow 100 ns to end. Numbers: figures/sub4/"
                 "sub4-analysis.csv.", fontsize=7.6, color=MUTED, linespacing=1.45)
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        out = out_dir / f"sub4-{'xl' if group == 'crosslinked' else 'noxl'}-{tag}.png"
        fig.savefig(out, dpi=200, facecolor=SURFACE)
        plt.close(fig)
        written.append(out)

    S_LABEL = "Lysine spacing S (pentads between lysines)"
    for group in GROUP_TITLE:
        for c in CUTOFFS_A:
            t = f"{c:g}A"
            figure(group, f"contacts-{t}-inter-share",
                   f"Share of K–K contacts that are between chains, {c:g} Å",
                   "Inter-chain share of contacts", "lysine_spacing_S", S_LABEL,
                   [(f"inter_share_of_hits_{t}", f"inter_share_sem_{t}", "o", None)], pct=True,
                   caption=f"Share of all K–K pairs within {c:g} Å (pair-frames) that join two "
                           "different chains; bar = 20-block SEM.")
            figure(group, f"contacts-{t}-p-pair",
                   f"Probability that a given K–K pair is in contact, {c:g} Å",
                   "P(pair within cutoff)", "lysine_spacing_S", S_LABEL,
                   [(f"p_pair_intra_{t}", None, "o", "same chain"),
                    (f"p_pair_inter_{t}", None, "s", "different chains")], logy=True,
                   caption="hits / (pairs × frames). Same-chain pairs < 4 residues apart and pairs "
                           "with both lysines on the surface are left out.")
        figure(group, "rg", "Radius of gyration per chain", "Rg (nm)", "residues_per_chain",
               "Residues per chain", [("rg_nm", "rg_sem_nm", "o", None)],
               caption="Chains unwrapped along their bonds, beads unweighted; mean over chains; "
                       "bar = 20-block SEM.")
        figure(group, "lys-sasa", "Lysine solvent exposure", "SASA per lysine (nm²)",
               "lysine_spacing_S", S_LABEL,
               [("lys_sasa_nm2", "lys_sasa_sem_nm2", "o", "all lysines"),
                ("lys_free_sasa_nm2", "lys_free_sasa_sem_nm2", "s", "not surface-bonded")],
               caption="Coarse-grained SASA (bead radius σ/2, probe 0.14 nm), every 10th frame; "
                       "neighbouring chains bury, the grafting wall does not.")
        figure(group, "chain-sasa", "Whole-chain solvent exposure, per residue",
               "SASA per residue (nm²)", "residues_per_chain", "Residues per chain",
               [("chain_sasa_per_residue_nm2", "chain_sasa_per_residue_sem_nm2", "o", None)],
               caption="SASA of one whole chain / its residues, every 80th frame; bar = "
                       "20-block SEM.")
    figure("crosslinked", "bridge-fraction", "Crosslinks that bridge two chains",
           "Bridges / all crosslinks formed", "lysine_spacing_S", S_LABEL,
           [("xl_bridge_fraction", None, "o", None)], pct=True,
           caption="From crosslink_events.csv: a bridge joins two chains, a loop joins one chain "
                   "to itself.")
    for p in written:
        print(f"wrote {p}")
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--runs", default="sub4-*", help="glob over simulations/ names")
    parser.add_argument("--jobs", type=int, default=25)
    args = parser.parse_args(argv)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        plot(read_csv(OUT_CSV), OUT_DIR)
        return 0

    names = sorted(n for n in registry() if fnmatch.fnmatch(n, args.runs))
    ready = [n for n in names if has_trajectory(n)]
    skipped = [n for n in names if n not in ready]
    print(f"{len(ready)} runs with a trajectory, {len(skipped)} without (skipped)")
    with Pool(args.jobs) as pool:
        new = pool.map(analyse, ready)
    old = {r["run"]: r for r in read_csv(OUT_CSV)} if OUT_CSV.is_file() else {}
    rows = {**old, **{r["run"]: r for r in new if r}}
    rows = [rows[k] for k in sorted(rows)]
    write_csv(rows, OUT_CSV)
    (OUT_DIR / "sub4-analysis-columns.txt").write_text(COLUMNS)
    print(f"wrote {OUT_CSV} ({len(rows)} runs)")
    plot(rows, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
