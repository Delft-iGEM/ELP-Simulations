"""Step 8 - pick the runs for the videos and write the contacts behind their metrics.

Renders nothing and only reads the trajectories. Writes to analysis/export/:
  video_runs.md                         one entry per run: why it was chosen, absolute paths
                                        (DCD, top.pdb, surface_state.json, settings), the frame
                                        where production starts, the total number of frames
  video_runs/<Name>_contacts.csv        every frozen-definition contact, every frame
  video_runs/<Name>_pins.csv            the surface-pinned lysines
  video_runs/<Name>_top.pdb, <Name>_surface_state.json   copies, so the four files sit together

Contacts use the frozen definitions (settings.yaml): lysine CA beads, 12 A, minimum image in
x, y and z, free lysines only (neither one pinned), intra-chain pairs at least 4 residues apart.
Each file is checked frame by frame against the extracted time series the metrics come from.

Selection rules (deterministic; ties go to the lower replicate number, then to the name):
  1 base case, replicate 1
  2 lowest grafting density among finished runs
  3 highest grafting density among finished runs
  4 longest chain and widest spacing: largest coded L + coded S (frozen linear coding)
  5 high vs low bridge-to-loop ratio at similar density: among pairs of runs from different
    design points whose densities differ by at most 5 %, the pair with the largest difference
    in log(bridge_loop_ratio_12A)
  6 wet-lab constructs (Monoblock, Triblock-32, Triblock-64, Tetrablock): only if they are in
    runs_metrics.csv - they are not SubM4 runs, so this is reported, never invented

    analysis/scripts/run.sh analysis/scripts/08_video_runs.py [--workers 4]
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, ROOT, RANGES, SIMS, TIMESERIES  # noqa: E402

S = yaml.safe_load((ANALYSIS / "settings.yaml").read_text())
CUT = float(S["response"]["cutoff_A"])
MIN_SEP = int(S["response"]["min_seq_separation_res"])
T0 = int(S["equilibration"]["discard_frames"])
OUT = ANALYSIS / "export"
VID = OUT / "video_runs"
CONSTRUCTS = ["monoblock", "triblock-32", "triblock-64", "tetrablock"]
SIMILAR_C = 0.05


def select(m: pd.DataFrame) -> list[dict]:
    m = m.sort_values(["replicate", "Name"], na_position="last")
    picks = []

    def add(row, why, slot):
        picks.append({"slot": slot, "Name": row.Name, "slug": row.slug, "why": why})

    b = m[(m.anchor_id == "BASE") & (m.replicate == 1)]
    if len(b) != 1:
        raise RuntimeError("base case replicate 1 not found exactly once")
    add(b.iloc[0], "Base case (BASE anchor, replicate 1).", "1")

    for slot, fn, word in [("2", "idxmin", "lowest"), ("3", "idxmax", "highest")]:
        v = getattr(m.concentration, fn[3:])()
        tied = m[np.isclose(m.concentration, v, rtol=1e-9)]
        r = tied.iloc[0]
        why = f"The {word} grafting density of all finished runs (C = {v:.4f} chains/nm²)."
        if len(tied) > 1:
            why += (f" Same simulated density as {', '.join(tied.Name.iloc[1:])}; ties go to the lower "
                    "replicate number, then the name.")
        if bool(r.get("qc_drift_primary", False)):
            why += (f" Note: flagged for drift on the primary response (|z| = {abs(r.inter_free_12A_drift_z):.1f};"
                    " kept in the analysis, README §5).")
        add(r, why, slot)

    def cod(col, p):
        lo, hi = RANGES[p]
        return (m[col] - lo) / (hi - lo) * 2 - 1

    score = cod("chain_length", "chain_length") + cod("spacing_pct", "spacing_pct")
    r = m.loc[score.idxmax()]
    longest = m.loc[m.chain_length.idxmax()]
    widest = m.loc[m.spacing_pct.idxmax()]
    why = (f"Longest chain combined with widest lysine spacing (largest coded L + coded S): L = {r.chain_length:.0f} "
           f"repeats, spacing {r.spacing_pct:.1f} % ({r.spacing_repeats:.0f} repeats).")
    if longest.slug != r.slug or widest.slug != r.slug:
        why += (f" No run has both maxima: the longest chain is {longest.Name} (L = {longest.chain_length:.0f}, "
                f"spacing {longest.spacing_pct:.1f} %), the widest spacing is {widest.Name} "
                f"({widest.spacing_pct:.1f} %).")
    add(r, why, "4")

    a = m.reset_index(drop=True)
    lc, lb = np.log(a.concentration.to_numpy()), np.log(a.bridge_loop_ratio_12A.to_numpy())
    best = None
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a.design_point[i] == a.design_point[j] or abs(lc[i] - lc[j]) > np.log(1 + SIMILAR_C):
                continue
            d = abs(lb[i] - lb[j])
            if best is None or d > best[0]:
                best = (d, i, j)
    _, i, j = best
    hi_, lo_ = (a.loc[i], a.loc[j]) if lb[i] > lb[j] else (a.loc[j], a.loc[i])
    for r, word, other in [(hi_, "high", lo_), (lo_, "low", hi_)]:
        add(r, f"{word.capitalize()} bridge-to-loop ratio ({r.bridge_loop_ratio_12A:.2f}) at similar density to "
               f"{other.Name} (ratio {other.bridge_loop_ratio_12A:.2f}); densities {r.concentration:.4f} vs "
               f"{other.concentration:.4f}. Chosen as the pair with the largest ratio difference among runs "
               f"whose densities differ by ≤ {SIMILAR_C:.0%}.", "5" + word[0])
    return picks


def constructs_report(m: pd.DataFrame) -> str:
    in_metrics = [c for c in CONSTRUCTS if m.Name.str.lower().str.startswith(c).any()]
    folders = {c: sorted(p.name for p in SIMS.glob(f"{c}*") if p.is_dir()) for c in CONSTRUCTS}
    readable = {c: [f for f in fs if (SIMS / f / "runtime" / "top.pdb").is_file()] for c, fs in folders.items()}
    lines = []
    if in_metrics:
        lines.append(f"Construct runs found in runs_metrics.csv: {in_metrics} (not expected - check).")
    else:
        lines.append("**None selected.** No run of Monoblock, Triblock-32, Triblock-64 or Tetrablock is in "
                     "`runs_metrics.csv`: every SubM4 run is a (VPGIG)ₙ chain with four VPGKG.")
    tot = sum(len(v) for v in folders.values())
    nread = sum(len(v) for v in readable.values())
    lines.append(f"Folders with these construct names do exist from earlier submodules ({tot} under `simulations/`: "
                 + ", ".join(f"{c} {len(folders[c])}" for c in CONSTRUCTS)
                 + f"), but {nread} of them have a readable `runtime/` (trajectory): they hold only `prepare.py` "
                 "and `analysis/` summaries. They were run under other protocols (box, length, burn-in) and are "
                 "not part of this dataset, so nothing was chosen for them.")
    return "\n\n".join(lines)


def contacts(slug: str) -> dict:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import self_distance_array
    warnings.filterwarnings("ignore")
    rt = SIMS / slug / "runtime"
    u = mda.Universe(str(rt / "top.pdb"), str(rt / f"{slug}.dcd"))
    lys = u.select_atoms("name CA and resname LYS")
    seg = lys.segindices
    first = np.array([s.residues[0].resindex for s in u.segments])
    within = lys.resindices - first[seg]
    chain = np.asarray(lys.chainIDs) if hasattr(lys, "chainIDs") else np.asarray(lys.segids)
    resid = np.asarray(lys.resids)
    if not np.array_equal(resid, within + 1):
        raise RuntimeError(f"{slug}: PDB resid is not the 1-based position in the chain")
    st = json.loads((rt / "surface_state.json").read_text())
    sites = [int(x) for x in st["site_beads"]]
    pinned = np.isin(lys.indices, [sites[int(i)] for i in st["pinned"]])
    iu, ju = np.triu_indices(lys.n_atoms, k=1)
    same = seg[iu] == seg[ju]
    keep = (~same | (np.abs(within[iu] - within[ju]) >= MIN_SEP)) & ~(pinned[iu] | pinned[ju])
    iu, ju, same = iu[keep], ju[keep], same[keep]

    rows, n_inter, n_intra = [], [], []
    for f, ts in enumerate(u.trajectory):
        hit = self_distance_array(lys.positions, box=ts.dimensions)[keep] < CUT
        h = np.flatnonzero(hit)
        n_inter.append(int((~same[h]).sum()))
        n_intra.append(int(same[h].sum()))
        if len(h):
            rows.append(np.c_[np.full(len(h), f), same[h], iu[h], ju[h]])
    a = np.vstack(rows) if rows else np.zeros((0, 4), int)
    k = a[:, 2:].astype(int)
    out = pd.DataFrame({"frame": a[:, 0].astype(int), "type": np.where(a[:, 1] == 1, "intra", "inter"),
                        "chain_A": chain[k[:, 0]], "resid_A": resid[k[:, 0]],
                        "chain_B": chain[k[:, 1]], "resid_B": resid[k[:, 1]],
                        "index_A": lys.indices[k[:, 0]], "index_B": lys.indices[k[:, 1]]})
    k_index = pd.Series(within).groupby(seg).rank(method="first").astype(int).to_numpy()
    pins = pd.DataFrame({"chain": chain, "resid": resid, "k_index": k_index, "index": lys.indices})[pinned]

    # check against the time series the metrics were computed from
    z = np.load(TIMESERIES / f"{slug}.npz")
    kB = list(z["defs"]).index("B_free_lysines_only")
    c12 = int(np.flatnonzero(z["cutoffs_A"] == CUT)[0])
    ok = (np.array_equal(z["counts"][:, kB, c12, 1], n_inter) and np.array_equal(z["counts"][:, kB, c12, 0], n_intra))
    return {"slug": slug, "contacts": out, "pins": pins, "n_frames": len(u.trajectory),
            "frames_match_npz": ok, "inter_prod_mean": float(np.mean(n_inter[T0:])),
            "intra_prod_mean": float(np.mean(n_intra[T0:])), "n_lys": lys.n_atoms}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    m = pd.read_csv(ANALYSIS / "runs_metrics.csv")
    m = m[(m.data_source == "trajectory") & ~m.excluded]
    picks = select(m)
    VID.mkdir(parents=True, exist_ok=True)
    slugs = list(dict.fromkeys(p["slug"] for p in picks))
    with ProcessPoolExecutor(args.workers) as ex:
        res = {r["slug"]: r for r in ex.map(contacts, slugs)}

    mi = m.set_index("slug")
    cook = pd.read_csv(ANALYSIS / "stats" / "run_diagnostics.csv").set_index("slug")
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    md = ["# SubM4 runs for the videos", "",
          f"Written by `analysis/scripts/08_video_runs.py` (repository HEAD `{commit}`). Nothing was rendered.",
          "",
          f"Frozen definitions used for the contact files: lysine CA beads, **{CUT:g} Å**, minimum image in x, y "
          f"and z, **free lysines only** (neither lysine surface-pinned), intra-chain pairs ≥ {MIN_SEP} residues apart. "
          f"Production starts at frame **{T0}** (100 ns burn-in); the metrics average frames {T0} onwards.",
          "",
          "Files per run in `video_runs/`: `<Name>_contacts.csv` (columns `frame, type, chain_A, resid_A, chain_B, "
          "resid_B, index_A, index_B`; one row per frame and contact, all frames including burn-in; frame is the "
          "0-based DCD frame index; resid is the residue number in `top.pdb`, 1-based within the chain; index is "
          "the 0-based atom index = VMD `index`), `<Name>_pins.csv` (`chain, resid, k_index, index`; K1-K4 counted "
          "from the N-terminus), and copies `<Name>_top.pdb`, `<Name>_surface_state.json`. The DCDs are not copied.",
          "",
          "Units: DCD and PDB in Å; z is the surface normal (wall at z = 19 Å). Chains are written whole.",
          ""]
    for p in picks:
        r, row = res[p["slug"]], mi.loc[p["slug"]]
        rt = SIMS / p["slug"] / "runtime"
        name = p["Name"]
        r["contacts"].to_csv(VID / f"{name}_contacts.csv", index=False)
        r["pins"].to_csv(VID / f"{name}_pins.csv", index=False)
        shutil.copy(rt / "top.pdb", VID / f"{name}_top.pdb")
        shutil.copy(rt / "surface_state.json", VID / f"{name}_surface_state.json")
        agree = (r["frames_match_npz"] and np.isclose(r["inter_prod_mean"], row.inter_free_12A)
                 and np.isclose(r["intra_prod_mean"], row.intra_free_12A))
        if not agree:
            raise RuntimeError(f"{name}: contact file does not reproduce the metrics")
        why = p["why"]
        if p["slug"] in cook.index and cook.loc[p["slug"], "qc_cooks_gt_4_over_n"]:
            why += (f" Note: influential in the frozen model (Cook's D = {cook.loc[p['slug'], 'cooks_d']:.2f}, "
                    f"studentised residual {cook.loc[p['slug'], 'resid_studentized']:.1f}; kept, README §5).")
        md += [f"## {p['slot']}. `{name}`", "", why, "",
               f"C = {row.concentration:.4f} chains/nm², L = {row.chain_length:.0f} repeats, spacing "
               f"{row.spacing_pct:.1f} % ({row.spacing_repeats:.0f} repeats), I = {row.ionic_strength:.3f} M. "
               f"Inter-chain free contacts {row.inter_free_12A:.3f}/frame, intra {row.intra_free_12A:.3f}/frame, "
               f"bridge-to-loop ratio {row.bridge_loop_ratio_12A:.2f}; {len(r['pins'])} of {r['n_lys']} lysines pinned.",
               "",
               "| | |", "|---|---|",
               f"| DCD | `{rt / (p['slug'] + '.dcd')}` |",
               f"| (resolves to) | `{(rt / (p['slug'] + '.dcd')).resolve()}` |",
               f"| top.pdb | `{rt / 'top.pdb'}` |",
               f"| surface_state.json | `{rt / 'surface_state.json'}` |",
               f"| settings | `{rt / 'config.yaml'}` (CALVADOS/OpenMM, nm); surface: `{rt / 'surface.yaml'}`; "
               f"planner: `{SIMS / p['slug'] / 'prepare.py'}` |",
               f"| production starts at frame | {T0} |",
               f"| total frames | {r['n_frames']} |",
               f"| contacts file check | every frame matches the extracted counts; production means "
               f"{r['inter_prod_mean']:.4f} inter / {r['intra_prod_mean']:.4f} intra = `runs_metrics.csv` |",
               ""]
    md += ["## 6. Wet-lab constructs", "", constructs_report(m), ""]
    (OUT / "video_runs.md").write_text("\n".join(md))
    for p in picks:
        print(p["slot"], p["Name"], "ok" if res[p["slug"]]["frames_match_npz"] else "MISMATCH")


if __name__ == "__main__":
    main()
