"""Step 2c - per-run metrics over the production window, plus convergence QC.

Reads the frozen definitions in analysis/settings.yaml (refuses to run while the burn-in
is not frozen) and the extracted time series. Writes
  analysis/runs_metrics.csv          one row per run
  analysis/per_lysine.csv            one row per lysine per run (crosslink-degree distribution)
  analysis/timeseries_primary.csv.gz per-frame primary response for every run (long format)
  analysis/reactive_anchors.csv      teammate's reactive anchors A10-A19 (excluded from the contact fits):
                                     team-A contacts and the crosslinks formed (loops / bridges)

Teammate runs (status summary_only) have no time series. They get one row in
runs_metrics.csv with inter_teamA_12A = inter-chain hits / production frames from their
crosslink_contacts_K.txt (same notebook rule, verified identical on our own runs); every
column that needs the trajectory (free-only response, SEs, QC on drift/equilibration) is
empty, so model.load() leaves them out of every fit except the team-A one.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, ROOT, SIMS, TIMESERIES, parse_contacts_summary  # noqa: E402

S = yaml.safe_load((ANALYSIS / "settings.yaml").read_text())
T0 = S["equilibration"]["discard_frames"]
NB = S["equilibration"]["n_blocks_for_sem"]
if T0 is None:
    sys.exit("settings.yaml: equilibration.discard_frames is not frozen yet")
DRIFT_Z = 3.0


def block_sem(x, nb=NB):
    blocks = np.array_split(np.asarray(x, float), nb)
    m = np.array([b.mean() for b in blocks])
    return m.std(ddof=1) / np.sqrt(nb)


def drift(x):
    a, b = np.array_split(np.asarray(x, float), 2)
    se = np.hypot(block_sem(a, NB // 2), block_sem(b, NB // 2))
    d = b.mean() - a.mean()
    return d, (d / se if se > 0 else (0.0 if d == 0 else np.inf))


def run_metrics(slug: str):
    z = np.load(TIMESERIES / f"{slug}.npz")
    defs, cuts = list(z["defs"]), z["cutoffs_A"]
    kB, kA = defs.index("B_free_lysines_only"), defs.index("A_drop_both_pinned")
    c = {v: int(np.flatnonzero(cuts == v)[0]) for v in (8.0, 12.0, 16.0)}
    n_ch = int(z["n_chains"])
    sl = slice(T0, None)
    cnt = z["counts"][sl].astype(float)

    def m_se(x):
        return x.mean(), block_sem(x)

    out = {"slug": slug, "production_frames": cnt.shape[0], "discard_frames": T0}
    prim = cnt[:, kB, c[12.0], 1]
    out["inter_free_12A"], out["inter_free_12A_se"] = m_se(prim)
    out["inter_free_12A_per_chain"] = out["inter_free_12A"] / n_ch
    out["inter_free_12A_frac_zero_frames"] = float((prim == 0).mean())
    out["inter_free_12A_drift"], out["inter_free_12A_drift_z"] = drift(prim)
    T0s = S["equilibration"]["sensitivity_discard_frames"]
    late = z["counts"][T0s:, kB, c[12.0], 1].astype(float)
    out["inter_free_12A_from200ns"], out["inter_free_12A_from200ns_se"] = late.mean(), block_sem(late)
    for name, k, cc in [("inter_free_8A", kB, 8.0), ("inter_free_16A", kB, 16.0), ("inter_teamA_12A", kA, 12.0)]:
        out[name], out[name + "_se"] = m_se(cnt[:, k, c[cc], 1])
    intra = cnt[:, kB, c[12.0], 0]
    out["intra_free_12A"], out["intra_free_12A_se"] = m_se(intra)
    out["bridge_loop_ratio_12A"] = out["inter_free_12A"] / out["intra_free_12A"] if out["intra_free_12A"] > 0 else np.nan

    lc = z["largest_cluster"][sl, kB, c[12.0]].astype(float) / n_ch
    out["largest_cluster_frac"], out["largest_cluster_frac_se"] = m_se(lc)

    # distinct free inter-chain pairs that were ever within 12 A in production (block resolution)
    blk = int(z["block"])
    first_blk = int(np.ceil(T0 / blk))
    ever = z["pair_block_hits"][first_blk:].sum(0) > 0
    sel = z["pair_keep_B"] & ~z["pair_same_chain"]
    out["inter_free_pairs_ever"] = int((ever & sel).sum())
    out["inter_free_pairs_total"] = int(sel.sum())

    h = z["z_quantiles_A"][sl] / 10.0          # nm
    out["h95_nm"], out["h95_nm_se"] = m_se(h[:, 2])
    out["h95_nm_sd_time"] = h[:, 2].std(ddof=1)
    out["h50_nm"] = h[:, 0].mean()
    out["zmean_nm"] = h[:, 4].mean()
    out["h95_drift_z"] = drift(h[:, 2])[1]

    rg = z["rg_A"][sl] / 10.0
    ree = z["ree_A"][sl] / 10.0
    for nm, a in [("rg_nm", rg), ("ree_nm", ree)]:
        out[nm], out[nm + "_se"] = m_se(a.mean(1))
        out[nm + "_sd"] = a.std(ddof=1)                       # over chains and frames
        out[nm + "_sd_between_chains"] = a.mean(0).std(ddof=1)
    out["rg_drift_z"] = drift(rg.mean(1))[1]

    deg = z["deg_inter"][sl, kB].astype(float).mean(0)       # per lysine, mean partners at 12 A
    lys = pd.DataFrame({"slug": slug, "chain": z["lys_chain"], "residue_in_chain": z["lys_within_res"],
                        "pinned": z["lys_pinned"], "mean_inter_partners_12A": deg})
    lys["k_index"] = lys.groupby("chain").cumcount() + 1     # K1..K4 from the N-terminus
    lys.loc[lys.pinned, "mean_inter_partners_12A"] = np.nan  # not in the reactive pool
    ts = pd.DataFrame({"slug": slug, "frame": np.arange(len(z["counts"])),
                       "inter_free_12A": z["counts"][:, kB, c[12.0], 1],
                       "intra_free_12A": z["counts"][:, kB, c[12.0], 0]})
    return out, lys, ts


def summary_rows(inv: pd.DataFrame) -> pd.DataFrame:
    """Team-definition-A response for runs known only from the notebook's contact summary."""
    rows = []
    for _, r in inv.iterrows():
        cs = parse_contacts_summary(SIMS / r.slug / "analysis" / "crosslink_contacts_K.txt")
        assert cs["first_frame"] == T0 and cs["cutoff_A"] == 12.0, r.slug
        rows.append({"slug": r.slug, "production_frames": cs["frames_analysed"], "discard_frames": T0,
                     "inter_teamA_12A": cs["hits_inter"] / cs["frames_analysed"],
                     "intra_teamA_12A": cs["hits_intra"] / cs["frames_analysed"],
                     "inter_teamA_pairs_ever": cs["pairs_ever_inter"], "inter_teamA_pairs_total": cs["pairs_inter"]})
    return inv.merge(pd.DataFrame(rows), on="slug")


KEEP = ["Name", "slug", "batch", "anchor_id", "replicate", "concentration", "chain_length",
        "spacing_repeats", "spacing_pct", "ionic_strength", "lead_repeats", "trail_repeats",
        "Sequence", "concentration_planned", "concentration_rel_change", "meta_n_residues",
        "meta_box_x_nm", "meta_box_z_nm", "meta_box_headroom_nm", "meta_walltime_used_fraction"]


def main():
    inv = pd.read_csv(ANALYSIS / "run_inventory.csv")
    fin = inv[inv.status == "finished"]
    res = [run_metrics(s) for s in fin.slug]
    m = pd.DataFrame([r[0] for r in res])
    df = fin[KEEP].merge(m, on="slug")
    df.insert(2, "data_source", "trajectory")
    df["design_point"] = np.where(df.anchor_id.notna(), df.anchor_id, df.slug)

    # QC flags (thresholds frozen in settings.yaml before the response was summarised)
    eq = pd.read_csv(ANALYSIS / "qc" / "equilibration_per_run.csv").pivot(index="slug", columns="observable", values="t0_frame")
    df = df.merge(eq.add_prefix("t0_detected_").reset_index(), on="slug")
    df["qc_drift_primary"] = df.inter_free_12A_drift_z.abs() > DRIFT_Z
    df["qc_drift_rg"] = df.rg_drift_z.abs() > DRIFT_Z
    df["qc_drift_h95"] = df.h95_drift_z.abs() > DRIFT_Z
    df["qc_burnin_exceeds_discard"] = df[[c for c in df if c.startswith("t0_detected_")]].max(1) > T0
    df["qc_plan_box_changed"] = df.concentration_rel_change.abs() > 0.01
    df["qc_box_headroom_lt_4nm"] = df.meta_box_headroom_nm < 4.0
    df["qc_any"] = df.filter(like="qc_").any(axis=1)
    df["excluded"] = False
    df["exclusion_reason"] = ""

    # Teammate runs: only the checks a summary allows (frames, planning, observed headroom).
    tm = summary_rows(inv[inv.status == "summary_only"])[KEEP + ["production_frames", "discard_frames",
                                                                  "inter_teamA_12A", "intra_teamA_12A",
                                                                  "inter_teamA_pairs_ever", "inter_teamA_pairs_total"]]
    tm.insert(2, "data_source", "teammate_summary")
    tm["design_point"] = tm.slug
    tm["qc_plan_box_changed"] = tm.concentration_rel_change.abs() > 0.01
    tm["qc_box_headroom_lt_4nm"] = tm.meta_box_headroom_nm < 4.0
    tm["qc_any"] = tm.qc_plan_box_changed | tm.qc_box_headroom_lt_4nm
    tm["checks_not_possible"] = "equilibration, drift, burn-in (no time series)"
    tm["excluded"] = False
    tm["exclusion_reason"] = ""
    df = pd.concat([df, tm], ignore_index=True)
    df.to_csv(ANALYSIS / "runs_metrics.csv", index=False)

    rx = inv[(inv.status == "excluded") & inv.crosslink_distance.notna()
             & inv.slug.map(lambda s: (SIMS / s / "analysis" / "crosslink_contacts_K.txt").is_file())]
    rx = summary_rows(rx)[KEEP + ["crosslink_distance", "crosslink_prob", "inter_teamA_12A", "intra_teamA_12A",
                                  "inter_teamA_pairs_ever"]]
    # Crosslinks formed during the reactive runs, from the teammate's table (his commit 0d88132;
    # the crosslink event files themselves are on /scratch/tnartey).
    xl = pd.read_csv(ROOT / "figures" / "sub4" / "crosslinks-submodule4-anchor.csv", sep=";", decimal=",",
                     encoding="utf-8-sig").rename(columns={"run": "slug", "formed": "xl_formed", "loops": "xl_loops",
                                                           "bridges": "xl_bridges"})
    rx = rx.merge(xl[["slug", "xl_formed", "xl_loops", "xl_bridges"]], on="slug", how="left", validate="1:1")
    assert rx.xl_formed.notna().all() and (rx.xl_loops + rx.xl_bridges == rx.xl_formed).all()
    free_lys = 64 - 19                                       # 19 of 64 lysines pinned (preattached 0.3)
    rx["xl_max"] = free_lys // 2                             # valence 1: at most 22 bonds
    rx["xl_conversion"] = rx.xl_formed / rx.xl_max
    rx["xl_bridge_frac"] = rx.xl_bridges / rx.xl_formed
    rx.to_csv(ANALYSIS / "reactive_anchors.csv", index=False)
    pd.concat([r[1] for r in res]).to_csv(ANALYSIS / "per_lysine.csv", index=False)
    pd.concat([r[2] for r in res]).to_csv(ANALYSIS / "timeseries_primary.csv.gz", index=False)
    print(f"{(df.data_source == 'trajectory').sum()} runs from trajectories + "
          f"{(df.data_source == 'teammate_summary').sum()} from teammate summaries; "
          f"{len(rx)} reactive anchors written separately. QC flags:")
    print(df.filter(like="qc_").fillna(False).astype(bool).sum())


if __name__ == "__main__":
    main()
