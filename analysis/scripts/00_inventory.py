"""Step 1 - inventory: every run-sheet row matched to its output, with a status and reason.

Writes analysis/run_inventory.csv. Reads only run sheets, metadata.csv, timing.yaml and
the DCD header; no trajectory frames and no responses. For teammate runs without a
runtime/ folder it reads load_summary.txt and the frame window of crosslink_contacts_K.txt
(not the counts) to decide whether the summary is usable.
"""
from __future__ import annotations

import pandas as pd

from common import (ANALYSIS, LOCAL_SHEETS, SIMS, TEAMMATE_SHEETS, anchor_fields, batch_of, parse_contacts_summary,
                    parse_load_summary, parse_sequence, read_sheet, slugify)

EXPECTED_FRAMES = 5714
META_KEYS = ["status", "n_residues", "n_frames_saved", "run_legs", "walltime_used_fraction",
             "box_x_nm", "box_z_nm", "z_max_nm", "box_headroom_nm", "box_tall_enough",
             "spacing_nm", "concentration_chains_per_nm2", "lattice_matches_box", "margin_nm",
             "ionic_strength_M", "lysines_total", "surface_bonds", "surface_bonds_run",
             "surface_chains_with_1", "surface_chains_with_2", "surface_chains_with_3plus",
             "crosslinking", "surface_seed", "total_time_ns"]

sheets = [read_sheet(p) for p in LOCAL_SHEETS + TEAMMATE_SHEETS]
rows = pd.concat(sheets, ignore_index=True)
rows["slug"] = rows["Name"].map(slugify)
assert not rows["slug"].duplicated().any(), "a run name appears in two sheets"

records = []
for _, r in rows.iterrows():
    rec = {"Name": r["Name"], "slug": r["slug"], "sheet": r["sheet"], "batch": batch_of(r["Name"]),
           **anchor_fields(r["Name"]), **parse_sequence(r["Sequence"]),
           "Sequence": r["Sequence"], "concentration_planned": r["concentration"],
           "ionic_strength": r["ionic strength (M)"], "steps": r["steps"], "mode": r["mode"],
           "preattached_fraction": r["preattached fraction"], "walltime": r["walltime"],
           "crosslink_distance": r["crosslink distance"], "crosslink_prob": r["crosslink prob"]}
    meta_path = SIMS / r["slug"] / "runtime" / "metadata.csv"
    reactive = isinstance(r["crosslink distance"], str) and r["crosslink distance"].strip() != ""
    if meta_path.is_file():
        meta = pd.read_csv(meta_path, index_col="key")["value"].reindex(META_KEYS)
        rec.update({f"meta_{k}": v for k, v in meta.items()})
        dcd = SIMS / r["slug"] / "runtime" / f"{r['slug']}.dcd"
        frames = 0 if pd.isna(meta["n_frames_saved"]) else int(float(meta["n_frames_saved"]))   # empty while a run is unfinished
        conc = float(meta["concentration_chains_per_nm2"])
        rec["concentration"] = conc
        rec["concentration_rel_change"] = conc / r["concentration"] - 1
        problems = []
        if not dcd.is_file():
            problems.append("no DCD")
        if meta["status"] != "completed":
            problems.append(f"status {meta['status']}")
        if frames < EXPECTED_FRAMES:
            problems.append(f"{frames}/{EXPECTED_FRAMES} frames")
        if pd.notna(meta["run_legs"]) and float(meta["run_legs"]) > 1:
            problems.append(f"{meta['run_legs']} legs (resumed after walltime)")
        if meta["status"] == "completed" and str(meta["box_tall_enough"]) != "True":
            problems.append("box not tall enough")
        if meta["status"] == "completed" and (str(meta["lattice_matches_box"]) != "True" or float(meta["margin_nm"]) != 0):
            problems.append("plan_box changed the lattice/box")
        if abs(rec["concentration_rel_change"]) > 0.01:
            problems.append("concentration changed by >1% at planning")
        if reactive:
            problems.append("reactive K-K crosslinking was on")
        rec["status"] = "finished" if not problems else ("excluded" if reactive else "incomplete")
        rec["reason"] = "; ".join(problems) or "completed, 5714 frames, single leg, box checks pass"
    elif (SIMS / r["slug"] / "analysis" / "crosslink_contacts_K.txt").is_file():
        # Teammate run: trajectory not readable, the notebook's per-run summaries are.
        an = SIMS / r["slug"] / "analysis"
        ld = parse_load_summary(an / "load_summary.txt")
        cs = parse_contacts_summary(an / "crosslink_contacts_K.txt")
        rec.update({"meta_n_residues": ld["n_residues"], "meta_box_x_nm": ld["box_x_nm"],
                    "meta_box_z_nm": ld["box_z_nm"], "meta_z_max_nm": ld["z_max_nm"],
                    "meta_box_headroom_nm": ld["box_headroom_nm"], "meta_n_frames_saved": cs["frames_total"]})
        conc = ld["n_chains"] / (ld["box_x_nm"] * ld["box_y_nm"])   # lattice tiles the box, no margin
        rec["concentration"] = conc
        rec["concentration_rel_change"] = conc / r["concentration"] - 1
        problems = []
        if ld["n_residues"] != 5 * rec["chain_length"]:
            problems.append(f"{ld['n_residues']} residues/chain, sheet says {5 * rec['chain_length']}")
        if cs["frames_total"] < EXPECTED_FRAMES:
            problems.append(f"{cs['frames_total']}/{EXPECTED_FRAMES} frames")
        if cs["cutoff_A"] != 12.0 or cs["first_frame"] != 1428:
            problems.append(f"summary at {cs['cutoff_A']} A from frame {cs['first_frame']}, not 12 A from 1428")
        if abs(rec["concentration_rel_change"]) > 0.01:
            problems.append("concentration changed by >1% at planning")
        if reactive:
            problems.append(f"reactive K-K crosslinking on (distance {r['crosslink distance']} nm, "
                            f"prob {r['crosslink prob']}) - protocol requires it off")
        rec["status"] = ("excluded" if reactive else "incomplete") if problems else "summary_only"
        rec["reason"] = "; ".join(problems + ["trajectory on /scratch/tnartey (permission denied); "
                                              "only the notebook summaries (team definition A, 12 A) are available"])
    else:
        rec["concentration"] = r["concentration"]
        rec["status"] = "missing"
        why = ["no simulations/<name>/ folder"]
        if reactive:
            why.append(f"reactive K-K crosslinking on (distance {r['crosslink distance']} nm, "
                       f"prob {r['crosslink prob']}) - protocol requires it off; not comparable "
                       f"even if the data were available")
            rec["status"] = "excluded"
        rec["reason"] = "; ".join(why)
    records.append(rec)

inv = pd.DataFrame(records)
inv.to_csv(ANALYSIS / "run_inventory.csv", index=False)
print(inv.groupby(["batch", "status"]).size().unstack(fill_value=0))
print(f"\n{len(inv)} rows; analysable: {(inv.status == 'finished').sum()} from trajectories, "
      f"{(inv.status == 'summary_only').sum()} from teammate summaries (team definition A only)")
print("\nPlanned in the protocol but in no sheet available here: sub4-main-* (100 rows), "
      "sub4-extra5 rows 15-20 (6 rows).")
