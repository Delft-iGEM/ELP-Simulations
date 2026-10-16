"""Step 0 - inventory and design table: every sheet row matched to its run, design checked, QC.

Writes analysis_lysplace/run_inventory.csv (one row per run: sheet, names, design class, K
positions, design descriptors, metadata and QC). Reads the sheets, molecules.fasta, metadata.csv,
crosslink.yaml and surface.yaml; no trajectory frames and no responses.

    .venv/bin/python analysis_lysplace/scripts/00_inventory.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import MOD, SIMS, descriptors, name_positions, parse_fasta, read_sheets  # noqa: E402

EXPECTED_FRAMES = 5714
META = ["status", "n_frames_saved", "run_legs", "walltime_used_fraction", "n_residues", "box_x_nm", "box_z_nm",
        "box_headroom_nm", "box_tall_enough", "concentration_chains_per_nm2", "lattice_matches_box", "margin_nm",
        "ionic_strength_M", "temperature_K", "lysines_total", "surface_bonds", "surface_bonds_run",
        "crosslinking", "surface_seed", "z_wall_nm", "time_per_frame_ps"]


def main():
    sheets = read_sheets()
    rows = []
    for _, r in sheets.iterrows():
        rt = SIMS / r.slug / "runtime"
        rec = {"Name": r.Name, "slug": r.slug, "sheet": r.sheet, "design_class": r.design_class,
               "crosslink_distance_nm": float(r["crosslink distance"]), "crosslink_prob": float(r["crosslink prob"]),
               "mode": r["mode"], "preattached_fraction": float(r["preattached fraction"])}
        problems = []
        fa = parse_fasta(rt / "molecules.fasta")
        want = name_positions(r.Name)
        got = [b["K"] for b in fa["blocks"]]
        if want != got:
            problems.append(f"name says K at {want}, sequence has {got}")
        rec.update({"length": fa["length"], "blocks_repeats": "-".join(str(b["n_repeats"]) for b in fa["blocks"]),
                    "M_positions": "|".join("+".join(map(str, b["M"])) or "0" for b in fa["blocks"]),
                    "K_positions": "|".join("+".join(map(str, b["K"])) or "0" for b in fa["blocks"]),
                    "K_residues": "+".join(map(str, fa["k_res"]))})
        rec.update(descriptors(fa["k_res"], fa["blocks"], fa["length"]))
        meta = pd.read_csv(rt / "metadata.csv", index_col="key")["value"].reindex(META)
        rec.update({f"meta_{k}": v for k, v in meta.items()})
        xl = yaml.safe_load((rt / "crosslink.yaml").read_text())
        rec.update({f"xl_{k}": xl.get(k) for k in ("distance", "valence", "prob", "start_step", "min_span")})
        if meta["status"] != "completed":
            problems.append(f"status {meta['status']}")
        if int(float(meta["n_frames_saved"])) != EXPECTED_FRAMES:
            problems.append(f"{meta['n_frames_saved']}/{EXPECTED_FRAMES} frames")
        if float(meta["run_legs"]) > 1:
            problems.append(f"{meta['run_legs']} legs")
        if str(meta["crosslinking"]) != "True":
            problems.append("crosslinking off")
        if int(float(meta["surface_bonds"])) != 19 or int(float(meta["surface_bonds_run"] or 0)) != 0:
            problems.append(f"surface bonds {meta['surface_bonds']} (run {meta['surface_bonds_run']})")
        if str(meta["box_tall_enough"]) != "True" or str(meta["lattice_matches_box"]) != "True":
            problems.append("box / lattice check failed")
        if int(xl.get("valence", 0)) != 1:
            problems.append(f"valence {xl.get('valence')}")
        rec["qc_ok"] = not problems
        rec["qc_reason"] = "; ".join(problems) or "completed, 5714 frames, single leg, design matches sequence"
        rows.append(rec)
    inv = pd.DataFrame(rows)
    inv.to_csv(MOD / "run_inventory.csv", index=False)
    print(inv.groupby(["design_class", "structure"]).size().unstack(fill_value=0))
    print(f"\n{len(inv)} runs, {inv.qc_ok.sum()} pass QC")
    if (~inv.qc_ok).any():
        print(inv.loc[~inv.qc_ok, ["Name", "qc_reason"]].to_string(index=False))
    const = ["length", "blocks_repeats", "M_positions", "meta_concentration_chains_per_nm2", "meta_ionic_strength_M",
             "meta_temperature_K", "meta_lysines_total", "xl_distance", "xl_prob", "xl_start_step", "xl_min_span"]
    print("\nconstant across runs:", {c: inv[c].unique().tolist() for c in const})


if __name__ == "__main__":
    main()
