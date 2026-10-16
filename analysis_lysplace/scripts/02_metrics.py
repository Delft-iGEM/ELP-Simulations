"""Step 2 - per-run, per-lysine and per-crosslink tables (definitions frozen in settings.yaml).

Writes
  analysis_lysplace/runs_metrics.csv   one row per run: design, network (end of run and 100 ns),
                                       timing, heights of bridges / loops / lysines / RGD, chain size
  analysis_lysplace/per_lysine.csv     one row per lysine per run: position, state, partner, height,
                                       and the lysine-model predictors
  analysis_lysplace/per_bond.csv       one row per K-K crosslink: kind, chains, residues, formation
                                       time, height over production
  analysis_lysplace/profiles.npz       production-averaged z profiles per run (for the figures)

    analysis/scripts/run.sh analysis_lysplace/scripts/02_metrics.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import MOD, SIMS, TIMESERIES, parse_fasta  # noqa: E402

S = yaml.safe_load((MOD / "settings.yaml").read_text())
T0 = int(S["structure"]["discard_frames"])
MAXB = int(S["network"]["max_bonds"])
XL_START_NS = 10.0                       # crosslink start_step 1 000 000 x 0.01 ps (checked in the inventory)


def network(bonds: pd.DataFrame, n_ch: int) -> dict:
    """Topology of the chain graph (chains = nodes, bridges = edges, surface not a node)."""
    br = bonds[bonds.kind == "inter"]
    parent = list(range(n_ch))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    pairs = set()
    for a, b in zip(br.chain_i, br.chain_j):
        a, b = int(a), int(b)
        pairs.add((min(a, b), max(a, b)))
        parent[find(a)] = find(b)
    roots = [find(i) for i in range(n_ch)]
    sizes = np.bincount(roots, minlength=n_ch)
    n_comp = len(set(roots))
    nb, nl = len(br), int((bonds.kind == "intra").sum())
    return {"n_bridges": nb, "n_loops": nl, "n_formed": nb + nl, "bridge_frac": nb / (nb + nl) if nb + nl else np.nan,
            "n_pairs_bridged": len(pairs), "secondary_loops": nb - len(pairs),
            "largest_cluster_frac": sizes.max() / n_ch, "n_components": n_comp,
            "cycle_rank": nb - n_ch + n_comp, "conversion_available": (nb + nl) / MAXB}


def run_tables(r: pd.Series):
    rt = SIMS / r.slug / "runtime"
    z = np.load(TIMESERIES / f"{r.slug}.npz")
    n_ch, n_res = int(z["n_chains"]), int(z["n_res"])
    dt_ns = float(z["time_per_frame_ps"]) / 1000.0
    ev = pd.read_csv(rt / "crosslink_events.csv")
    bonds = ev[ev.kind.isin(["inter", "intra"])].reset_index(drop=True)
    bonds["t_ns"] = bonds.time_ps / 1000.0 - XL_START_NS
    f100 = int(np.ceil(100.0 / dt_ns)) - 1                   # last frame at or before 100 ns
    out = {"slug": r.slug, **network(bonds, n_ch),
           **{f"{k}_100ns": v for k, v in network(bonds[bonds.frame <= f100], n_ch).items()
              if k in ("n_bridges", "n_loops", "n_pairs_bridged", "largest_cluster_frac")}}
    br = bonds[bonds.kind == "inter"].sort_values("t_ns")
    if len(br):
        out["t50_bridges_ns"] = float(br.t_ns.iloc[int(np.ceil(0.5 * len(br))) - 1])
        out["t90_bridges_ns"] = float(br.t_ns.iloc[int(np.ceil(0.9 * len(br))) - 1])
        out["t_last_bond_ns"] = float(bonds.t_ns.max())

    # heights over production (bond-frame weighted means)
    h = z["h_A"][T0:] / 10.0
    out["h50_nm"], out["h95_nm"] = h[:, 0].mean(), h[:, 1].mean()
    npf, zb, zrel, ab = z["bond_prod_frames"], z["bond_z_mean_A"] / 10, z["bond_z_rel_h95"], z["bond_frac_above_h50"]
    for lab, sel in [("bridge", z["bond_bridge"]), ("loop", ~z["bond_bridge"])]:
        w = npf[sel]
        ok = w > 0
        if ok.any():
            out[f"z_{lab}_nm"] = float(np.average(zb[sel][ok], weights=w[ok]))
            out[f"z_{lab}_rel_h95"] = float(np.average(zrel[sel][ok], weights=w[ok]))
            out[f"{lab}_frac_above_h50"] = float(np.average(ab[sel][ok], weights=w[ok]))
    rgd = z["rgd_A"][T0:]
    out["z_rgd_nm"], out["rgd_frac_above_h50"] = rgd[:, 0].mean() / 10, rgd[:, 1].mean()
    out["rg_nm"], out["ree_nm"] = z["rg_A"][T0:].mean() / 10, z["ree_A"][T0:].mean() / 10

    # per-lysine table
    fa = parse_fasta(rt / "molecules.fasta")
    lys = z["lys_index"]
    chain, res = lys // n_res, lys % n_res
    pinned = z["lys_pinned"]
    state = np.where(pinned, "pinned", "unreacted").astype(object)
    partner = np.full(len(lys), -1)
    t_form = np.full(len(lys), np.nan)
    pos = {int(a): i for i, a in enumerate(lys)}
    for a, b, kind, t in zip(z["bond_a"], z["bond_b"], z["bond_bridge"], bonds.t_ns):
        for x, y in ((a, b), (b, a)):
            i = pos[int(x)]
            state[i] = "bridge" if kind else "loop"
            partner[i], t_form[i] = int(y), t
    blk = np.zeros(len(lys), int)
    pib = np.zeros(len(lys), int)
    for bi, b in enumerate(fa["blocks"]):
        for p in b["K"]:
            sel = res == b["start"] + (p - 1) * 5 + 3
            blk[sel], pib[sel] = bi + 1, p
    if (blk == 0).any():
        raise RuntimeError(f"{r.slug}: a K bead outside the blocks")
    L = pd.DataFrame({"slug": r.slug, "chain": chain, "res_index": res, "resid": res + 1, "atom_index": lys,
                      "block": blk, "pos_in_block": pib, "pinned": pinned, "state": state,
                      "partner_atom": partner, "t_form_ns": t_form,
                      "z_mean_nm": z["lys_z_mean_A"] / 10, "frac_above_h50": z["lys_frac_above_h50"]})
    L["k_order"] = L.groupby("chain").res_index.rank(method="first").astype(int)
    L["rel_pos"] = L.res_index / (n_res - 1)
    # predictors of the lysine model (sequence neighbours on the same chain)
    feats = []
    for c, g in L.groupby("chain"):
        rr, pp = g.res_index.to_numpy(), g.pinned.to_numpy()
        d = np.abs(rr[:, None] - rr[None, :]).astype(float)
        np.fill_diagonal(d, np.inf)
        near = d.min(1)
        nearest_pinned = np.array([pp[np.isclose(d[i], near[i])].any() for i in range(len(rr))])
        feats.append(pd.DataFrame({"sep_nearest_res": near, "n_k_within_30": (d < 30).sum(1),
                                   "nearest_pinned": nearest_pinned}, index=g.index))
    L = L.join(pd.concat(feats))
    out["n_unreacted"] = int((L.state == "unreacted").sum())
    for st in ("pinned", "bridge", "loop", "unreacted"):
        s = L[L.state == st]
        out[f"z_lys_{st}_nm"] = s.z_mean_nm.mean() if len(s) else np.nan

    B = pd.DataFrame({"slug": r.slug, "kind": np.where(z["bond_bridge"], "bridge", "loop"),
                      "chain_i": bonds.chain_i, "resid_i": bonds.resid_i, "chain_j": bonds.chain_j.astype(int),
                      "resid_j": bonds.resid_j.astype(int), "span_res": bonds.span, "frame": bonds.frame,
                      "t_ns": bonds.t_ns, "prod_frames": npf, "z_mean_nm": zb, "z_rel_h95": zrel,
                      "frac_above_h50": ab})
    prof = {k: z[k] for k in z.files if k.startswith("prof_")}
    return out, L, B, prof


def main():
    inv = pd.read_csv(MOD / "run_inventory.csv")
    inv = inv[inv.qc_ok]
    res = [run_tables(r) for _, r in inv.iterrows()]
    m = inv.merge(pd.DataFrame([x[0] for x in res]), on="slug")
    m["one_per_block"] = m.structure == "one_per_block"
    L = pd.concat([x[1] for x in res], ignore_index=True)
    B = pd.concat([x[2] for x in res], ignore_index=True)
    # consistency: per-lysine states reproduce the network counts
    cnt = L.groupby("slug").state.value_counts().unstack(fill_value=0)
    chk = m.set_index("slug").reindex(cnt.index)
    assert (cnt.bridge == 2 * chk.n_bridges).all() and (cnt.loop == 2 * chk.n_loops).all(), "state counts"
    assert (cnt.pinned == 19).all(), "19 pinned lysines per run"
    m.to_csv(MOD / "runs_metrics.csv", index=False)
    L.merge(m[["slug", "Name", "design_class"]], on="slug").to_csv(MOD / "per_lysine.csv", index=False)
    B.merge(m[["slug", "Name", "design_class"]], on="slug").to_csv(MOD / "per_bond.csv", index=False)
    edges = res[0][3]["prof_edges_A"]
    np.savez_compressed(MOD / "profiles.npz", slugs=m.slug.to_numpy(), edges_A=edges,
                        **{k: np.vstack([x[3][k] for x in res]) for k in res[0][3] if k != "prof_edges_A"})
    print(f"{len(m)} runs, {len(L)} lysines, {len(B)} crosslinks")
    print(m.groupby("design_class")[["n_bridges", "n_loops", "n_pairs_bridged", "largest_cluster_frac",
                                     "z_bridge_nm", "h95_nm"]].median().round(2))


if __name__ == "__main__":
    main()
