"""Step 1 - one read of every trajectory: heights, bond and lysine positions in z, chain size.

Writes analysis_lysplace/timeseries/<slug>.npz. Per-frame series cover all frames; production
sums (frames >= discard_frames in settings.yaml) are accumulated here so the DCDs are read once.

Units: DCD / PDB in Angstrom, heights stored in Angstrom above the wall (z - z_wall).
Atom index of a residue = chain * 357 + resid - 1 (event files use 1-based resid per chain).

    .venv/bin/python analysis_lysplace/scripts/01_extract.py [--workers 8] [--only SLUG ...]
    sbatch analysis_lysplace/scripts/extract.sbatch
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from common import MOD, ROOT, SIMS, TIMESERIES, parse_fasta  # noqa: E402

S = yaml.safe_load((MOD / "settings.yaml").read_text())
T0 = int(S["structure"]["discard_frames"])
ZBIN_A = 2.0


def bonds_table(rt: Path, n_res: int) -> pd.DataFrame:
    e = pd.read_csv(rt / "crosslink_events.csv")
    x = e[e.kind.isin(["inter", "intra"])].copy()
    x["a"] = (x.chain_i * n_res + x.resid_i - 1).astype(int)
    x["b"] = (x.chain_j * n_res + x.resid_j - 1).astype(int)
    return x[["frame", "kind", "a", "b"]].reset_index(drop=True)


def extract(slug: str) -> str:
    import MDAnalysis as mda
    warnings.filterwarnings("ignore")
    rt = SIMS / slug / "runtime"
    meta = pd.read_csv(rt / "metadata.csv", index_col="key")["value"]
    u = mda.Universe(str(rt / "top.pdb"), str(rt / f"{slug}.dcd"))
    n_frames, n_atoms = len(u.trajectory), u.atoms.n_atoms
    n_ch = len(u.segments)
    fa = parse_fasta(rt / "molecules.fasta")
    n_res = fa["length"]
    if n_atoms != n_ch * n_res:
        raise RuntimeError(f"{slug}: {n_atoms} atoms != {n_ch} x {n_res}")
    if "".join(mda.lib.util.convert_aa_code(r) for r in u.segments[0].residues.resnames) != fa["seq"]:
        raise RuntimeError(f"{slug}: topology sequence differs from molecules.fasta")
    z_wall = float(meta["z_wall_nm"]) * 10.0

    st = json.loads((rt / "surface_state.json").read_text())
    sites = np.array([int(b) for b in st["site_beads"]])
    lys = np.array(sorted(c * n_res + k for c in range(n_ch) for k in fa["k_res"]))
    if not np.array_equal(np.sort(sites), lys):
        raise RuntimeError(f"{slug}: surface sites are not the K beads")
    pinned = np.isin(lys, sites[[int(i) for i in st["pinned"]]])
    bt = bonds_table(rt, n_res)
    if not (bt.a.isin(lys).all() and bt.b.isin(lys).all()):
        raise RuntimeError(f"{slug}: a crosslink is not between two K beads")
    if np.isin(np.r_[bt.a, bt.b], lys[pinned]).any():
        raise RuntimeError(f"{slug}: a pinned lysine crosslinked")
    rgd = np.array([c * n_res + i for c in range(n_ch) for i in fa["rgd_res"]])

    # By one-letter code: the three-letter column has placeholder rows (Z = VAL, X = LYS, MW -2).
    mw = pd.read_csv(ROOT / "residues_CALVADOS2.csv").set_index("one")["MW"]
    m_ch = np.tile(mw.reindex(list(fa["seq"])).to_numpy(float), (n_ch, 1))
    if not (m_ch > 0).all():
        raise RuntimeError(f"{slug}: non-positive residue mass")
    m_all = m_ch.ravel()

    edges = np.arange(0.0, float(u.dimensions[2]) + ZBIN_A, ZBIN_A)
    nb = len(edges) - 1
    prof = {k: np.zeros(nb) for k in ("mass", "lys", "bridge", "loop", "rgd")}
    h = np.zeros((n_frames, 2), np.float32)                       # h50, h95
    rgd_t = np.zeros((n_frames, 2), np.float32)                   # mean z, fraction above h50
    rg = np.zeros((n_frames, n_ch), np.float32)
    ree = np.zeros((n_frames, n_ch), np.float32)
    nbond = len(bt)
    bz_sum = np.zeros(nbond)
    bz_rel = np.zeros(nbond)                                      # sum of z / h95
    b_above = np.zeros(nbond)
    b_n = np.zeros(nbond)
    lz_sum = np.zeros(len(lys))
    l_above = np.zeros(len(lys))
    fr0 = bt.frame.to_numpy()
    is_bridge = (bt.kind == "inter").to_numpy()
    ia, ib = bt.a.to_numpy(), bt.b.to_numpy()

    for f, ts in enumerate(u.trajectory):
        x = u.atoms.positions
        z = x[:, 2] - z_wall
        o = np.argsort(z)
        cw = np.cumsum(m_all[o]) / m_all.sum()
        h50, h95 = z[o][np.searchsorted(cw, [0.5, 0.95])]
        h[f] = h50, h95
        zr = z[rgd]
        rgd_t[f] = zr.mean(), (zr > h50).mean()
        p = x.reshape(n_ch, n_res, 3)
        com = (p * m_ch[..., None]).sum(1) / m_ch.sum(1)[:, None]
        rg[f] = np.sqrt((m_ch * ((p - com[:, None]) ** 2).sum(2)).sum(1) / m_ch.sum(1))
        ree[f] = np.linalg.norm(p[:, -1] - p[:, 0], axis=1)
        if f < T0:
            continue
        prof["mass"] += np.histogram(z, edges, weights=m_all)[0]
        prof["lys"] += np.histogram(z[lys], edges)[0]
        prof["rgd"] += np.histogram(zr, edges)[0]
        zl = z[lys]
        lz_sum += zl
        l_above += zl > h50
        on = fr0 <= f
        if on.any():
            zm = 0.5 * (z[ia[on]] + z[ib[on]])
            bz_sum[on] += zm
            bz_rel[on] += zm / h95
            b_above[on] += zm > h50
            b_n[on] += 1
            prof["bridge"] += np.histogram(zm[is_bridge[on]], edges)[0]
            prof["loop"] += np.histogram(zm[~is_bridge[on]], edges)[0]

    n_prod = n_frames - T0
    with np.errstate(invalid="ignore", divide="ignore"):
        np.savez_compressed(
            TIMESERIES / f"{slug}.npz",
            n_frames=n_frames, n_chains=n_ch, n_res=n_res, discard_frames=T0, z_wall_A=z_wall,
            time_per_frame_ps=float(meta["time_per_frame_ps"]), box_A=u.dimensions[:3],
            h_A=h, rgd_A=rgd_t, rg_A=rg, ree_A=ree,
            lys_index=lys, lys_pinned=pinned, lys_z_mean_A=lz_sum / n_prod, lys_frac_above_h50=l_above / n_prod,
            bond_frame=fr0, bond_bridge=is_bridge, bond_a=ia, bond_b=ib, bond_prod_frames=b_n,
            bond_z_mean_A=bz_sum / b_n, bond_z_rel_h95=bz_rel / b_n, bond_frac_above_h50=b_above / b_n,
            prof_edges_A=edges, **{f"prof_{k}": v for k, v in prof.items()},
        )
    return f"{slug}: {n_frames} frames, {nbond} crosslinks ({is_bridge.sum()} bridges)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    inv = pd.read_csv(MOD / "run_inventory.csv")
    slugs = a.only or inv.loc[inv.qc_ok, "slug"].tolist()
    TIMESERIES.mkdir(parents=True, exist_ok=True)
    todo = [s for s in slugs if a.overwrite or not (TIMESERIES / f"{s}.npz").is_file()]
    print(f"{len(todo)} of {len(slugs)} runs to extract with {a.workers} workers", flush=True)
    failed = []
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(extract, s): s for s in todo}
        for fut in as_completed(futs):
            try:
                print(fut.result(), flush=True)
            except Exception as e:  # noqa: BLE001 - report and carry on with the batch
                failed.append(futs[fut])
                print(f"FAILED {futs[fut]}: {e!r}", flush=True)
    print(f"done; {len(failed)} failed: {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
