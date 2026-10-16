"""Step 2a - one pass over every finished trajectory: per-frame time series, no burn-in applied.

Writes analysis/timeseries/<slug>.npz. Everything that depends on a choice still open
(cutoff, pair definition, equilibration) is stored for every option, so freezing the
definitions afterwards never requires re-reading the trajectories, and nothing here
summarises the response per run.

Units: the DCD and top.pdb are in Angstrom (MDAnalysis), CALVADOS/OpenMM settings in nm.
z is the surface normal (one-sided wall at z_wall_nm, pins on z_wall + tether).

    python analysis/scripts/01_extract_timeseries.py [--workers 8] [--only SLUG ...]
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

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, ROOT, SIMS, TIMESERIES  # noqa: E402

CUTOFFS_A = np.array([8.0, 10.0, 12.0, 14.0, 16.0])   # Angstrom; 12 is the team convention
REF_CUTOFF_A = 12.0                                    # for the per-pair / per-lysine detail arrays
DEFS = ("A_drop_both_pinned", "B_free_lysines_only")
MIN_SEQ_SEP = 4          # residues; intra pairs closer along the chain are dropped (never binding here)
BLOCK = 100              # frames per block for the block-summed arrays (7 ns)
ZBIN_A = 2.0             # density-profile bin width (Angstrom)


def _masses():
    t = pd.read_csv(ROOT / "residues_CALVADOS2.csv")
    return dict(zip(t["three"], t["MW"]))


def _largest_cluster(adj_rows, n_chains):
    """Largest connected component size of the chain graph, for a batch of frames.

    adj_rows: (n_frames, n_pairs) bool hits over the chain pairs (i<j)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    iu, ju = np.triu_indices(n_chains, k=1)
    out = np.empty(len(adj_rows), dtype=np.int8)
    for f, hits in enumerate(adj_rows):
        if not hits.any():
            out[f] = 1
            continue
        g = csr_matrix((np.ones(hits.sum()), (iu[hits], ju[hits])), shape=(n_chains, n_chains))
        _, lab = connected_components(g, directed=False)
        out[f] = np.bincount(lab).max()
    return out


def extract(slug: str) -> str:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import self_distance_array
    warnings.filterwarnings("ignore")

    rt = SIMS / slug / "runtime"
    out = TIMESERIES / f"{slug}.npz"
    meta = pd.read_csv(rt / "metadata.csv", index_col="key")["value"]
    u = mda.Universe(str(rt / "top.pdb"), str(rt / f"{slug}.dcd"))
    n_frames = len(u.trajectory)
    chains = u.segments
    n_chains = len(chains)
    z_wall_A = float(meta["z_wall_nm"]) * 10.0

    # ---- lysines, chain membership, pins --------------------------------------------
    lys = u.select_atoms("name CA and resname LYS")
    seg = lys.segindices
    first = np.array([s.residues[0].resindex for s in chains])
    within = lys.resindices - first[seg]
    st = json.loads((rt / "surface_state.json").read_text())
    sites = [int(b) for b in st["site_beads"]]
    if sorted(sites) != sorted(lys.indices.tolist()):
        raise RuntimeError(f"{slug}: surface site beads are not the LYS beads")
    pinned_beads = {sites[int(i)] for i in st["pinned"]}
    if any(not e.get("initial") for e in st.get("events", [])):
        raise RuntimeError(f"{slug}: surface bonds formed during the run (expected none)")
    pinned = np.isin(lys.indices, list(pinned_beads))

    iu, ju = np.triu_indices(lys.n_atoms, k=1)
    same = seg[iu] == seg[ju]
    sep_ok = ~same | (np.abs(within[iu] - within[ju]) >= MIN_SEQ_SEP)
    both_p = pinned[iu] & pinned[ju]
    any_p = pinned[iu] | pinned[ju]
    keep = {DEFS[0]: sep_ok & ~both_p, DEFS[1]: sep_ok & ~any_p}

    # chain-pair index of every lysine pair, for the cluster graph
    ci, cj = np.minimum(seg[iu], seg[ju]), np.maximum(seg[iu], seg[ju])
    chain_pair_id = np.full((n_chains, n_chains), -1)
    tri = np.triu_indices(n_chains, k=1)
    chain_pair_id[tri] = np.arange(len(tri[0]))
    pair_cp = np.where(same, -1, chain_pair_id[ci, cj])

    # ---- per-chain bead masses -------------------------------------------------------
    mw = _masses()
    m_all = np.array([mw[r] for r in u.atoms.resnames], dtype=float)
    seg_atoms = [s.atoms.indices for s in chains]
    n_res = len(seg_atoms[0])

    n_c, n_d = len(CUTOFFS_A), len(DEFS)
    counts = np.zeros((n_frames, n_d, n_c, 2), dtype=np.int16)         # [..., 0]=intra, 1=inter
    largest = np.zeros((n_frames, n_d, n_c), dtype=np.int8)
    deg_inter = np.zeros((n_frames, n_d, lys.n_atoms), dtype=np.int8)   # at REF_CUTOFF_A
    n_blocks = int(np.ceil(n_frames / BLOCK))
    pair_block_hits = np.zeros((n_blocks, len(iu)), dtype=np.int16)    # at REF_CUTOFF_A, all pairs
    rg = np.zeros((n_frames, n_chains), dtype=np.float32)
    ree = np.zeros((n_frames, n_chains), dtype=np.float32)
    zq = np.zeros((n_frames, 5), dtype=np.float32)                       # q50, q90, q95, q99, mean
    zmax = np.zeros(n_frames, dtype=np.float32)
    zedges = np.arange(0.0, float(u.dimensions[2]) + ZBIN_A, ZBIN_A)
    zprof_block = np.zeros((n_blocks, len(zedges) - 1), dtype=np.float32)
    inter_hits_cp = np.zeros((n_frames, n_d, n_c, len(tri[0])), dtype=bool)
    ref_c = int(np.flatnonzero(CUTOFFS_A == REF_CUTOFF_A)[0])

    for f, ts in enumerate(u.trajectory):
        x = u.atoms.positions
        d = self_distance_array(lys.positions, box=ts.dimensions)
        hit_all = d[None, :] < CUTOFFS_A[:, None]                          # (n_c, n_pairs)
        for k, name in enumerate(DEFS):
            h = hit_all & keep[name][None, :]
            counts[f, k, :, 0] = (h & same).sum(1)
            counts[f, k, :, 1] = (h & ~same).sum(1)
            inter = h & ~same
            for c in range(n_c):
                if inter[c].any():
                    inter_hits_cp[f, k, c, pair_cp[inter[c]]] = True
            hr = inter[ref_c]
            deg_inter[f, k] = np.bincount(iu[hr], minlength=lys.n_atoms) + np.bincount(ju[hr], minlength=lys.n_atoms)
        pair_block_hits[f // BLOCK] += hit_all[ref_c]

        # height above the wall: mass-weighted quantiles of all beads
        z = x[:, 2] - z_wall_A
        o = np.argsort(z)
        cw = np.cumsum(m_all[o]) / m_all.sum()
        zq[f, :4] = z[o][np.searchsorted(cw, [0.5, 0.9, 0.95, 0.99])]
        zq[f, 4] = np.average(z, weights=m_all)
        zmax[f] = z.max()
        zprof_block[f // BLOCK] += np.histogram(x[:, 2], bins=zedges, weights=m_all)[0]

        # chain size (chains are written whole: max bond 4.6 A against a 3.8 A bond)
        for c, idx in enumerate(seg_atoms):
            p = x[idx]
            w = m_all[idx]
            com = (p * w[:, None]).sum(0) / w.sum()
            rg[f, c] = np.sqrt((w * ((p - com) ** 2).sum(1)).sum() / w.sum())
            ree[f, c] = np.linalg.norm(p[-1] - p[0])

    for k in range(n_d):
        for c in range(n_c):
            largest[:, k, c] = _largest_cluster(inter_hits_cp[:, k, c], n_chains)

    np.savez_compressed(
        out,
        cutoffs_A=CUTOFFS_A, defs=np.array(DEFS), ref_cutoff_A=REF_CUTOFF_A, block=BLOCK,
        counts=counts, largest_cluster=largest, deg_inter=deg_inter,
        pair_block_hits=pair_block_hits, pair_i=iu, pair_j=ju, pair_same_chain=same,
        pair_keep_A=keep[DEFS[0]], pair_keep_B=keep[DEFS[1]],
        lys_chain=seg, lys_within_res=within, lys_pinned=pinned,
        rg_A=rg, ree_A=ree, z_quantiles_A=zq, z_quantile_levels=np.array([0.5, 0.9, 0.95, 0.99, -1]),
        zmax_A=zmax, zprof_block=zprof_block, zprof_edges_A=zedges, z_wall_A=z_wall_A,
        n_chains=n_chains, n_res=n_res, n_frames=n_frames,
        time_per_frame_ps=float(meta["time_per_frame_ps"]), box_A=u.dimensions[:3],
    )
    return f"{slug}: {n_frames} frames, {lys.n_atoms} K ({pinned.sum()} pinned)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    inv = pd.read_csv(ANALYSIS / "run_inventory.csv")
    slugs = a.only or inv.loc[inv.status == "finished", "slug"].tolist()
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


if __name__ == "__main__":
    main()
