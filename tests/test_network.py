"""tools.network on hand-built graphs with known answers, plus two finished runs.

Run from the worktree root:
    /home/tnartey/ELP-Simulations/.venv/bin/python -m pytest tests/test_network.py -q

The hand-built cases are the smallest networks that separate the metrics from
each other: an inter bond (one surface-to-surface strand), a primary loop
(a cycle that is not a strand), a trifunctional star (a real junction), a
free-mode chain on one pin (everything dangling), and a chain pinned twice
(a strand that needs no junction at all). The real-run case reads finished
trajectories' event files from the main checkout, read-only, and pins the
numbers quoted in docs/NETWORK.md.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tools.network import (
    SURFACE, K_B, active_edges, analyse_runtime, build_graph, compare, junction_graph,
    kT_per_nm3_pa, lysine_resids, read_events, summary_lines, topology_metrics, write_outputs,
)

MAIN_SIMS = Path("/home/tnartey/ELP-Simulations/simulations")
SCRATCH = Path("/tmp/claude-686540/-home-tnartey-ELP-Simulations/"
               "79dd0f10-f2dd-47ea-9bd0-3dfd10be601a/scratchpad")


def bond(ci, ri, cj, rj, kind=None, **extra):
    """One K-K event row in the crosslink_events.csv schema (origin optional)."""
    row = {"frame": "1", "time_ps": "1.0", "resid_i": str(ri), "chain_i": str(ci),
           "resid_j": str(rj), "chain_j": str(cj),
           "kind": kind or ("intra" if ci == cj else "inter"),
           "span": str(abs(ri - rj)) if ci == cj else "", "x": "0", "y": "0", "z": "0",
           "pymol_i": "", "pymol_j": ""}
    row.update(extra)
    return row


def pin(ci, ri, origin="initial"):
    return {"frame": "0", "time_ps": "0.0", "resid_i": str(ri), "chain_i": str(ci),
            "resid_j": "", "chain_j": "", "kind": "surface", "span": "",
            "x": "0", "y": "0", "z": "0", "pymol_i": "", "pymol_j": "", "origin": origin}


def seq_with_k(n: int, resids: list[int], brush: bool = True) -> str:
    """Length-n sequence of G with K at the given 1-based resids (Z anchor if brush)."""
    s = ["G"] * n
    if brush:
        s[0] = "Z"
    for r in resids:
        s[r - 1] = "K"
    return "".join(s)


# ---------------------------------------------------------------------------

def test_lysine_resids_are_one_based():
    assert lysine_resids("ZPGKGVPGKG") == [4, 9]
    assert lysine_resids("KGG") == [1]


def test_kT_constant_is_computed_not_typed():
    assert math.isclose(kT_per_nm3_pa(293.15), K_B * 293.15 / 1e-27)
    assert 4.0e6 < kT_per_nm3_pa(293.15) < 4.1e6      # ≈ 4.047 MPa


def test_a_two_chains_one_inter_bond_is_one_surface_to_surface_strand():
    n, seq = 50, seq_with_k(50, [25])
    G = build_graph([bond(0, 25, 1, 25)], 2, n, "brush", seq)
    m = topology_metrics(G, box_xy_nm=(10.0, 10.0), height_nm=5.0)

    assert m["n_crosslinks"] == 1 and m["n_inter"] == 1 and m["n_intra"] == 0
    assert m["n_anchors"] == 2 and m["n_surface_bonds"] == 0
    assert m["fraction_chains_connected_to_surface"] == 1.0
    assert m["fraction_chains_surface_via_crosslinks"] == 1.0   # each reaches S via the other's anchor
    assert m["fraction_chains_two_routes_to_surface"] == 1.0    # own anchor + partner's anchor
    # S–anchor–K(0)=K(1)–anchor–S is one strand through a degree-2 point, not two.
    assert m["n_strands_active"] == 1
    assert m["n_junctions_active"] == 0
    assert m["cycle_rank"] == 1
    assert m["functionality_distribution"] == {"4": 1}
    assert m["effective_functionality_distribution"] == {"2": 1}
    # Beyond resid 25 on each chain nothing is held: 25 of 49 segments per chain.
    assert math.isclose(m["dangling_end_fraction"], 2 * 25 / (2 * 49))
    assert math.isclose(m["active_contour_fraction"], 2 * 24 / (2 * 49))
    assert m["n_strands_pruned"] == 2 and m["n_strands_rent"] == 2
    kT = kT_per_nm3_pa() / 1e3
    assert math.isclose(m["G_affine_kPa"], 1 * kT / 500.0)
    # The surface is a fixed (affine) junction, so a surface-to-surface strand
    # loses nothing to junction fluctuations: phantom == affine here.
    assert math.isclose(m["G_phantom_kPa"], m["G_affine_kPa"])
    assert math.isclose(m["strands_per_nm2"], 0.01)


def test_b_primary_loop_is_a_cycle_but_not_a_strand():
    seq = seq_with_k(50, [20, 30])
    G = build_graph([bond(0, 20, 0, 30)], 1, 50, "brush", seq)
    m = topology_metrics(G, box_xy_nm=(5.0, 5.0), height_nm=5.0)

    assert m["n_intra"] == 1 and m["n_primary_loops"] == 1 and m["n_secondary_loops"] == 0
    assert m["primary_loop_spans"] == [10] and m["primary_loop_fraction"] == 1.0
    assert m["cycle_rank"] == 1                       # the loop is a real cycle …
    assert m["n_strands_active"] == 0                 # … but carries no load
    assert m["n_junctions_active"] == 0
    assert m["effective_functionality_distribution"] == {"0": 1}
    assert m["inactive_contour_fraction"] == 1.0
    # The loop segment is not a dangling *end*, and the N-end is the anchor:
    # only the 20 segments after K30 (of 49) dangle.
    assert math.isclose(m["dangling_end_fraction"], 20 / 49)
    # The degree-≤1 pruning cannot see a pendant loop (its junction has degree 3):
    # it keeps the bridge and the loop, and RENT then subtracts the loop.
    assert m["n_strands_pruned"] == 2 and m["n_strands_rent"] == 1
    assert m["G_affine_kPa"] == 0.0 and m["G_rent_kPa"] > 0.0

    JG = junction_graph(G)
    j = JG.graph["junction_of"][(0, 20)]
    assert j == JG.graph["junction_of"][(0, 30)]
    assert JG.number_of_edges(j, j) == 1              # the primary loop is a self-loop
    assert active_edges(JG) == set()


def test_c_three_chain_star_is_one_trifunctional_junction():
    seq = seq_with_k(50, [25])
    rows = [bond(0, 25, 1, 25), bond(1, 25, 2, 25)]  # valence-2 lysine on chain 1
    G = build_graph(rows, 3, 50, "brush", seq)
    m = topology_metrics(G, box_xy_nm=(10.0, 10.0), height_nm=10.0)

    assert m["n_junctions"] == 1 and m["functionality_distribution"] == {"6": 1}
    assert m["effective_functionality_distribution"] == {"3": 1}
    assert m["n_strands_active"] == 3 and m["n_junctions_active"] == 1
    assert m["cycle_rank"] == 2
    assert m["fraction_chains_two_routes_to_surface"] == 1.0
    kT = kT_per_nm3_pa() / 1e3
    assert math.isclose(m["G_affine_kPa"], 3 * kT / 1000.0)
    assert math.isclose(m["G_phantom_kPa"], 2 * kT / 1000.0)
    assert m["largest_cluster_chains"] == 3 and m["n_chain_clusters"] == 1


def test_d_free_mode_chain_on_one_pin_is_all_dangling():
    seq = seq_with_k(50, [10, 25, 40], brush=False)
    G = build_graph([pin(0, 25)], 1, 50, "free", seq)
    m = topology_metrics(G, box_xy_nm=(5.0, 5.0), height_nm=5.0)

    assert m["n_anchors"] == 0 and m["n_surface_bonds"] == 1 and m["n_crosslinks"] == 0
    assert m["fraction_chains_connected_to_surface"] == 1.0
    assert m["fraction_chains_surface_via_crosslinks"] == 0.0
    assert m["routes_to_surface_per_chain"] == {"0": 1}
    assert m["fraction_chains_two_routes_to_surface"] == 0.0
    assert m["n_sol_chains"] == 0
    assert m["dangling_end_fraction"] == 1.0 and m["inactive_contour_fraction"] == 1.0
    assert m["n_strands_active"] == 0 and m["cycle_rank"] == 0
    assert m["G_affine_kPa"] == 0.0
    # Chain ends are degree-1 nodes in free mode (no anchor), so the graph has
    # S + 3 sites + 2 ends.
    assert G.number_of_nodes() == 6
    assert G.nodes[(0, 1)]["kind"] == "end" and G.nodes[(0, 50)]["kind"] == "end"


def test_chain_pinned_twice_is_an_active_strand_without_a_junction():
    seq = seq_with_k(50, [10, 40], brush=False)
    G = build_graph([pin(0, 10), pin(0, 40, origin="run")], 1, 50, "preattached", seq)
    m = topology_metrics(G, box_xy_nm=(5.0, 5.0), height_nm=5.0)
    assert m["n_strands_active"] == 1 and m["n_junctions_active"] == 0
    assert m["cycle_rank"] == 1
    assert math.isclose(m["active_contour_fraction"], 30 / 49)
    assert m["routes_to_surface_per_chain"] == {"0": 2}


def test_non_reactive_brush_run_has_anchors_only():
    seq = seq_with_k(40, [20])
    G = build_graph([], 4, 40, "brush", seq)
    m = topology_metrics(G)
    assert m["n_crosslinks"] == 0 and m["n_anchors"] == 4
    assert m["fraction_chains_connected_to_surface"] == 1.0
    assert m["fraction_chains_surface_via_crosslinks"] == 0.0
    assert m["n_strands_active"] == 0 and m["cycle_rank"] == 0
    assert m["G_affine_kPa"] is None and m["strands_per_nm2"] is None   # no box given
    assert math.isclose(m["dangling_end_fraction"], 1.0)


def test_sol_cluster_is_not_connected():
    # Two free-mode chains bonded to each other but pinned nowhere: a floating cluster.
    seq = seq_with_k(50, [25], brush=False)
    G = build_graph([bond(0, 25, 1, 25)], 2, 50, "free", seq)
    m = topology_metrics(G)
    assert m["n_sol_chains"] == 2
    assert m["fraction_chains_connected_to_surface"] == 0.0
    assert m["n_components"] == 2
    assert m["routes_to_surface_per_chain"] == {"0": 0, "1": 0}


def test_missing_origin_and_kind_columns_are_tolerated():
    seq = seq_with_k(50, [25])
    row = bond(0, 25, 1, 25)
    del row["kind"]
    assert "origin" not in row
    G = build_graph([row], 2, 50, "brush", seq)
    xl = [d for _, _, d in G.edges(data=True) if d["kind"] == "crosslink"]
    assert len(xl) == 1 and xl[0]["bond"] == "inter" and xl[0]["origin"] is None


def test_secondary_loop_is_not_primary():
    # K20–K40 bonded on one chain, with K30 in between bonded to another chain.
    seq = seq_with_k(60, [20, 30, 40])
    rows = [bond(0, 20, 0, 40), bond(0, 30, 1, 30)]
    m = topology_metrics(build_graph(rows, 2, 60, "brush", seq))
    assert m["n_intra"] == 1 and m["n_primary_loops"] == 0 and m["n_secondary_loops"] == 1


def test_compare_table_lists_every_run():
    seq = seq_with_k(50, [25])
    m1 = topology_metrics(build_graph([bond(0, 25, 1, 25)], 2, 50, "brush", seq), box_xy_nm=(10, 10), height_nm=5)
    m2 = topology_metrics(build_graph([], 2, 50, "brush", seq))
    table = compare([("with-bond", m1), ("bare", m2)])
    lines = table.splitlines()
    assert lines[0].startswith("run") and "G_aff" in lines[0]
    assert lines[2].startswith("with-bond") and lines[3].startswith("bare")
    assert "—" in lines[3]           # no box → no modulus, shown as a dash not a crash


# ---------------------------------------------------------------------------
# (e) finished runs from the main checkout — read-only

REAL_RUNS = {
    # name: (n_crosslinks, n_intra, n_primary_loops, nu, mu, cycle_rank)
    "xl-dense-reactive": (40, 17, 5, 70, 35, 40),
    "does-link-work2": (23, 12, 12, 17, 6, 23),
}


def _out_dir(name: str) -> Path:
    worktree_analysis = Path("simulations") / name / "analysis"
    if worktree_analysis.is_dir():
        return worktree_analysis
    out = SCRATCH / name
    out.mkdir(parents=True, exist_ok=True)
    return out


@pytest.mark.parametrize("name,expected", list(REAL_RUNS.items()))
def test_e_real_runs(name, expected):
    runtime = MAIN_SIMS / name / "runtime"
    if not (runtime / "crosslink_events.csv").is_file() or not (runtime / "metadata.csv").is_file():
        pytest.skip(f"{runtime} not available on this machine")
    rows = read_events(runtime / "crosslink_events.csv")
    assert "origin" not in rows[0]           # these files predate the column

    G, m = analyse_runtime(runtime)
    n_xl, n_intra, n_prim, nu, mu, cyc = expected
    assert m["n_crosslinks"] == n_xl and m["n_intra"] == n_intra
    assert m["n_primary_loops"] == n_prim
    assert m["n_strands_active"] == nu and m["n_junctions_active"] == mu
    assert m["cycle_rank"] == cyc
    assert m["mode"] == "brush" and m["n_sol_chains"] == 0
    assert m["fraction_chains_connected_to_surface"] == 1.0
    assert set(m["functionality_distribution"]) == {"4"}    # bifunctional linker, no terminal K
    assert m["G_affine_kPa"] > m["G_phantom_kPa"] >= 0
    assert m["height_source"].startswith("ideal coil")     # old metadata has no z-extent

    out = _out_dir(name)
    write_outputs(m, out / "network_summary.txt", out / "network.json")
    back = json.loads((out / "network.json").read_text())
    assert back["n_strands_active"] == nu
    assert any("active strands nu" in line for line in summary_lines(m))
