"""Network topology of a crosslinked ELP coating — what the crosslink events add up to.

``runtime/crosslink_events.csv`` is a list of bonds: lysine i on chain a met
lysine j on chain b (or the surface) and reacted. On its own that list answers
"how many" and "intra or inter", which is what ``crosslink_summary.txt``
reports. It does not answer the question the project is asking — *what kind of
network did these bonds make, and would it be stiff?* — because stiffness is a
property of the connected structure, not of the bond count: a chain bonded
back onto itself contributes nothing, a cluster of chains that never reaches the
surface is sol, and a strand only carries load if it sits between two points
that are themselves held.

This module turns the event list into a graph and asks that structural
question with the classical rubber-elasticity bookkeeping:

* Miller & Macosko (Macromolecules 1976, 9, 206): an *elastically effective*
  junction has at least three arms that lead independently to the infinite
  network; an *elastically effective strand* connects two such junctions.
* Flory's phantom network (Proc. R. Soc. A 1976, 351, 351) and the affine
  model: ``G_aff = ν kT / V``, ``G_ph = (ν − μ) kT / V`` with ν effective strands
  and μ effective junctions per volume.
* Zhong, Wang, Johnson et al. (Science 2016, 353, 1264) — real elastic network
  theory (RENT): primary loops (a strand bonded back to its own junction) are
  the dominant defect and each one removes elastically effective material.
* Lange, Schwenke, Kurakazu et al. (Macromolecules 2011, 44, 9666) on the
  defect census — dangling ends, loops, sol fraction — of model networks.

The surface plays the role of the "infinite network": in a grafted coating
every chain is held by the substrate, so a strand is under load when the
coating is sheared iff it lies on a path that leaves the surface and comes back
to it. ``SURFACE`` is one node of infinite functionality that is never pruned.

Two graphs
----------
``build_graph`` gives the *site graph*: every lysine (and, in brush mode, the
residue-0 anchor) is a node carrying its chain and residue, the backbone
between consecutive sites is a ``strand`` edge with its residue count, each
K–K bond is a ``crosslink`` edge and each K–surface bond (or brush anchor) a
``surface`` edge to ``SURFACE``. That is the faithful picture and what the
tests build by hand.

The elastic bookkeeping is done on the *junction graph* (``junction_graph``):
every crosslink and surface edge is contracted — the two lysines of a bond are
one junction, a pinned lysine *is* the surface — so the only edges left are
backbone strands, exactly the objects that carry load. Contraction preserves
the cycle rank, and it makes the defects fall out geometrically: a primary loop
becomes a self-loop on its junction, and a junction whose only remaining arms
are two backbone strands has degree 2 and is not a junction at all — precisely
the RENT statement that a primary loop demotes an f=4 junction to f=2.

A strand is *active* iff it lies on a simple cycle through ``SURFACE``, i.e. it
belongs to a biconnected block that contains ``SURFACE``. That single criterion
is Miller–Macosko's "leads to the infinite network in at least two independent
ways" made exact for a finite graph: a pendant loop, a dangling end, a bridge
from the surface to a structure that is held nowhere else, and a cluster that
floats free are all excluded, without needing the degree-≤1 pruning heuristic
(which cannot see a pendant loop, because its junction has degree 3). The
pruning count is still reported, together with the literature-style
"subtract the primary loops" correction, so the two can be compared.

Read before believing the moduli
--------------------------------
* The bonds come from a distance criterion on coarse-grained beads, at one
  cutoff and one solvent quality (``tools.crosslink`` says why). The topology
  is meaningful; the kinetics are not; the conversion is set by the cutoff.
* ``G`` needs a volume, and a coating's height is ill-defined: the default
  takes the trajectory's bead z-extent from metadata when it is there and
  falls back to the ideal-coil size. Prefer the per-area numbers
  (strands per nm²) when comparing coatings of the same footprint, and treat
  every G as an order of magnitude for *ranking sequences*, not a prediction.
* Entanglements, finite extensibility, the excluded-volume contribution and
  the (weak) elasticity of pendant material are all ignored, as they are in
  the phantom model.
* The runs are small (9–16 chains). Percolation statistics from a 3x3 grafting
  lattice are what they are; use the metrics to compare sequences run under
  identical conditions.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import networkx as nx
from networkx.algorithms.connectivity import local_node_connectivity

SURFACE = "SURFACE"

# Boltzmann constant, J/K (CODATA 2018, exact). kT/nm^3 in Pa is
# k_B * T / 1e-27; at 293.15 K that is 4.047e6 Pa — not 4.07e6 (that would be
# T ≈ 295 K) — so the constant is computed, not typed.
K_B = 1.380649e-23
DEFAULT_T_K = 293.15
NM3_TO_M3 = 1e-27

SUMMARY_FILENAME = "network_summary.txt"
JSON_FILENAME = "network.json"

# Edge kinds in the site graph.
STRAND, CROSSLINK, SURFACE_BOND = "strand", "crosslink", "surface"


def kT_per_nm3_pa(temperature_K: float = DEFAULT_T_K) -> float:
    """kT / nm^3 in pascal — the unit every modulus estimate is a multiple of."""
    return K_B * temperature_K / NM3_TO_M3


# ---------------------------------------------------------------------------
# Reading inputs
# ---------------------------------------------------------------------------

def lysine_resids(sequence: str) -> list[int]:
    """1-based residue numbers (PDB resSeq, as in the events file) of every K.

    ``tools.crosslink.SiteInfo.resid`` is ``residue.resSeq``, which CALVADOS
    numbers from 1, so sequence index ``i`` is resid ``i + 1``. The brush anchor
    'Z' at index 0 is resid 1.
    """
    return [i + 1 for i, aa in enumerate(sequence) if aa.upper() == "K"]


def read_events(path: Path | str) -> list[dict[str, str]]:
    """Rows of a crosslink_events.csv, or [] if the file is absent (non-reactive run).

    Older files lack the ``origin`` column (and any file may gain columns);
    everything downstream reads by name with ``.get`` so that is harmless.
    """
    path = Path(path)
    if not path.is_file():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# The site graph
# ---------------------------------------------------------------------------

def build_graph(rows: Iterable[Mapping[str, Any]], n_chains: int, n_residues: int,
                mode: str = "brush", sequence: str | None = None,
                sites: Sequence[int] | None = None) -> nx.MultiGraph:
    """Site graph of a run: lysines, anchors, chain ends and SURFACE.

    Nodes are ``(chain, resid)`` tuples (0-based chain, 1-based resid) with
    attributes ``chain``, ``resid`` and ``kind`` in {'anchor', 'site', 'end'},
    plus the single node ``SURFACE``. Edges carry ``kind``:

    * ``strand`` — backbone between consecutive nodes of a chain, with
      ``chain`` and ``n_res`` (residues along the contour, ``resid_b − resid_a``);
    * ``crosslink`` — a K–K bond from the events (``intra``/``inter`` kept as
      ``bond``), with ``frame``, ``time_ps`` and ``origin`` when present;
    * ``surface`` — a K–surface bond, or the brush anchor at residue 0.

    ``sites`` overrides the reactive residues (default: the K's of ``sequence``);
    any residue that appears in the events is added as a site regardless, so a
    run with an explicit ``crosslink_selection`` still builds. Chain ends that
    are not sites become ``end`` nodes of degree 1, so a dangling terminal
    segment is an edge with a residue count like any other strand.

    A MultiGraph, because a primary loop (two adjacent sites of one chain bonded
    together) is a strand and a crosslink between the same two nodes, and that
    parallel pair is exactly the cycle the loop metrics need to see.
    """
    rows = list(rows)
    if sites is None:
        if sequence is None:
            raise ValueError("build_graph needs either `sequence` or `sites`.")
        sites = lysine_resids(sequence)
    brush = mode == "brush"

    # Which residues are nodes on each chain: the sites, any resid the events
    # mention, the anchor (brush) and the two termini.
    per_chain: dict[int, set[int]] = {c: set(sites) for c in range(n_chains)}
    site_kind: dict[tuple[int, int], str] = {}
    for c in range(n_chains):
        for r in sites:
            site_kind[(c, r)] = "site"
    for row in rows:
        for suffix in ("i", "j"):
            c, r = _int(row.get(f"chain_{suffix}")), _int(row.get(f"resid_{suffix}"))
            if c is None or r is None:
                continue
            per_chain.setdefault(c, set()).add(r)
            site_kind.setdefault((c, r), "site")
    for c in per_chain:
        if brush:
            per_chain[c].add(1)
            site_kind[(c, 1)] = "anchor"
        else:
            per_chain[c].add(1)
            site_kind.setdefault((c, 1), "end")
        per_chain[c].add(n_residues)
        site_kind.setdefault((c, n_residues), "end")

    G = nx.MultiGraph(n_chains=n_chains, n_residues=n_residues, mode=mode)
    G.add_node(SURFACE, kind="surface")
    for c, resids in per_chain.items():
        ordered = sorted(resids)
        for r in ordered:
            G.add_node((c, r), chain=c, resid=r, kind=site_kind[(c, r)])
        for a, b in zip(ordered, ordered[1:]):
            G.add_edge((c, a), (c, b), kind=STRAND, chain=c, n_res=b - a)
        if brush:
            G.add_edge((c, 1), SURFACE, kind=SURFACE_BOND, chain=c, origin="anchor")

    for row in rows:
        kind = (row.get("kind") or "").strip()
        ci, ri = _int(row.get("chain_i")), _int(row.get("resid_i"))
        if ci is None or ri is None:
            continue
        attrs = {"frame": _int(row.get("frame")), "time_ps": row.get("time_ps"),
                 "origin": row.get("origin") or None}
        if kind == SURFACE_BOND:
            G.add_edge((ci, ri), SURFACE, kind=SURFACE_BOND, chain=ci, **attrs)
            continue
        cj, rj = _int(row.get("chain_j")), _int(row.get("resid_j"))
        if cj is None or rj is None:
            continue
        if not kind:
            kind = "intra" if ci == cj else "inter"
        G.add_edge((ci, ri), (cj, rj), kind=CROSSLINK, bond=kind, **attrs)
    return G


# ---------------------------------------------------------------------------
# The junction graph — crosslinks and pins contracted, only strands remain
# ---------------------------------------------------------------------------

def junction_graph(G: nx.MultiGraph) -> nx.MultiGraph:
    """Contract every crosslink and surface edge of a site graph.

    Nodes are ``SURFACE`` (everything pinned to it, anchors included) and
    ``('J', k)`` junctions with a ``members`` list of site-graph nodes; an
    unbonded lysine is a one-member junction of degree 2 and a chain end a
    one-member node of degree 1. Edges are the backbone strands with their
    ``chain``, ``resid_a``, ``resid_b`` and ``n_res``. ``JG.graph['junction_of']``
    maps every site-graph node to its junction.
    """
    parent: dict[Any, Any] = {n: n for n in G.nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        # Keep SURFACE as its own root so the merged group is recognisable.
        if rb == SURFACE:
            ra, rb = rb, ra
        parent[rb] = ra

    for u, v, data in G.edges(data=True):
        if data.get("kind") in (CROSSLINK, SURFACE_BOND):
            union(u, v)

    groups: dict[Any, list] = defaultdict(list)
    for n in G.nodes:
        groups[find(n)].append(n)

    JG = nx.MultiGraph(**G.graph)
    junction_of: dict[Any, Any] = {}
    k = 0
    for root, members in groups.items():
        if root == SURFACE:
            node = SURFACE
        else:
            node = ("J", k)
            k += 1
        chains = sorted({G.nodes[m]["chain"] for m in members if m != SURFACE})
        JG.add_node(node, members=sorted(m for m in members if m != SURFACE),
                    chains=chains, kind="surface" if node == SURFACE else "junction")
        for m in members:
            junction_of[m] = node

    for u, v, data in G.edges(data=True):
        if data.get("kind") != STRAND:
            continue
        ra, rb = sorted((G.nodes[u]["resid"], G.nodes[v]["resid"]))
        JG.add_edge(junction_of[u], junction_of[v], chain=data["chain"],
                    resid_a=ra, resid_b=rb, n_res=data["n_res"])
    JG.graph["junction_of"] = junction_of
    return JG


def active_edges(JG: nx.MultiGraph) -> set[tuple[Any, Any, Any]]:
    """Edge keys ``(u, v, key)`` of the strands that carry load.

    A strand is active iff it lies on a simple cycle through SURFACE, which is
    the same as lying in a biconnected block that contains SURFACE. Blocks are
    defined for simple graphs, so every edge is first subdivided into
    ``u – m1 – m2 – v``; that turns a self-loop into a triangle and parallel
    edges into a hexagon without changing which edges share a cycle.
    """
    if SURFACE not in JG:
        return set()
    H = nx.Graph()
    H.add_nodes_from(JG.nodes)
    marker: dict[Any, tuple] = {}
    for i, (u, v, key) in enumerate(JG.edges(keys=True)):
        m1, m2 = ("m", i, 1), ("m", i, 2)
        H.add_edge(u, m1)
        H.add_edge(m1, m2)
        H.add_edge(m2, v)
        marker[m1] = (u, v, key)
    active = set()
    for block in nx.biconnected_components(H):
        # A block with ≤ 2 nodes is a bridge; bigger blocks contain a cycle.
        if SURFACE in block and len(block) > 2:
            active.update(marker[m] for m in block if m in marker)
    return active


def _prune(JG: nx.MultiGraph) -> nx.MultiGraph:
    """Miller–Macosko style pruning: drop degree ≤ 1 nodes (never SURFACE) until stable."""
    P = JG.copy()
    while True:
        drop = [n for n, d in P.degree() if d <= 1 and n != SURFACE]
        if not drop:
            return P
        P.remove_nodes_from(drop)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _chain_nodes(G: nx.MultiGraph) -> dict[int, list]:
    out: dict[int, list] = defaultdict(list)
    for n, data in G.nodes(data=True):
        if n != SURFACE:
            out[data["chain"]].append(n)
    return out


def _percolation(G: nx.MultiGraph, chains: Sequence[int]) -> tuple[dict[int, bool], dict[int, int]]:
    """Per chain: reaches SURFACE without its own pins?  Number of node-disjoint routes?

    The first question removes the chain's own ``surface`` edges (anchor and
    pinned lysines) and asks whether any of its nodes still reaches SURFACE —
    i.e. whether the *crosslinks* tie it into something held. The second
    contracts each chain to one node, puts a dummy 'pin' node on every surface
    edge (so no chain is adjacent to SURFACE and every attachment is its own
    node), and takes the node connectivity between chain and SURFACE: the
    number of node-disjoint paths, which is the number of independent ways the
    chain is held (own anchor counts as one). Node connectivity is a max-flow
    problem — O(V·E) per chain with networkx's default algorithm, negligible at
    the sizes here but the one thing in this module that is not linear.
    """
    by_chain = _chain_nodes(G)
    via_xl: dict[int, bool] = {}
    for c in chains:
        H = nx.Graph()
        H.add_nodes_from(G.nodes)
        for u, v, data in G.edges(data=True):
            if data.get("kind") == SURFACE_BOND and data.get("chain") == c:
                continue
            H.add_edge(u, v)
        via_xl[c] = any(nx.has_path(H, n, SURFACE) for n in by_chain.get(c, []))

    P = nx.Graph()
    P.add_node(SURFACE)
    for c in chains:
        P.add_node(("C", c))
    for i, (u, v, data) in enumerate(G.edges(data=True)):
        kind = data.get("kind")
        if kind == STRAND:
            continue
        if kind == SURFACE_BOND:
            site = u if v == SURFACE else v
            pin = ("pin", i)
            P.add_edge(("C", G.nodes[site]["chain"]), pin)
            P.add_edge(pin, SURFACE)
        elif kind == CROSSLINK:
            cu, cv = G.nodes[u]["chain"], G.nodes[v]["chain"]
            if cu != cv:
                P.add_edge(("C", cu), ("C", cv))
    routes: dict[int, int] = {}
    for c in chains:
        node = ("C", c)
        if not nx.has_path(P, node, SURFACE):
            routes[c] = 0
        else:
            routes[c] = int(local_node_connectivity(P, node, SURFACE))
    return via_xl, routes


def topology_metrics(G: nx.MultiGraph, *, box_xy_nm: tuple[float, float] | None = None,
                     height_nm: float | None = None, height_source: str = "given",
                     temperature_K: float = DEFAULT_T_K) -> dict[str, Any]:
    """Every metric documented in docs/NETWORK.md, as one flat JSON-friendly dict.

    ``box_xy_nm`` and ``height_nm`` are only needed for the moduli; without
    them the counts and fractions are still returned and the G's are None.
    """
    n_chains = G.graph["n_chains"]
    n_res = G.graph["n_residues"]
    mode = G.graph.get("mode", "brush")
    chains = list(range(n_chains))
    by_chain = _chain_nodes(G)

    # ---- bond census ------------------------------------------------------
    xl_edges = [(u, v, d) for u, v, d in G.edges(data=True) if d.get("kind") == CROSSLINK]
    surf_edges = [(u, v, d) for u, v, d in G.edges(data=True)
                  if d.get("kind") == SURFACE_BOND and d.get("origin") != "anchor"]
    n_anchor = sum(1 for _, _, d in G.edges(data=True)
                   if d.get("kind") == SURFACE_BOND and d.get("origin") == "anchor")
    intra = [(u, v, d) for u, v, d in xl_edges if G.nodes[u]["chain"] == G.nodes[v]["chain"]]
    inter = [(u, v, d) for u, v, d in xl_edges if G.nodes[u]["chain"] != G.nodes[v]["chain"]]

    bonds_per_chain = Counter({c: 0 for c in chains})
    for u, v, _ in xl_edges:
        bonds_per_chain[G.nodes[u]["chain"]] += 1
        if G.nodes[v]["chain"] != G.nodes[u]["chain"]:
            bonds_per_chain[G.nodes[v]["chain"]] += 1
    surface_per_chain = Counter({c: 0 for c in chains})
    for u, v, _ in surf_edges:
        site = u if v == SURFACE else v
        surface_per_chain[G.nodes[site]["chain"]] += 1

    # Attached residues per chain: anything that is bonded to something
    # (crosslink, pin or anchor). A site with only backbone edges is free.
    attached: dict[int, set[int]] = {c: set() for c in chains}
    for u, v, d in G.edges(data=True):
        if d.get("kind") == STRAND:
            continue
        for n in (u, v):
            if n != SURFACE:
                attached[G.nodes[n]["chain"]].add(G.nodes[n]["resid"])

    # ---- loops ------------------------------------------------------------
    primary_spans, secondary_spans = [], []
    for u, v, _ in intra:
        c = G.nodes[u]["chain"]
        a, b = sorted((G.nodes[u]["resid"], G.nodes[v]["resid"]))
        inside = [r for r in attached[c] if a < r < b]
        (primary_spans if not inside else secondary_spans).append(b - a)
    primary_spans.sort()
    n_xl = len(xl_edges)

    # ---- dangling ends (chain segments beyond the last attachment) ---------
    # In brush mode the anchor at resid 1 is an attachment, so a brush chain never
    # has a dangling N-end; a chain with no attachment at all is dangling entirely.
    contour_per_chain = n_res - 1          # backbone segments in one chain
    dangling_segments = 0
    for c in chains:
        if not attached[c]:
            dangling_segments += contour_per_chain
            continue
        dangling_segments += (min(attached[c]) - 1) + (n_res - max(attached[c]))
    total_contour = n_chains * contour_per_chain

    # ---- junction graph and the elastic census ----------------------------
    JG = junction_graph(G)
    active = active_edges(JG)
    active_deg: Counter = Counter()
    active_contour = 0
    for u, v, key in active:
        active_deg[u] += 1
        active_deg[v] += 1            # a self-loop on SURFACE counts twice, as it should
        active_contour += JG.edges[u, v, key]["n_res"]
    n_pass_through = sum(1 for n, d in active_deg.items() if n != SURFACE and d == 2)
    nu = len(active) - n_pass_through
    mu = sum(1 for n, d in active_deg.items() if n != SURFACE and d >= 3)

    junctions = [n for n in JG.nodes if n != SURFACE and JG.degree(n) >= 3]
    functionality = Counter(JG.degree(n) for n in junctions)
    eff_functionality = Counter(active_deg.get(n, 0) for n in junctions)

    n_components = nx.number_connected_components(G)
    # Chain clusters through crosslinks alone (SURFACE would join everything in brush mode).
    CG = nx.Graph()
    CG.add_nodes_from(chains)
    for u, v, _ in inter:
        CG.add_edge(G.nodes[u]["chain"], G.nodes[v]["chain"])
    clusters = sorted((len(cc) for cc in nx.connected_components(CG)), reverse=True)

    surf_component = nx.node_connected_component(G, SURFACE)
    connected = {c: any(n in surf_component for n in by_chain.get(c, [])) for c in chains}
    # Sol = not held by the substrate at all, whether alone or in a floating cluster.
    n_sol = sum(1 for c in chains if not connected[c])
    via_xl, routes = _percolation(G, chains)

    # Pruning count (the textbook heuristic) and its RENT correction, for comparison.
    pruned = _prune(JG)
    n_strands_pruned = pruned.number_of_edges()
    n_junctions_pruned = sum(1 for n, d in pruned.degree() if n != SURFACE and d >= 3)
    nu_rent = n_strands_pruned - len(primary_spans)

    cycle_rank = G.number_of_edges() - G.number_of_nodes() + n_components

    # ---- moduli -----------------------------------------------------------
    kT_nm3 = kT_per_nm3_pa(temperature_K)
    area = box_xy_nm[0] * box_xy_nm[1] if box_xy_nm else None
    volume = area * height_nm if (area and height_nm) else None

    def modulus_kpa(count: int) -> float | None:
        if not volume:
            return None
        return max(count, 0) * kT_nm3 / volume / 1e3

    def frac(num: float, den: float) -> float:
        return num / den if den else 0.0

    return {
        "mode": mode,
        "n_chains": n_chains,
        "n_residues": n_res,
        "n_sites_per_chain": sum(1 for n, d in G.nodes(data=True)
                                 if n != SURFACE and d["kind"] == "site" and d["chain"] == 0),
        "n_crosslinks": n_xl,
        "n_intra": len(intra),
        "n_inter": len(inter),
        "n_surface_bonds": len(surf_edges),
        "n_anchors": n_anchor,
        "bonds_per_chain": {str(c): bonds_per_chain[c] for c in chains},
        "bonds_per_chain_distribution": {str(k): v for k, v in sorted(Counter(bonds_per_chain.values()).items())},
        "surface_bonds_per_chain_distribution": {str(k): v for k, v in sorted(Counter(surface_per_chain.values()).items())},
        "n_primary_loops": len(primary_spans),
        "n_secondary_loops": len(secondary_spans),
        "primary_loop_fraction": frac(len(primary_spans), n_xl),
        "primary_loop_spans": primary_spans,
        "small_loop_span_min": primary_spans[0] if primary_spans else None,
        "small_loop_span_median": (primary_spans[len(primary_spans) // 2] if primary_spans else None),
        "n_components": n_components,
        "n_chain_clusters": len(clusters),
        "largest_cluster_chains": clusters[0] if clusters else 0,
        "cluster_sizes": clusters,
        "fraction_chains_connected_to_surface": frac(sum(connected.values()), n_chains),
        "fraction_chains_surface_via_crosslinks": frac(sum(via_xl.values()), n_chains),
        "fraction_chains_two_routes_to_surface": frac(sum(1 for r in routes.values() if r >= 2), n_chains),
        "routes_to_surface_per_chain": {str(c): routes[c] for c in chains},
        "n_sol_chains": n_sol,
        "dangling_end_segments": dangling_segments,
        "dangling_end_fraction": frac(dangling_segments, total_contour),
        "active_contour_fraction": frac(active_contour, total_contour),
        "inactive_contour_fraction": 1.0 - frac(active_contour, total_contour),
        "cycle_rank": cycle_rank,
        "n_strands_active": nu,
        "n_junctions_active": mu,
        "n_junctions": len(junctions),
        "functionality_distribution": {str(k): v for k, v in sorted(functionality.items())},
        "effective_functionality_distribution": {str(k): v for k, v in sorted(eff_functionality.items())},
        "n_strands_pruned": n_strands_pruned,
        "n_junctions_pruned": n_junctions_pruned,
        "n_strands_rent": nu_rent,
        "strands_per_chain": frac(nu, n_chains),
        "strands_per_nm2": (nu / area) if area else None,
        "junctions_per_nm2": (mu / area) if area else None,
        "box_x_nm": box_xy_nm[0] if box_xy_nm else None,
        "box_y_nm": box_xy_nm[1] if box_xy_nm else None,
        "height_nm": height_nm,
        "height_source": height_source if height_nm else None,
        "volume_nm3": volume,
        "temperature_K": temperature_K,
        "kT_per_nm3_MPa": kT_nm3 / 1e6,
        "G_affine_kPa": modulus_kpa(nu),
        "G_phantom_kPa": modulus_kpa(nu - mu),
        "G_rent_kPa": modulus_kpa(nu_rent),
    }


# ---------------------------------------------------------------------------
# Running on a simulation folder
# ---------------------------------------------------------------------------

def coating_height(meta: Mapping[str, str], height_nm: float | None = None) -> tuple[float | None, str]:
    """Height to turn the footprint into a volume, and where it came from.

    Explicit argument > trajectory bead extent (``z_max_nm − z_min_nm`` from
    ``sim metadata``) > ideal-coil size ``0.38 √N`` (``chain_coil_size_nm``).
    The last is a placeholder so a G is printed at all; it is *not* a brush
    height and the summary says which one was used.
    """
    from tools.metadata import as_float

    if height_nm:
        return float(height_nm), "given"
    z_max, z_min = as_float(meta, "z_max_nm"), as_float(meta, "z_min_nm")
    if not math.isnan(z_max) and not math.isnan(z_min) and z_max > z_min:
        return z_max - z_min, "trajectory z-extent (z_max_nm - z_min_nm)"
    coil = as_float(meta, "chain_coil_size_nm")
    if not math.isnan(coil) and coil > 0:
        return coil, "ideal coil size 0.38*sqrt(N) (no z-extent in metadata)"
    return None, "unknown"


def analyse_runtime(runtime_dir: Path | str, *, height_nm: float | None = None,
                    events_path: Path | str | None = None) -> tuple[nx.MultiGraph, dict[str, Any]]:
    """Site graph and metrics for a prepared/finished run's runtime/ folder.

    Reads metadata.csv (n_residues, nmol, sequence, box, mode, temperature),
    molecules.fasta as a fallback for the sequence, and crosslink_events.csv
    (absent → a non-reactive run: anchors only). ``mode`` comes from
    metadata when written there, else from surface.yaml, else brush.
    """
    from tools.metadata import as_float, read_metadata
    from tools.surface import mode_of

    runtime_dir = Path(runtime_dir)
    meta = read_metadata(runtime_dir)
    sequence = meta.get("sequence") or ""
    if not sequence:
        fasta = runtime_dir / "molecules.fasta"
        if fasta.is_file():
            sequence = "".join(l.strip() for l in fasta.read_text().splitlines()
                               if l and not l.startswith(">"))
    if not sequence:
        raise FileNotFoundError(f"No sequence in {runtime_dir}/metadata.csv or molecules.fasta")
    n_res = int(as_float(meta, "n_residues", len(sequence)))
    n_chains = int(as_float(meta, "nmol", 0))
    if not n_chains:
        raise ValueError(f"nmol missing from {runtime_dir}/metadata.csv — run `sim metadata` first")
    mode = meta.get("mode") or mode_of(runtime_dir)

    rows = read_events(events_path if events_path else runtime_dir / "crosslink_events.csv")
    G = build_graph(rows, n_chains, n_res, mode=mode, sequence=sequence)

    bx, by = as_float(meta, "box_x_nm"), as_float(meta, "box_y_nm")
    box = (bx, by) if not (math.isnan(bx) or math.isnan(by)) else None
    h, source = coating_height(meta, height_nm)
    T = as_float(meta, "temperature_K", DEFAULT_T_K)
    metrics = topology_metrics(G, box_xy_nm=box, height_nm=h, height_source=source,
                               temperature_K=T if not math.isnan(T) else DEFAULT_T_K)
    metrics["sim_name"] = meta.get("sim_name") or runtime_dir.parent.name
    metrics["events_file"] = str(events_path or runtime_dir / "crosslink_events.csv")
    metrics["n_events_read"] = len(rows)
    return G, metrics


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if math.isnan(value):
            return "—"
        return f"{value:.{digits}g}" if abs(value) < 1e-2 or abs(value) >= 1e4 else f"{value:.{digits}f}"
    return str(value)


def summary_lines(m: Mapping[str, Any]) -> list[str]:
    """Human-readable report of ``topology_metrics`` output."""
    L = []
    L.append(f"network topology — {m.get('sim_name', '')}  ({m['mode']} mode, "
             f"{m['n_chains']} chains x {m['n_residues']} residues, "
             f"{m['n_sites_per_chain']} lysines per chain)")
    L.append("")
    L.append("bonds")
    L.append(f"  crosslinks (K-K):        {m['n_crosslinks']}  "
             f"({m['n_intra']} intra, {m['n_inter']} inter)")
    L.append(f"  surface bonds (K-pin):   {m['n_surface_bonds']}   anchors: {m['n_anchors']}")
    L.append(f"  K-K bonds per chain:     {m['bonds_per_chain_distribution']}   (bonds: chains)")
    L.append(f"  primary loops:           {m['n_primary_loops']} of {m['n_crosslinks']} "
             f"(fraction {_fmt(m['primary_loop_fraction'])}); secondary loops {m['n_secondary_loops']}")
    L.append(f"  primary loop spans:      {m['primary_loop_spans']}")
    L.append("")
    L.append("connectivity")
    L.append(f"  components (with SURFACE): {m['n_components']};  chain clusters via inter bonds: "
             f"{m['n_chain_clusters']}, largest {m['largest_cluster_chains']} chains")
    L.append(f"  chains connected to surface:            {_fmt(m['fraction_chains_connected_to_surface'])}")
    L.append(f"  ... via crosslinks (not own pins):      {_fmt(m['fraction_chains_surface_via_crosslinks'])}")
    L.append(f"  ... with >= 2 independent routes:       {_fmt(m['fraction_chains_two_routes_to_surface'])}")
    L.append(f"  sol chains (not connected to surface):  {m['n_sol_chains']}")
    L.append("")
    L.append("elastic census (SURFACE = the infinite network)")
    L.append(f"  cycle rank (independent cycles):        {m['cycle_rank']}")
    L.append(f"  active strands nu:                      {m['n_strands_active']}   "
             f"({_fmt(m['strands_per_chain'])} per chain, {_fmt(m['strands_per_nm2'])} per nm^2)")
    L.append(f"  active junctions mu:                    {m['n_junctions_active']} of {m['n_junctions']} junctions")
    L.append(f"  functionality f  {{f: junctions}}:        {m['functionality_distribution']}")
    L.append(f"  effective f_eff  {{f_eff: junctions}}:    {m['effective_functionality_distribution']}")
    L.append(f"  pruning count (deg<=1 removed):         {m['n_strands_pruned']} strands, "
             f"{m['n_junctions_pruned']} junctions;  minus primary loops (RENT): {m['n_strands_rent']}")
    L.append(f"  dangling-end contour fraction:          {_fmt(m['dangling_end_fraction'])}")
    L.append(f"  elastically inactive contour fraction:  {_fmt(m['inactive_contour_fraction'])}")
    L.append("")
    L.append("modulus estimates (order of magnitude, for ranking only)")
    L.append(f"  box {_fmt(m['box_x_nm'])} x {_fmt(m['box_y_nm'])} nm, height {_fmt(m['height_nm'])} nm "
             f"[{m['height_source'] or 'none'}] -> V = {_fmt(m['volume_nm3'])} nm^3;  "
             f"kT/nm^3 = {_fmt(m['kT_per_nm3_MPa'])} MPa at {m['temperature_K']} K")
    L.append(f"  G_affine  = nu kT/V:                    {_fmt(m['G_affine_kPa'])} kPa")
    L.append(f"  G_phantom = (nu - mu) kT/V:             {_fmt(m['G_phantom_kPa'])} kPa")
    L.append(f"  G_rent    = (pruned - loops) kT/V:      {_fmt(m['G_rent_kPa'])} kPa")
    L.append("")
    L.append("Topology is meaningful, kinetics are not; every count depends on the reaction")
    L.append("cutoff. See docs/NETWORK.md for definitions, assumptions and limits.")
    return L


COMPARE_COLUMNS = [
    ("run", "sim_name", "s"),
    ("chains", "n_chains", "d"),
    ("xl", "n_crosslinks", "d"),
    ("intra", "n_intra", "d"),
    ("inter", "n_inter", "d"),
    ("surf", "n_surface_bonds", "d"),
    ("loops", "n_primary_loops", "d"),
    ("nu", "n_strands_active", "d"),
    ("mu", "n_junctions_active", "d"),
    ("cyc", "cycle_rank", "d"),
    ("viaXL", "fraction_chains_surface_via_crosslinks", "f"),
    ("2route", "fraction_chains_two_routes_to_surface", "f"),
    ("dangl", "dangling_end_fraction", "f"),
    ("inact", "inactive_contour_fraction", "f"),
    ("nu/nm2", "strands_per_nm2", "g"),
    ("G_aff", "G_affine_kPa", "g"),
    ("G_ph", "G_phantom_kPa", "g"),
]


def compare(runs: Sequence[tuple[str, Mapping[str, Any]]]) -> str:
    """Fixed-width table of the headline metrics for several runs (G in kPa)."""
    def cell(value, kind):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "—"
        if kind == "d":
            return str(int(value))
        if kind == "f":
            return f"{value:.2f}"
        if kind == "g":
            return f"{value:.3g}"
        return str(value)

    table = []
    for name, m in runs:
        row = []
        for header, key, kind in COMPARE_COLUMNS:
            row.append(name if key == "sim_name" else cell(m.get(key), kind))
        table.append(row)
    headers = [h for h, _, _ in COMPARE_COLUMNS]
    widths = [max(len(h), *(len(r[i]) for r in table)) for i, h in enumerate(headers)]
    fmt = "  ".join(f"{{:<{w}}}" if i == 0 else f"{{:>{w}}}" for i, w in enumerate(widths))
    lines = [fmt.format(*headers), fmt.format(*("-" * w for w in widths))]
    lines += [fmt.format(*r) for r in table]
    lines.append("")
    lines.append("nu/mu = active strands/junctions; cyc = cycle rank; viaXL = chains reaching the surface "
                 "through crosslinks; 2route = chains with >= 2 node-disjoint routes to the surface; "
                 "dangl/inact = contour fractions; G in kPa.")
    return "\n".join(lines)


def write_outputs(metrics: Mapping[str, Any], summary_path: Path | None, json_path: Path | None) -> None:
    """network_summary.txt next to the trajectory, network.json in analysis/."""
    if summary_path is not None:
        Path(summary_path).write_text("\n".join(summary_lines(metrics)) + "\n")
    if json_path is not None:
        Path(json_path).write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
