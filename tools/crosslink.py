"""Reactive crosslinking during the run — opt-in, off by default.

Lysine beads that come within a reaction distance of each other are bonded
*while the simulation is running*, so a formed crosslink tethers its chains and
changes what happens next. That is the whole point: detecting contacts after the
fact (as a post-hoc first-passage script does) counts encounters that a real
network would have prevented, because in reality the first bond stops the chains
from wandering off to make the next one.

Off unless asked for. With ``crosslink_distance = None`` — the default — this
module is never imported into the run path, no force is added to the system and
no random numbers are drawn, so a non-reactive run is bit-for-bit what it was
before this feature existed. ``template/run.py`` picks the path:

    settings = load_settings(runtime_dir)
    if settings is None:
        sim.simulate()                      # untouched CALVADOS path
    else:
        run_reactive(sim, settings)

How the bonding works
---------------------
OpenMM can update the parameters of an *existing* bond in a live Context but
cannot add a new one, and rebuilding the Simulation per reaction would cost
seconds and throw away integrator state. So every possible site-site pair is
added up front as a harmonic bond of stiffness zero — exactly zero energy and
zero force, and no effect on the nonbonded interactions, which the beads keep
whether or not they are crosslinked — and reacting a pair means turning its
stiffness up. Switching a stiff bond on instantaneously injects k*(r-r0)^2 of
energy and can spike the temperature to NaN, so the stiffness is ramped in over
``crosslink_ramp_steps``.

The pair list is O(n^2) in the number of reactive sites: 64 sites (16 chains x 4
lysines) is 2016 bonds, and the block-runs batch (16 chains x 16 lysines = 256
sites) is 32640 — a few hundred kB of bond parameters, re-uploaded whole by
``updateParametersInContext`` a handful of times per reaction, which is nothing
next to a 1000-step chunk of MD. Intra-chain pairs closer than ``min_span``
residues are left out (they can never be a meaningful loop). Above roughly 500
sites (125k bonds) this needs replacing with a pair set restricted at setup by
a generous distance cutoff — you get a warning at that point.

What the numbers mean, and what they don't
------------------------------------------
* CALVADOS has no explicit crosslinker, no solvent and no activation barrier.
  The reaction is a pure distance criterion on coarse-grained beads. The
  *sequence and topology* of events is meaningful; the *rate* is not. Do not
  report event times as kinetics.
* CALVADOS2 has no temperature-dependent hydrophobicity, so it cannot reproduce
  the LCST transition that drives coacervation in real ELP crosslinking. Every
  run is effectively at a single solvent quality. Use HPS-T if the transition
  matters to the question.
* ``crosslink_distance`` dominates every result. Sweep it (0.5, 0.8, 1.2 nm) and
  report the sensitivity rather than one number from one cutoff.
* Contacts present in the first frames are artefacts of the initial placement,
  not encounters the dynamics produced. Set ``crosslink_start_step`` past
  equilibration (the analyze notebook's Rg/RMSD check says where that is) and
  state the value used.
* ``crosslink_valence = 1`` models a bifunctional crosslinker (glutaraldehyde at
  1:1). THPP is trifunctional; Kawamoto et al. (Macromolecules 2015, 48, 8980)
  found A2+B4 networks far more prone to cyclic defects than A2+B3, so the
  valence materially changes the loop fraction this reports.

Units are nanometres throughout, matching OpenMM and the box-planning code —
*not* the angstroms an MDAnalysis-based post-hoc script works in. The startup
line prints both so the two cannot be confused.

Surface attachment through lysines (``tools.surface``, the ``free`` and
``preattached`` modes) runs in the same loop (`run_reactive`) and writes to the
same event record: rows with ``kind = surface`` are lysine-surface bonds, and
the ``origin`` column separates bonds present at initialisation from ones formed
during the run. A lysine is in exactly one of three states — free, bonded to the
surface, or bonded to another lysine — so the two reactors mask each other out.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
from yaml import safe_dump, safe_load

# Written by prepare.py next to the trajectory; read by run.py.
CROSSLINK_FILENAME = "crosslink.yaml"
# Outputs, all in the same runtime/ folder as the .dcd.
EVENTS_FILENAME = "crosslink_events.csv"
SUMMARY_FILENAME = "crosslink_summary.txt"
STATE_FILENAME = "crosslink_state.json"

# CALVADOS's Langevin timestep (calvados/sim.py) — fixed, not configurable.
DT_PS = 0.01

# How many parameter updates a stiffness ramp is spread over. The ramp only
# advances when the loop stops to update parameters, so if the loop only ever
# stopped every `check_every` steps a 500-step ramp inside a 1000-step check
# interval would jump straight to full stiffness — a single-step ramp, which is
# what the ramp exists to avoid. While a bond is ramping the loop therefore
# steps in smaller chunks, which changes nothing about the dynamics (stepping
# 100 + 100 is identical to stepping 200) and everything about how gently the
# bond switches on.
RAMP_UPDATES = 5

# Above this many sites the all-pairs bond list stops being free — see the
# module docstring.
PAIR_WARN_SITES = 500

# `kind` is intra / inter for lysine-lysine bonds and surface for a
# lysine-surface bond (then resid_j, chain_j, span and pymol_j are blank and
# x, y, z is the pin). `origin` is "initial" for a bond the run started with
# and "run" for one it formed. Consumers that read by column name are
# unaffected by the extra column; a movie script should skip pymol_j when blank.
#
# `frame` is the 0-based index of the saved frame nearest the event, as mdtraj
# and MDAnalysis number frames. OpenMM's DCDReporter writes its first frame
# after `n_save` steps, not at step 0, so frame k holds step (k + 1) * n_save
# and an event at MD step s is nearest frame round(s / n_save) - 1. A bond
# present at initialisation (step 0) predates every frame and is put on
# frame 0. PyMOL states are 1-based: state = frame + 1 (the pymol_i / pymol_j
# bead indices are already 1-based). `time_ps` is exact (step x 0.01 ps) and is
# what the tooling sorts and filters on; `frame` is only for movies.
CSV_COLUMNS = [
    "frame", "time_ps", "resid_i", "chain_i", "resid_j", "chain_j",
    "kind", "span", "x", "y", "z", "pymol_i", "pymol_j", "origin",
]


class SiteInfo(NamedTuple):
    """One reactive bead, with what the event log needs to describe it."""

    index: int          # bead index in the OpenMM system / mdtraj topology, 0-based
    chain: int          # chain index, 0-based
    resid: int          # residue number within the chain, 1-based (PDB resSeq)
    res_in_chain: int   # residue index within the chain, 0-based (0 = the anchor)
    resname: str
    atomname: str


@dataclass(frozen=True)
class CrosslinkSettings:
    """Everything that defines a reactive run. See the module docstring."""

    distance: float                 # nm — reaction distance between two sites
    valence: int = 1                # max bonds per site
    prob: float = 1.0               # P(react | within cutoff at a check)
    check_every: int = 1000         # MD steps between reaction checks
    k: float = 2000.0               # kJ/mol/nm^2, final bond stiffness
    r0: float = 0.6                 # nm, bond equilibrium length
    ramp_steps: int = 500           # steps over which k rises 0 -> k
    selection: tuple[int, ...] | None = None   # explicit bead indices, or None
    start_step: int = 0             # ignore reactions before this step
    seed: int | None = None         # reaction RNG seed (None -> from the run's)
    # Minimum |residue_i - residue_j| for two sites on the SAME chain to be a
    # candidate pair. Neighbouring CA beads sit one bond length (0.38 nm) apart
    # and i, i+2 at most 0.76 nm, both inside any sensible reaction distance
    # (0.5-1.2 nm), so with min_span <= 2 a sequence such as ...KK... would
    # "crosslink" at the first check for no reason but the backbone, using up
    # valence that a real bifunctional crosslinker could not: it cannot bridge
    # two lysines whose side chains are held side by side into a loop that
    # constrains anything. Pairs closer than min_span along the chain are not
    # even dormant bonds. Inter-chain pairs are never affected. 0 keeps every
    # pair (the behaviour before this setting existed).
    min_span: int = 3

    def __post_init__(self) -> None:
        if self.distance is None or self.distance <= 0:
            raise ValueError("crosslink_distance must be a positive number of nanometres.")
        if self.distance > 5.0:
            # Not fatal — the trivial-reaction test deliberately uses a huge one.
            print(f"⚠  crosslink_distance is {self.distance} nm. That is very large for a "
                  f"bead-bead reaction distance; did you mean angstroms? "
                  f"({self.distance} A = {self.distance / 10} nm)")
        if self.valence < 1:
            raise ValueError("crosslink_valence must be at least 1.")
        if not 0.0 < self.prob <= 1.0:
            raise ValueError("crosslink_prob must be in (0, 1].")
        if self.check_every < 1:
            raise ValueError("crosslink_check_every must be at least 1 step.")
        if self.k <= 0:
            raise ValueError("crosslink_k must be positive.")
        if self.r0 <= 0:
            raise ValueError("crosslink_r0 must be positive.")
        if self.ramp_steps < 0:
            raise ValueError("crosslink_ramp_steps cannot be negative.")
        if self.start_step < 0:
            raise ValueError("crosslink_start_step cannot be negative.")
        if self.min_span < 0:
            raise ValueError("crosslink_min_span cannot be negative.")

    @property
    def distance_angstrom(self) -> float:
        return self.distance * 10.0

    def banner(self) -> str:
        """The startup line. Prints both units so nm/A cannot be mixed up."""
        return (f"crosslinking: ON, reaction distance {self.distance:.2f} nm "
                f"({self.distance_angstrom:.1f} A), valence {self.valence}, p={self.prob}, "
                f"min intra-chain span {self.min_span} residues")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["selection"] = list(self.selection) if self.selection is not None else None
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CrosslinkSettings":
        """Build from a crosslink.yaml / state-file dict. Keys added since the file
        was written (e.g. ``min_span``) fall back to their defaults, unknown keys
        are ignored, so old files still load."""
        known = {f: data[f] for f in cls.__dataclass_fields__ if f in data}
        selection = known.get("selection")
        if selection is not None:
            known["selection"] = tuple(int(i) for i in selection)
        return cls(**known)


def write_settings(runtime_dir: Path, settings: CrosslinkSettings | None) -> Path | None:
    """Record the reactive settings next to the trajectory, or remove a stale file.

    Nothing is written for a non-reactive run: the *absence* of crosslink.yaml is
    what run.py reads as "off", so an old file left behind by a regenerated run
    would silently turn crosslinking on.
    """
    path = Path(runtime_dir) / CROSSLINK_FILENAME
    if settings is None:
        path.unlink(missing_ok=True)
        return None
    with open(path, "w") as f:
        safe_dump(settings.to_dict(), f, sort_keys=True)
    return path


def load_settings(runtime_dir: Path) -> CrosslinkSettings | None:
    """The run's crosslink settings, or None when the run is not reactive."""
    path = Path(runtime_dir) / CROSSLINK_FILENAME
    if not path.is_file():
        return None
    with open(path) as f:
        data = safe_load(f) or {}
    if not data or data.get("distance") is None:
        return None
    return CrosslinkSettings.from_dict(data)


def all_sites(topology) -> list[SiteInfo]:
    """Every bead in the mdtraj topology, described the way the event log needs."""
    records = []
    for chain in topology.chains:
        for res_in_chain, residue in enumerate(chain.residues):
            for atom in residue.atoms:
                records.append(SiteInfo(
                    index=atom.index,
                    chain=chain.index,
                    resid=int(residue.resSeq),
                    res_in_chain=res_in_chain,
                    resname=str(residue.name),
                    atomname=str(atom.name),
                ))
    return records


def find_sites(topology, selection: tuple[int, ...] | None = None,
               verbose: bool = True, drop_anchor: bool = True) -> list[SiteInfo]:
    """The reactive beads: lysines by default, or exactly `selection` if given.

    Auto-detection tries residue name LYS, then residue name K, then bead name K,
    in that order, and stops at the first rule that matches anything. If none
    match it raises with the names actually present rather than falling back to
    something arbitrary — silently crosslinking the wrong beads is worse than
    not starting.

    In brush mode (`drop_anchor=True`, the default) residue 0 of every chain is
    dropped whatever the rule says. That bead is the surface anchor: `sim new`
    tags it "Z", CALVADOS gives it zero effective mass and OpenMM holds it fixed,
    so bonding to it would tether a chain to an immovable point rather than to
    another chain. In the free/preattached modes residue 0 is an ordinary
    residue and is kept (`drop_anchor=False`).
    """
    records = all_sites(topology)

    if selection is not None:
        by_index = {r.index: r for r in records}
        missing = sorted(set(selection) - by_index.keys())
        if missing:
            raise ValueError(
                f"crosslink_selection contains bead indices that are not in the system: "
                f"{missing[:10]}{'...' if len(missing) > 10 else ''} "
                f"(the system has {len(records)} beads, 0-based)."
            )
        chosen = [by_index[i] for i in sorted(set(selection))]
        rule = f"explicit crosslink_selection ({len(chosen)} beads)"
    else:
        rules = (
            ("residue name LYS", lambda r: r.resname.upper() == "LYS"),
            ("residue name K", lambda r: r.resname.upper() == "K"),
            ("bead name K", lambda r: r.atomname.upper() == "K"),
        )
        chosen, rule = [], ""
        for label, match in rules:
            chosen = [r for r in records if match(r)]
            if chosen:
                rule = label
                break
        if not chosen:
            resnames = sorted({r.resname for r in records})
            atomnames = sorted({r.atomname for r in records})
            raise ValueError(
                "No lysine beads found to crosslink. Tried residue name LYS, residue name K "
                f"and bead name K.\n  residue names present: {', '.join(resnames)}\n"
                f"  bead names present: {', '.join(atomnames)}\n"
                "Pass crosslink_selection with explicit bead indices if the sites are named "
                "something else."
            )

    anchors = [r for r in chosen if r.res_in_chain == 0] if drop_anchor else []
    sites = [r for r in chosen if r.res_in_chain != 0] if drop_anchor else list(chosen)

    if verbose:
        print(f"crosslinking: {len(sites)} reactive sites, matched by {rule}")
        if anchors:
            print(f"crosslinking: dropped {len(anchors)} site(s) at residue 0 of their chain — "
                  f"those are the fixed surface anchors, not free lysines")
        if len(sites) > PAIR_WARN_SITES:
            n_pairs = len(sites) * (len(sites) - 1) // 2
            print(f"⚠  {len(sites)} reactive sites means {n_pairs} dormant bonds. This all-pairs "
                  f"setup is only cheap below ~{PAIR_WARN_SITES} sites; expect slow setup and "
                  f"high memory. Restrict the pair set by a setup-time distance cutoff instead.")
    if not sites:
        raise ValueError(
            "Every candidate reactive site is a fixed surface anchor (residue 0 of its chain), "
            "so there is nothing to crosslink."
        )
    return sites


def min_image_distances(pos: np.ndarray, idx_i: np.ndarray, idx_j: np.ndarray,
                        box_lengths: np.ndarray) -> np.ndarray:
    """Distances between bead pairs under the minimum-image convention, in nm.

    The box CALVADOS builds is periodic in all three dimensions — z included; the
    surface is a potential well at low z, not a hard boundary — so minimum-
    imaging z is consistent with the system as declared, not an approximation.
    In practice chains never come within L_z/2 of each other in z anyway.
    """
    d = pos[idx_i] - pos[idx_j]
    d -= box_lengths * np.round(d / box_lengths)
    return np.linalg.norm(d, axis=1)


def frame_of_step(step: int, n_save: int) -> int:
    """0-based index of the saved frame nearest MD step `step`.

    DCDReporter writes frame k at step (k + 1) * n_save (its first report comes
    after n_save steps, not at step 0), so this is round(step / n_save) - 1,
    floored at 0 for events that predate the first frame.
    """
    if not n_save:
        return 0
    return max(0, int(round(step / n_save)) - 1)


@dataclass
class CrosslinkEvent:
    step: int
    site_a: int
    site_b: int
    distance: float
    midpoint: tuple[float, float, float]

    def to_dict(self) -> dict[str, Any]:
        return {"step": self.step, "site_a": self.site_a, "site_b": self.site_b,
                "distance": self.distance, "midpoint": list(self.midpoint)}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "CrosslinkEvent":
        return cls(int(d["step"]), int(d["site_a"]), int(d["site_b"]),
                   float(d["distance"]), tuple(float(v) for v in d["midpoint"]))


class Crosslinker:
    """The reacting network: which pairs may bond, which have, and how stiff they are.

    Holds one HarmonicBondForce whose bonds are all the site-site pairs. A bond
    at stiffness 0 is dormant (no energy, no force); reacting a pair ramps its
    stiffness up to `settings.k`.
    """

    def __init__(self, settings: CrosslinkSettings, sites: list[SiteInfo], seed: int):
        self.settings = settings
        self.sites = sites
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        n = len(sites)
        # Upper-triangle pair list, in a fixed order: pairs[b] is the pair
        # belonging to bond index b, and that order has to be reproducible
        # because a restart rebuilds the force from scratch and re-applies the
        # saved reacted flags by bond index. Pairs on one chain closer than
        # `min_span` residues are dropped here, so they are not even dormant
        # bonds (see CrosslinkSettings.min_span); save_state records a hash of
        # the resulting list so a state written with a different one is refused.
        slot_a, slot_b = np.triu_indices(n, k=1)
        chain = np.array([s.chain for s in sites], dtype=int)
        res = np.array([s.res_in_chain for s in sites], dtype=int)
        too_close = (chain[slot_a] == chain[slot_b]) & \
            (np.abs(res[slot_a] - res[slot_b]) < settings.min_span)
        self.n_pairs_dropped_by_span = int(too_close.sum())
        slot_a, slot_b = slot_a[~too_close], slot_b[~too_close]
        self.pair_slot_a = slot_a
        self.pair_slot_b = slot_b
        bead = np.array([s.index for s in sites], dtype=int)
        self.site_beads = bead
        self.pair_bead_i = bead[slot_a]
        self.pair_bead_j = bead[slot_b]

        self.reacted = np.zeros(len(slot_a), dtype=bool)
        self.used = np.zeros(n, dtype=int)
        self.ramping: dict[int, int] = {}      # bond index -> step it formed
        self.events: list[CrosslinkEvent] = []
        self.force = None
        self.resume_step = 0                   # set by load_state on a restart

    # ---- setup ------------------------------------------------------------

    @property
    def n_pairs(self) -> int:
        return len(self.pair_bead_i)

    @property
    def pairs_hash(self) -> str:
        """Fingerprint of the bond list (bead i, bead j per bond index). Two
        Crosslinkers with the same hash number their bonds identically."""
        import hashlib

        h = hashlib.sha1()
        h.update(np.ascontiguousarray(self.pair_bead_i, dtype=np.int64).tobytes())
        h.update(np.ascontiguousarray(self.pair_bead_j, dtype=np.int64).tobytes())
        return h.hexdigest()

    def add_force(self, system) -> Any:
        """Add every pair as a dormant (zero-stiffness) harmonic bond.

        Called before the Context exists, so bonds restored from a checkpointed
        state are already at the right stiffness when the Context is built — no
        update-in-context needed, and nothing switches on abruptly.
        """
        from openmm import HarmonicBondForce
        from openmm.unit import kilojoule_per_mole, nanometer

        force = HarmonicBondForce()
        force.setName("crosslinks")
        # THE critical line. The reaction criterion is a *minimum-image*
        # distance (min_image_distances), so two lysines sitting 0.6 nm apart
        # across a periodic boundary are a legitimate pair and do get bonded.
        # Without this flag OpenMM evaluates that bond on the raw coordinate
        # difference instead — one box length, ~28 nm in the block runs — and a
        # k = 2000 kJ/mol/nm^2 bond stretched 27 nm past its 0.6 nm rest length
        # pulls with 5e4 kJ/mol/nm. The system is destroyed within a few
        # thousand steps.
        #
        # This is not hypothetical: it wrecked all 12 crosslinking runs of the
        # 2026-09-19 block batch (see docs/BUGS.md BUG-XL-7). Every force
        # CALVADOS builds itself sets this (calvados/interactions.py), and the
        # surface tethers get it via periodicdistance() in their expression;
        # this force was the one that did not.
        force.setUsesPeriodicBoundaryConditions(True)
        r0 = self.settings.r0 * nanometer
        unit_k = kilojoule_per_mole / nanometer ** 2
        for b in range(self.n_pairs):
            # A bond the checkpoint caught mid-ramp comes back at the stiffness
            # it had reached, not at full: handing it the final k here would be
            # the instantaneous switch-on the ramp exists to avoid.
            k = self.settings.k * self._ramp_fraction(b, self.resume_step) if self.reacted[b] else 0.0
            force.addBond(int(self.pair_bead_i[b]), int(self.pair_bead_j[b]), r0, k * unit_k)
        system.addForce(force)
        self.force = force
        return force

    # ---- the reaction -----------------------------------------------------

    def _ramp_fraction(self, bond: int, step: int) -> float:
        """How much of its final stiffness a reacted bond has earned by `step`."""
        formed_at = self.ramping.get(int(bond))
        if formed_at is None or self.settings.ramp_steps <= 0:
            return 1.0
        return min(1.0, max(0.0, (step - formed_at) / self.settings.ramp_steps))

    def _set_bond(self, bond: int, k: float) -> None:
        """Set bond `bond` to stiffness `k` (kJ/mol/nm^2) in the force object; the
        Context sees it after updateParametersInContext."""
        from openmm.unit import kilojoule_per_mole, nanometer

        self.force.setBondParameters(int(bond), int(self.pair_bead_i[bond]),
                                     int(self.pair_bead_j[bond]),
                                     self.settings.r0 * nanometer,
                                     k * kilojoule_per_mole / nanometer ** 2)

    def advance_ramps(self, step: int) -> bool:
        """Push every ramping bond to the stiffness its age entitles it to."""
        if not self.ramping:
            return False
        done = []
        for b in list(self.ramping):
            frac = self._ramp_fraction(b, step)
            self._set_bond(int(b), frac * self.settings.k)
            if frac >= 1.0:
                done.append(b)
        for b in done:
            del self.ramping[b]
        return True

    def check_reactions(self, context, step: int, blocked: np.ndarray | None = None) -> bool:
        """React every pair within the cutoff, reading positions from a live Context."""
        from openmm.unit import nanometer

        # enforcePeriodicBox=False keeps whole chains intact; the minimum image
        # is applied per pair below, which is what the distance criterion needs.
        state = context.getState(getPositions=True, enforcePeriodicBox=False)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
        return self.react(pos, np.diag(np.asarray(box)), step, blocked)   # orthorhombic box

    def react(self, pos: np.ndarray, box_lengths: np.ndarray, step: int,
              blocked: np.ndarray | None = None) -> bool:
        """React every pair that is within the cutoff, closest pairs first.

        Pure geometry, so the identical rule can be replayed over the frames of a
        finished trajectory (`post_hoc_events`) and compared against what the
        reactive run actually did. `blocked` marks sites (by slot) that are
        bonded to the surface and so can never crosslink.
        """
        settings = self.settings

        spare = self.used < settings.valence
        if blocked is not None:
            spare = spare & ~np.asarray(blocked, dtype=bool)
        live = (~self.reacted) & spare[self.pair_slot_a] & spare[self.pair_slot_b]
        candidates = np.flatnonzero(live)
        if candidates.size == 0:
            return False

        d = min_image_distances(pos, self.pair_bead_i[candidates],
                                self.pair_bead_j[candidates], box_lengths)
        within = d < settings.distance
        if not within.any():
            return False

        # Nearest first, so a lysine pairs with its closest available partner
        # rather than whichever one happens to come first in the pair list.
        order = np.argsort(d[within], kind="stable")
        bonds = candidates[within][order]
        dists = d[within][order]

        changed = False
        for b, dist in zip(bonds, dists):
            a_slot, b_slot = self.pair_slot_a[b], self.pair_slot_b[b]
            # Re-check: a site may have used up its valence earlier in this pass.
            if self.used[a_slot] >= settings.valence or self.used[b_slot] >= settings.valence:
                continue
            if settings.prob < 1.0 and self.rng.random() > settings.prob:
                continue
            self.reacted[b] = True
            self.used[a_slot] += 1
            self.used[b_slot] += 1
            if self.force is not None:
                # Post-hoc replay has no force to ramp; only a live run does.
                # Write the bond's starting stiffness now rather than waiting
                # for the next advance_ramps: with ramp_steps > 0 that is 0
                # (unchanged), but with ramp_steps = 0 — "switch on at once" —
                # it is the full k, and leaving it to the next chunk would
                # silently delay the bond by up to check_every steps.
                self.ramping[int(b)] = step
                self._set_bond(int(b), self.settings.k * self._ramp_fraction(int(b), step))
            i, j = self.pair_bead_i[b], self.pair_bead_j[b]
            # Midpoint of the minimum-image pair, not of the raw coordinates: the
            # positions come unwrapped (enforcePeriodicBox=False) and a pair that
            # met across the boundary would otherwise be marked mid-box.
            dij = pos[i] - pos[j]
            dij -= box_lengths * np.round(dij / box_lengths)
            midpoint = pos[j] + 0.5 * dij
            self.events.append(CrosslinkEvent(
                step=step, site_a=int(a_slot), site_b=int(b_slot),
                distance=float(dist),
                midpoint=(float(midpoint[0]), float(midpoint[1]), float(midpoint[2])),
            ))
            changed = True
        return changed

    # ---- restart safety ---------------------------------------------------

    def save_state(self, runtime_dir: Path, step: int) -> Path:
        """Persist the network next to the checkpoint.

        Without this a restarted run would rebuild every bond dormant and quietly
        produce a different network from the one it had already formed.
        """
        path = Path(runtime_dir) / STATE_FILENAME
        payload = {
            "step": int(step),
            "seed": int(self.seed),
            "settings": self.settings.to_dict(),
            "site_beads": [int(i) for i in self.site_beads],
            # The bond list itself, fingerprinted: `reacted` is a list of bond
            # indices and means nothing under a different pair list.
            "n_pairs": int(self.n_pairs),
            "pairs_hash": self.pairs_hash,
            "reacted": np.flatnonzero(self.reacted).tolist(),
            "used": self.used.tolist(),
            "ramping": {str(int(b)): int(t) for b, t in self.ramping.items()},
            "events": [e.to_dict() for e in self.events],
        }
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=1))
        tmp.replace(path)   # atomic: a crash mid-write must not eat the network
        return path

    def load_state(self, runtime_dir: Path) -> int | None:
        """Restore a persisted network. Returns the step it was saved at, or None.

        The saved sites must be the same beads in the same order, or the bond
        indices mean something different and the network would be reconstructed
        wrong; that is a hard error rather than a warning.
        """
        path = Path(runtime_dir) / STATE_FILENAME
        if not path.is_file():
            return None
        payload = json.loads(path.read_text())

        saved_beads = [int(i) for i in payload.get("site_beads", [])]
        if saved_beads != [int(i) for i in self.site_beads]:
            raise ValueError(
                f"{STATE_FILENAME} was written for a different set of reactive sites "
                f"({len(saved_beads)} beads) than this system has ({len(self.site_beads)}). "
                "Delete the runtime folder and start over rather than restarting into a "
                "network that means something else."
            )
        # Same sites is necessary but not sufficient: the bond indices in
        # `reacted` also depend on which pairs were kept (min_span). A state
        # file from before the hash existed was written with every pair kept,
        # so it is compatible exactly when this Crosslinker dropped nothing.
        saved_hash = payload.get("pairs_hash")
        n_all = len(self.site_beads) * (len(self.site_beads) - 1) // 2
        if saved_hash is not None:
            same_pairs = saved_hash == self.pairs_hash
        else:
            same_pairs = int(payload.get("n_pairs", n_all)) == self.n_pairs == n_all
        if not same_pairs:
            saved_span = payload.get("settings", {}).get("min_span", 0)
            raise ValueError(
                f"{STATE_FILENAME} was written for a different bond list "
                f"({payload.get('n_pairs', n_all)} pairs, min_span {saved_span}) than this run "
                f"builds ({self.n_pairs} pairs, min_span {self.settings.min_span}), so its bond "
                "indices would be reapplied to the wrong pairs. Restart with the same "
                "crosslink_min_span, or delete the runtime folder and start over.")
        saved_settings = payload.get("settings", {})
        for key in ("distance", "valence", "r0", "k"):
            if key in saved_settings and saved_settings[key] != getattr(self.settings, key):
                print(f"⚠  crosslink_{key} changed since the checkpoint "
                      f"({saved_settings[key]} -> {getattr(self.settings, key)}); bonds already "
                      f"formed keep the geometry they were formed with.")

        self.reacted[:] = False
        self.reacted[np.array(payload["reacted"], dtype=int)] = True
        self.used = np.array(payload["used"], dtype=int)
        self.events = [CrosslinkEvent.from_dict(e) for e in payload.get("events", [])]
        # `used` is derivable from `reacted`; if the two disagree the file is
        # corrupt (or hand-edited) and the valence bookkeeping would be wrong
        # for the rest of the run. Refuse rather than carry on.
        counts = np.bincount(self.pair_slot_a[self.reacted], minlength=len(self.used)) + \
            np.bincount(self.pair_slot_b[self.reacted], minlength=len(self.used))
        if len(self.used) != len(self.site_beads) or not np.array_equal(counts, self.used):
            raise ValueError(f"{STATE_FILENAME} is inconsistent: per-site bond counts from "
                             f"'reacted' do not match 'used'. The file is corrupt; do not restart "
                             "from it.")
        if len(self.events) != int(self.reacted.sum()):
            raise ValueError(f"{STATE_FILENAME} is inconsistent: {len(self.events)} events but "
                             f"{int(self.reacted.sum())} reacted bonds.")
        step = int(payload.get("step", 0))
        # A bond caught mid-ramp resumes from the stiffness it had reached:
        # `ramping` stores the step it formed at, and _ramp_fraction(b, step)
        # measured from there is the fraction it had reached, so the value is
        # kept as is (clamped to the save step in case of a corrupt file).
        self.ramping = {int(b): min(int(formed_at), step)
                        for b, formed_at in payload.get("ramping", {}).items()}
        # Re-seed from (seed, step) rather than replaying the draws the
        # pre-checkpoint reactions consumed: a plain default_rng(seed) here would
        # hand the restarted run the identical random stream it already used.
        # Deterministic given the checkpoint, which is what reproducibility needs.
        self.rng = np.random.default_rng([self.seed, step])
        self.resume_step = step
        return step

    # ---- outputs ----------------------------------------------------------

    def event_rows(self, n_save: int) -> list[dict[str, Any]]:
        """The event log, one row per crosslink, in `CSV_COLUMNS` order."""
        rows = []
        for event in self.events:
            a, b = self.sites[event.site_a], self.sites[event.site_b]
            intra = a.chain == b.chain
            # The frame a movie should flash on, not the MD step: reactions are
            # checked more often than frames are saved, so round to the nearest
            # saved frame — 0-based, and frame k holds step (k + 1) * n_save
            # (see the CSV_COLUMNS comment).
            frame = frame_of_step(event.step, n_save)
            x, y, z = event.midpoint
            rows.append({
                "frame": frame,
                "time_ps": round(event.step * DT_PS, 6),
                "resid_i": a.resid,
                "chain_i": a.chain,
                "resid_j": b.resid,
                "chain_j": b.chain,
                "kind": "intra" if intra else "inter",
                # Loop size along the backbone. Only defined within one chain;
                # blank (not 0) between chains, which is a different thing.
                "span": abs(a.resid - b.resid) if intra else "",
                "x": round(x, 4),
                "y": round(y, 4),
                "z": round(z, 4),
                # PyMOL's `index` selector is 1-based; OpenMM and MDAnalysis are 0-based.
                "pymol_i": a.index + 1,
                "pymol_j": b.index + 1,
                # Lysine-lysine bonds are only ever formed during the run.
                "origin": "run",
            })
        return rows

    def write_events(self, runtime_dir: Path, n_save: int) -> Path:
        """crosslink_events.csv, with the column names the post-hoc tooling expects."""
        path = Path(runtime_dir) / EVENTS_FILENAME
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(self.event_rows(n_save))
        return path

    def summary_lines(self, n_save: int, steps_done: int,
                      pinned: np.ndarray | None = None) -> list[str]:
        """The summary text. `pinned` (per slot, from the SurfaceBinder) marks
        lysines bonded to the surface, which can never crosslink; the conversion
        is then also given against the sites that were actually available."""
        rows = self.event_rows(n_save)
        n = len(rows)
        intra = [r for r in rows if r["kind"] == "intra"]
        inter = [r for r in rows if r["kind"] == "inter"]
        spans = sorted(int(r["span"]) for r in intra)
        n_sites = len(self.sites)
        capacity = n_sites * self.settings.valence
        n_pinned = int(np.sum(pinned)) if pinned is not None else 0
        n_avail = n_sites - n_pinned
        capacity_avail = n_avail * self.settings.valence

        lines = [
            f"crosslink summary — written {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
            "",
            f"steps simulated:       {steps_done}",
            f"reactive sites:        {n_sites}",
            f"candidate pairs:       {self.n_pairs}"
            + (f"  ({self.n_pairs_dropped_by_span} intra-chain pairs closer than "
               f"{self.settings.min_span} residues excluded)" if self.n_pairs_dropped_by_span else ""),
            f"crosslinks formed:     {n}",
            # Conversion is per *site*: each bond consumes one valence slot at
            # each end, so 2 x events out of n_sites x valence available slots.
            f"conversion:            {2 * n / capacity:.4f}  (2 x events / (sites x valence))"
            if capacity else "conversion:            n/a",
            f"saturated sites:       {int(np.sum(self.used >= self.settings.valence))} of {n_sites}",
        ]
        if pinned is not None:
            # A surface-bonded lysine is out of the pool for good, so the
            # all-sites conversion above understates how far the *available*
            # sites got. Both are reported; say which one you quote.
            lines += [
                f"surface-bonded sites:  {n_pinned} of {n_sites} (cannot crosslink; end of run)",
                f"conversion (available): {2 * n / capacity_avail:.4f}  "
                f"(2 x events / ((sites - surface-bonded) x valence))"
                if capacity_avail else "conversion (available): n/a",
                f"saturated (available): "
                f"{int(np.sum((self.used >= self.settings.valence) & ~np.asarray(pinned, dtype=bool)))} "
                f"of {n_avail}",
            ]
        lines += [
            "",
            f"intra-chain (loops):   {len(intra)}"
            + (f"  ({len(intra) / n:.3f} of all events)" if n else ""),
            f"inter-chain (bridges): {len(inter)}"
            + (f"  ({len(inter) / n:.3f} of all events)" if n else ""),
        ]
        if spans:
            median = spans[len(spans) // 2] if len(spans) % 2 else \
                (spans[len(spans) // 2 - 1] + spans[len(spans) // 2]) / 2
            lines.append(f"loop span (residues):  min {spans[0]}, median {median}, max {spans[-1]}")
        else:
            lines.append("loop span (residues):  no intra-chain events")

        lines += [
            "",
            "The intra fraction is the primary output: it is the primary-loop estimate that",
            "feeds real elastic network theory (RENT). Event *times* are not kinetics — see",
            "the caveats in tools/crosslink.py.",
            "",
            "settings",
            "--------",
            f"distance:        {self.settings.distance} nm ({self.settings.distance_angstrom} A)",
            f"valence:         {self.settings.valence}",
            f"prob:            {self.settings.prob}",
            f"check_every:     {self.settings.check_every} steps "
            f"({self.settings.check_every * DT_PS} ps)",
            f"start_step:      {self.settings.start_step}",
            f"k:               {self.settings.k} kJ/mol/nm^2",
            f"r0:              {self.settings.r0} nm",
            f"ramp_steps:      {self.settings.ramp_steps}",
            f"min_span:        {self.settings.min_span} residues (intra-chain)",
            f"selection:       {'auto (lysines)' if self.settings.selection is None else 'explicit'}",
            f"reaction seed:   {self.seed}",
        ]
        return lines

    def write_summary(self, runtime_dir: Path, n_save: int, steps_done: int) -> Path:
        path = Path(runtime_dir) / SUMMARY_FILENAME
        path.write_text("\n".join(self.summary_lines(n_save, steps_done)) + "\n")
        return path

    def write_outputs(self, runtime_dir: Path, n_save: int, steps_done: int) -> None:
        """Both output files. Called at every checkpoint, so a killed run keeps its log."""
        self.write_events(runtime_dir, n_save)
        self.write_summary(runtime_dir, n_save, steps_done)


def _chunk_sizes(remaining: int, check_every: int, ramp_steps: int, ramping: bool,
                 until: "tuple[int, ...] | list[int]" = ()) -> int:
    """How far to step before stopping to update bond parameters.

    Normally one check interval. While a bond is mid-ramp, smaller steps, so the
    ramp is resolved over several parameter updates instead of jumping to full
    stiffness at the next check (see RAMP_UPDATES). `until` lists the steps
    left to the next scheduled stops (next reaction check per reactor, next
    checkpoint): the chunk never runs past any of them, otherwise a check due
    at step N would be made at the first stop after N — up to check_every - 1
    steps late, and every later check would inherit the delay.
    """
    chunk = check_every
    if ramping and ramp_steps > 0:
        chunk = min(chunk, max(1, ramp_steps // RAMP_UPDATES))
    chunk = min(chunk, remaining, *until)
    return max(1, chunk)


def _reaction_seed(settings: CrosslinkSettings, sim) -> int:
    """Seed for the reaction RNG: the setting, else the run's, else a fresh draw.

    Always ends up a concrete number and is recorded in the summary, so which
    pairs reacted is reproducible even when OpenMM's own seed was left open.
    """
    if settings.seed is not None:
        return int(settings.seed)
    run_seed = getattr(sim, "random_number_seed", None)
    if run_seed is not None:
        return int(run_seed)
    return int(np.random.SeedSequence().entropy % (2 ** 63))


def write_event_record(runtime_dir: Path, n_save: int, steps_done: int,
                       xl: "Crosslinker | None", sb=None, n_chains: int | None = None) -> None:
    """crosslink_events.csv and crosslink_summary.txt for whichever reactors ran.

    One record for both: surface rows (kind = surface) and lysine-lysine rows
    share the schema; bonds present at t = 0 come first, the rest in step order.
    """
    runtime_dir = Path(runtime_dir)
    rows: list[dict[str, Any]] = []
    if sb is not None:
        rows += sb.event_rows(n_save)
    if xl is not None:
        rows += xl.event_rows(n_save)
    rows.sort(key=lambda r: (r["origin"] != "initial", float(r["time_ps"])))
    with open(runtime_dir / EVENTS_FILENAME, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    lines: list[str] = []
    if xl is not None:
        lines += xl.summary_lines(n_save, steps_done,
                                  pinned=sb.pinned if sb is not None else None)
    else:
        lines += [f"crosslink summary — written "
                  f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}",
                  "", f"steps simulated:       {steps_done}",
                  "lysine-lysine crosslinking: off (crosslink_distance = None)"]
    if sb is not None:
        lines += ["", *sb.summary_lines(n_chains if n_chains is not None else 0)]
    (runtime_dir / SUMMARY_FILENAME).write_text("\n".join(lines) + "\n")


def run_reactive(sim, settings: CrosslinkSettings | None, runtime_dir: Path | None = None,
                 surface=None) -> "Crosslinker | None":
    """Run `sim` with reactive crosslinking and/or surface attachment. Mirrors
    calvados Sim.simulate().

    Same integrator, platform, reporters, checkpointing and final PDBs as the
    stock path — the only difference is that stepping is broken into chunks so
    bond parameters can be updated between them. Chunking does not change the
    dynamics (stepping 500 + 500 is stepping 1000) and reporters fire on the
    simulation's own step counter, so the trajectory is unaffected by it.

    `settings` is the lysine-lysine crosslinking (None = off); `surface` is a
    tools.surface.SurfaceSettings for the free/preattached modes (None = brush).
    With `surface=None` this is exactly the loop it was before surface modes
    existed. At least one of the two must be given.

    Only the configuration `sim new` actually generates is supported; the
    equilibration and clock-time modes are refused rather than half-replicated.
    """
    import os

    import openmm
    from openmm import app, unit
    from tqdm import tqdm

    if settings is None and surface is None:
        raise ValueError("run_reactive needs crosslink settings, surface settings or both.")

    path = str(sim.path)
    runtime_dir = Path(runtime_dir) if runtime_dir is not None else Path(path)

    for flag in ("slab_eq", "box_eq", "bilayer_eq"):
        if getattr(sim, flag, False):
            raise NotImplementedError(
                f"Reactive crosslinking does not support config '{flag}'. That mode rebuilds "
                f"the Simulation mid-run, which would drop the crosslink force."
            )
    if getattr(sim, "runtime", 0) > 0:
        raise NotImplementedError(
            "Reactive crosslinking needs a step count, not the clock-time mode "
            f"(config runtime = {sim.runtime} h). Set runtime: 0 and use steps."
        )

    brush = surface is None
    n_chains = int(sim.top.n_chains)

    xl = None
    if settings is not None:
        print(settings.banner(), flush=True)
        sites = find_sites(sim.top, settings.selection, drop_anchor=brush)
        seed = _reaction_seed(settings, sim)
        xl = Crosslinker(settings, sites, seed)

    sb = None
    if surface is not None:
        from tools.surface import STATE_FILENAME as SURFACE_STATE, SurfaceBinder

        print(surface.banner(), flush=True)
        # The two reactors share one site list (all lysines, residue 0 included)
        # so that "this lysine is bonded to the surface" and "this lysine is
        # crosslinked" index the same slots and can exclude each other.
        surface_sites = find_sites(sim.top, None, verbose=xl is None, drop_anchor=False)
        if xl is not None and settings.selection is None and \
                [s.index for s in xl.sites] != [s.index for s in surface_sites]:
            raise ValueError("the crosslink and surface site lists differ; this should not happen "
                             "when crosslink_selection is None.")
        if xl is not None and settings.selection is not None:
            raise NotImplementedError("crosslink_selection is not supported together with a "
                                      "surface mode: both reactors need the same lysine list.")
        sb = SurfaceBinder(surface, surface_sites, int(surface.seed))

    fcheck_in = f"{path}/{sim.frestart}"
    fcheck_out = f"{path}/restart.chk"
    restarting = os.path.isfile(fcheck_in) and sim.restart == "checkpoint"
    append = False

    if restarting:
        # Restoring a network is the whole reason crosslink_state.json exists: a
        # restart that rebuilt every bond dormant would silently make a
        # different network from the one the run had already formed.
        if xl is not None:
            saved_at = xl.load_state(runtime_dir)
            if saved_at is None:
                raise FileNotFoundError(
                    f"Restarting from {fcheck_in} but no {STATE_FILENAME} was found next to it. "
                    f"The crosslinks formed before the checkpoint would be lost and the run would "
                    f"continue with a different network. Delete the checkpoint to start over, or "
                    f"restore the state file."
                )
            print(f"crosslinking: resuming with {int(xl.reacted.sum())} bonds already formed "
                  f"(state saved at step {saved_at})")
        if sb is not None:
            saved_at = sb.load_state(runtime_dir)
            if saved_at is None:
                raise FileNotFoundError(
                    f"Restarting from {fcheck_in} but no {SURFACE_STATE} was found next to it. "
                    f"The surface bonds would be forgotten and the chains would float away. "
                    f"Delete the checkpoint to start over, or restore the state file."
                )
            print(f"surface: resuming with {sb.n_pinned} lysines bonded to the surface "
                  f"(state saved at step {saved_at})")
    else:
        # A stale state file from a previous attempt must not resurrect bonds
        # into a run that is starting from scratch.
        if xl is not None:
            stale = runtime_dir / STATE_FILENAME
            if stale.is_file():
                print(f"crosslinking: ignoring stale {STATE_FILENAME} (not restarting from a checkpoint)")
                stale.unlink()
        if sb is not None:
            stale = runtime_dir / SURFACE_STATE
            if stale.is_file():
                print(f"surface: ignoring stale {SURFACE_STATE} (not restarting from a checkpoint)")
                stale.unlink()
            plan = getattr(sim, "surface_plan", None)
            if plan is None:
                raise RuntimeError("no initial surface plan on the Sim object — prepare.py's "
                                   "build_sim must call tools.surface.place_chains for "
                                   "free/preattached modes.")
            sb.apply_plan(plan, n_chains)

    if xl is not None:
        xl.add_force(sim.system)
        print(f"crosslinking: {xl.n_pairs} dormant bonds added, checks every "
              f"{settings.check_every} steps ({settings.check_every * DT_PS:.1f} ps)"
              + (f", ignored before step {settings.start_step}" if settings.start_step else ""),
              flush=True)
    if sb is not None:
        sb.add_force(sim.system)
        print(f"surface: {sb.n_sites} tethers added ({sb.n_pinned} active), cap {sb.n_max} lysines, "
              f"checks every {surface.check_every} steps ({surface.check_every * DT_PS:.1f} ps)"
              + (f", ignored before step {surface.start_step}" if surface.start_step else ""),
              flush=True)

    # ---- everything below mirrors calvados.sim.Sim.simulate() --------------
    if sim.restart == "pdb" and os.path.isfile(fcheck_in):
        pdb = app.pdbfile.PDBFile(fcheck_in)
    else:
        pdb = app.pdbfile.PDBFile(sim.pdb_cg)

    integrator = openmm.openmm.LangevinMiddleIntegrator(
        sim.temp * unit.kelvin, sim.friction_coeff / unit.picosecond, DT_PS * unit.picosecond)
    if sim.random_number_seed is not None:
        integrator.setRandomNumberSeed(sim.random_number_seed)
    print(integrator.getFriction(), integrator.getTemperature())

    platform = openmm.Platform.getPlatformByName(sim.platform)
    if sim.platform == "CPU":
        simulation = app.simulation.Simulation(pdb.topology, sim.system, integrator, platform,
                                               dict(Threads=str(sim.threads)))
    else:
        if os.environ.get("CUDA_VISIBLE_DEVICES") is None:
            platform.setPropertyDefaultValue("DeviceIndex", str(sim.gpu_id))
        simulation = app.simulation.Simulation(pdb.topology, sim.system, integrator, platform)
    print("Running on", platform.getName())

    if restarting:
        if not os.path.isfile(f"{path}/{sim.sysname:s}.dcd"):
            raise Exception(f"Did not find {path}/{sim.sysname:s}.dcd trajectory to append to!")
        append = True
        print(f"Reading check point file {fcheck_in}")
        print(f"Appending trajectory to {path}/{sim.sysname:s}.dcd")
        simulation.loadCheckpoint(fcheck_in)
        # The checkpoint and the network state(s) are written one after the
        # other; a kill in between leaves coordinates from one step and bonds
        # from another. currentStep comes from the checkpoint itself, so the
        # three can be compared — and a mismatch is a wrong network, not a
        # warning.
        for name, reactor in (("crosslink", xl), ("surface", sb)):
            if reactor is not None and reactor.resume_step != simulation.currentStep:
                raise RuntimeError(
                    f"{name} state was saved at step {reactor.resume_step} but the checkpoint "
                    f"{fcheck_in} is at step {simulation.currentStep}. They do not belong "
                    "together; restore a matching pair or start over.")
    else:
        if sim.restart == "pdb":
            print(f"Reading in system configuration {sim.frestart}")
        elif sim.restart == "checkpoint":
            print(f"No checkpoint file {sim.frestart} found: Starting from new system configuration")
        else:
            print("Starting from new system configuration")

        if os.path.isfile(f"{path}/{sim.sysname:s}.dcd"):
            stamp = datetime.now().strftime("%Y%d%m_%Hh%Mm%Ss")
            print(f"Backing up existing {path}/{sim.sysname:s}.dcd to "
                  f"{path}/backup_{sim.sysname:s}_{stamp}.dcd")
            os.rename(f"{path}/{sim.sysname:s}.dcd", f"{path}/backup_{sim.sysname:s}_{stamp}.dcd")
        print(f"Writing trajectory to new file {path}/{sim.sysname:s}.dcd")
        simulation.context.setPositions(pdb.positions)
        if sb is not None:
            _report_start(simulation, sb, "before minimisation")
        print("Minimizing energy.")
        simulation.minimizeEnergy()
        if sb is not None:
            _report_start(simulation, sb, "after minimisation", check=True)

    simulation.reporters.append(
        app.dcdreporter.DCDReporter(f"{path}/{sim.sysname:s}.dcd", sim.wfreq, append=append))
    simulation.reporters.append(
        app.statedatareporter.StateDataReporter(
            f"{path}/{sim.sysname}.log", sim.logfreq, step=True, speed=True, elapsedTime=True,
            potentialEnergy=sim.report_potential_energy, separator="\t", append=append))

    # ---- the reaction loop -------------------------------------------------
    total = int(sim.steps)
    n_save = int(sim.wfreq)
    checkpoint_every = max(1, total // 10)   # same cadence as the stock path
    next_checkpoint = checkpoint_every
    # Each reactor keeps its own check cadence; the loop stops at the sooner.
    check_every = min(r.settings.check_every for r in (xl, sb) if r is not None)
    ramp_steps = min(r.settings.ramp_steps for r in (xl, sb) if r is not None)
    next_check_xl = simulation.currentStep + settings.check_every if xl is not None else None
    next_check_sb = simulation.currentStep + surface.check_every if sb is not None else None

    def save_all(step: int, done: int) -> None:
        simulation.saveCheckpoint(fcheck_out)
        # State first, outputs after: the checkpoint and the network it
        # belongs to have to land together.
        if xl is not None:
            xl.save_state(runtime_dir, step)
        if sb is not None:
            sb.save_state(runtime_dir, step)
        # `step` (the simulation's own counter, restart segments included), not
        # `done`: the events span the whole trajectory, so "steps simulated"
        # must too.
        write_event_record(runtime_dir, n_save, step, xl, sb, n_chains)

    print("STARTING SIMULATION", flush=True)
    done = 0
    with tqdm(total=total, mininterval=1) as bar:
        while done < total:
            ramping = bool(xl is not None and xl.ramping) or bool(sb is not None and sb.ramping)
            step = simulation.currentStep
            until = [next_checkpoint - done] + \
                [nc - step for nc in (next_check_xl, next_check_sb) if nc is not None]
            chunk = _chunk_sizes(total - done, check_every, ramp_steps, ramping, until)
            simulation.step(chunk)
            done += chunk
            step = simulation.currentStep

            # A check that falls due is always consumed (next_check advanced),
            # whether or not start_step has been reached; otherwise the chunk
            # clipping above would stop every single step until start_step.
            if sb is not None:
                changed = sb.advance_ramps(step)
                if step >= next_check_sb:
                    next_check_sb = step + surface.check_every
                    if step >= surface.start_step:
                        # Surface first: a lysine that binds the surface this
                        # check is not available to crosslink in the same check.
                        changed |= sb.check_reactions(
                            simulation.context, step,
                            blocked=(xl.used > 0) if xl is not None else None)
                if changed:
                    sb.force.updateParametersInContext(simulation.context)

            if xl is not None:
                changed = xl.advance_ramps(step)
                if step >= next_check_xl:
                    next_check_xl = step + settings.check_every
                    if step >= settings.start_step:
                        changed |= xl.check_reactions(
                            simulation.context, step,
                            blocked=sb.pinned if sb is not None else None)
                if changed:
                    xl.force.updateParametersInContext(simulation.context)

            if done >= next_checkpoint or done >= total:
                save_all(step, done)
                next_checkpoint = done + checkpoint_every
            bar.update(chunk)

    save_all(simulation.currentStep, done)

    stamp = datetime.now().strftime("%Y%d%m_%Hh%Mm%Ss")
    state_final = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    app.pdbreporter.PDBReporter(f"{path}/{sim.sysname}_{stamp}.pdb", 0).report(simulation, state_final)
    app.pdbreporter.PDBReporter(f"{path}/checkpoint.pdb", 0).report(simulation, state_final)

    if xl is not None:
        n = len(xl.events)
        intra = sum(1 for r in xl.event_rows(n_save) if r["kind"] == "intra")
        print(f"crosslinking: {n} crosslinks formed ({intra} intra, {n - intra} inter); "
              f"see {EVENTS_FILENAME} and {SUMMARY_FILENAME}", flush=True)
    if sb is not None:
        counts = sb.per_chain_counts(n_chains)
        print(f"surface: {sb.n_pinned} of {sb.n_sites} lysines bonded to the surface "
              f"({sb.n_pinned / sb.n_sites:.3f}); per chain min {int(counts.min())} / "
              f"max {int(counts.max())}"
              + (f"; cap reached at step {sb.cap_reached_step}" if sb.cap_reached_step is not None
                 else "; cap not reached")
              + f"; see {EVENTS_FILENAME} and {SUMMARY_FILENAME}", flush=True)
    return xl if xl is not None else sb


def _report_start(simulation, sb, label: str, check: bool = False) -> None:
    """Potential energy and pin placement at the start, so a bad construction is
    seen here rather than as a NaN a thousand steps in."""
    from openmm.unit import kilojoule_per_mole, nanometer

    from tools.surface import FORCE_GROUP

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    total = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    tether = simulation.context.getState(getEnergy=True, groups={FORCE_GROUP}) \
        .getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    pinned = np.flatnonzero(sb.pinned)
    off = np.linalg.norm(pos[sb.site_beads[pinned]] - sb.pin_xyz[pinned], axis=1) if len(pinned) else np.zeros(0)
    lowest = float(pos[:, 2].min())
    print(f"surface: potential energy {label}: {total:,.1f} kJ/mol total, "
          f"{tether:,.2f} kJ/mol in the {len(pinned)} surface tethers; pinned lysines at most "
          f"{off.max() if len(off) else 0.0:.4f} nm off their pins; lowest bead z = {lowest:.3f} nm "
          f"({lowest * 10:.1f} A)", flush=True)
    if not np.isfinite(total):
        raise RuntimeError(f"potential energy is not finite {label} — the initial configuration "
                           "is broken; not stepping.")
    if check and len(off) and off.max() > 0.5:
        raise RuntimeError(f"a pinned lysine ended {off.max():.3f} nm off its pin {label}; the "
                           "construction is being pulled apart — not stepping.")


def post_hoc_events(traj_path: Path, top_path: Path, settings: CrosslinkSettings,
                    n_save: int, seed: int = 0, verbose: bool = True,
                    drop_anchor: bool = True,
                    surface_rows: list[dict[str, Any]] | None = None) -> "Crosslinker":
    """Apply the same reaction rule to a *finished* trajectory, frame by frame.

    This is the comparison the reactive mode exists to make. It is the
    first-passage rule a post-hoc script applies: walk the saved frames, bond any
    pair that has come within the cutoff, and never let it come apart. The
    difference is that here the bond changes nothing — the chains carry on in the
    trajectory as if unbonded, free to drift apart and meet new partners.

    So a post-hoc count over a non-reactive run should come out *higher* than the
    reactive run's, because it counts encounters a real network would have
    prevented by tethering the chains at the first bond. The size of that gap is
    the point; a post-hoc conversion near 100% mostly says the cutoff is generous
    and the chains were never restrained.

    Resolution caveat: this only sees saved frames, while a reactive run checks
    every `check_every` steps. With `wfreq` much larger than `check_every` the
    post-hoc pass samples the encounter history more coarsely, which cuts the
    other way — it can miss brief encounters. Compare like with like where you
    can (`crosslink_check_every` = `wfreq`) and say which you used.

    `drop_anchor` is False for a free/preattached run, where residue 0 is an
    ordinary residue. `surface_rows` — the run's own surface-bond events, if it
    had any — keep a lysine that is bonded to the surface out of the pair pool
    from the step it bound, exactly as the live run did.
    """
    import mdtraj as md

    traj = md.load(str(traj_path), top=str(top_path))
    sites = find_sites(traj.topology, settings.selection, verbose=verbose, drop_anchor=drop_anchor)
    xl = Crosslinker(settings, sites, seed)     # force stays None: nothing to ramp

    if traj.unitcell_lengths is None:
        raise ValueError(f"{traj_path} has no box information, so the minimum image convention "
                         "cannot be applied.")

    # step at which each site became surface-bonded (inf = never)
    bound_at = np.full(len(sites), np.inf)
    if surface_rows:
        by_key = {(s.chain, s.resid): i for i, s in enumerate(sites)}
        for row in surface_rows:
            if row.get("kind") != "surface":
                continue
            slot = by_key.get((int(row["chain_i"]), int(row["resid_i"])))
            if slot is None:
                # A surface row naming a lysine this topology does not have is
                # never harmless: skipping it silently would let the replay
                # crosslink a lysine the live run had already bonded to the
                # surface — exactly the "bonded to the ground AND to another
                # lysine" state this masking exists to prevent. It means the
                # events file and the trajectory are not from the same run
                # (regenerated topology, hand-edited CSV), so stop.
                known = sorted(by_key)
                raise ValueError(
                    f"surface event for chain {row['chain_i']} resid {row['resid_i']} matches no "
                    f"lysine in {top_path} (it has {known[:8]}{'...' if len(known) > 8 else ''}). "
                    f"The events file and the trajectory are not from the same run; replaying "
                    f"them together would crosslink a surface-bonded lysine."
                )
            # time_ps is step x 0.01 rounded; divide back and round to the
            # integer step, or float error (e.g. 1000000.0000000001) could
            # leave the lysine unblocked at the very frame it bound.
            bound_at[slot] = min(bound_at[slot], int(round(float(row["time_ps"]) / DT_PS)))

    for frame in range(traj.n_frames):
        # DCDReporter writes frame k at step (k + 1) * n_save — the first frame
        # is n_save steps in, not step 0 — so this is the step the live run's
        # start_step and surface-bond times must be compared against.
        step = (frame + 1) * n_save
        if step < settings.start_step:
            continue
        xl.react(np.asarray(traj.xyz[frame], dtype=np.float64),
                 np.asarray(traj.unitcell_lengths[frame], dtype=np.float64), step,
                 blocked=(bound_at <= step) if surface_rows else None)
    if verbose:
        rows = xl.event_rows(n_save)
        intra = sum(1 for r in rows if r["kind"] == "intra")
        print(f"post-hoc: {len(rows)} events over {traj.n_frames} frames "
              f"({intra} intra, {len(rows) - intra} inter)")
    return xl
