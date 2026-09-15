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
lysines) is 2016 bonds, which is nothing. Above roughly 500 sites (125k bonds)
this needs replacing with a pair set restricted at setup by a generous distance
cutoff — you get a warning at that point.

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

from tools.results import publish

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

CSV_COLUMNS = [
    "frame", "time_ps", "resid_i", "chain_i", "resid_j", "chain_j",
    "kind", "span", "x", "y", "z", "pymol_i", "pymol_j",
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

    @property
    def distance_angstrom(self) -> float:
        return self.distance * 10.0

    def banner(self) -> str:
        """The startup line. Prints both units so nm/A cannot be mixed up."""
        return (f"crosslinking: ON, reaction distance {self.distance:.2f} nm "
                f"({self.distance_angstrom:.1f} A), valence {self.valence}, p={self.prob}")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["selection"] = list(self.selection) if self.selection is not None else None
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CrosslinkSettings":
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
               verbose: bool = True) -> list[SiteInfo]:
    """The reactive beads: lysines by default, or exactly `selection` if given.

    Auto-detection tries residue name LYS, then residue name K, then bead name K,
    in that order, and stops at the first rule that matches anything. If none
    match it raises with the names actually present rather than falling back to
    something arbitrary — silently crosslinking the wrong beads is worse than
    not starting.

    Residue 0 of every chain is dropped whatever the rule says. That bead is the
    surface anchor: `sim new` tags it "Z", CALVADOS gives it zero effective mass
    and OpenMM holds it fixed, so bonding to it would tether a chain to an
    immovable point rather than to another chain.
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

    anchors = [r for r in chosen if r.res_in_chain == 0]
    sites = [r for r in chosen if r.res_in_chain != 0]

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
        # saved reacted flags by bond index.
        slot_a, slot_b = np.triu_indices(n, k=1)
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

    def advance_ramps(self, step: int) -> bool:
        """Push every ramping bond to the stiffness its age entitles it to."""
        from openmm.unit import kilojoule_per_mole, nanometer

        if not self.ramping:
            return False
        r0 = self.settings.r0 * nanometer
        unit_k = kilojoule_per_mole / nanometer ** 2
        done = []
        for b in list(self.ramping):
            frac = self._ramp_fraction(b, step)
            self.force.setBondParameters(int(b), int(self.pair_bead_i[b]),
                                         int(self.pair_bead_j[b]), r0,
                                         frac * self.settings.k * unit_k)
            if frac >= 1.0:
                done.append(b)
        for b in done:
            del self.ramping[b]
        return True

    def check_reactions(self, context, step: int) -> bool:
        """React every pair within the cutoff, reading positions from a live Context."""
        from openmm.unit import nanometer

        # enforcePeriodicBox=False keeps whole chains intact; the minimum image
        # is applied per pair below, which is what the distance criterion needs.
        state = context.getState(getPositions=True, enforcePeriodicBox=False)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
        return self.react(pos, np.diag(np.asarray(box)), step)   # orthorhombic box

    def react(self, pos: np.ndarray, box_lengths: np.ndarray, step: int) -> bool:
        """React every pair that is within the cutoff, closest pairs first.

        Pure geometry, so the identical rule can be replayed over the frames of a
        finished trajectory (`post_hoc_events`) and compared against what the
        reactive run actually did.
        """
        settings = self.settings

        spare = self.used < settings.valence
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
                self.ramping[int(b)] = step
            i, j = self.pair_bead_i[b], self.pair_bead_j[b]
            midpoint = 0.5 * (pos[i] + pos[j])
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
        step = int(payload.get("step", 0))
        # A bond caught mid-ramp resumes from the stiffness it had reached, so a
        # restart neither re-spikes it nor leaves it soft forever.
        self.ramping = {}
        for b, formed_at in payload.get("ramping", {}).items():
            elapsed = max(0, step - int(formed_at))
            self.ramping[int(b)] = step - elapsed
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
            # saved frame.
            frame = int(round(event.step / n_save)) if n_save else 0
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

    def summary_lines(self, n_save: int, steps_done: int) -> list[str]:
        rows = self.event_rows(n_save)
        n = len(rows)
        intra = [r for r in rows if r["kind"] == "intra"]
        inter = [r for r in rows if r["kind"] == "inter"]
        spans = sorted(int(r["span"]) for r in intra)
        n_sites = len(self.sites)
        capacity = n_sites * self.settings.valence

        lines = [
            f"crosslink summary — written {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
            "",
            f"steps simulated:       {steps_done}",
            f"reactive sites:        {n_sites}",
            f"candidate pairs:       {self.n_pairs}",
            f"crosslinks formed:     {n}",
            # Conversion is per *site*: each bond consumes one valence slot at
            # each end, so 2 x events out of n_sites x valence available slots.
            f"conversion:            {2 * n / capacity:.4f}  (2 x events / (sites x valence))"
            if capacity else "conversion:            n/a",
            f"saturated sites:       {int(np.sum(self.used >= self.settings.valence))} of {n_sites}",
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
            f"selection:       {'auto (lysines)' if self.settings.selection is None else 'explicit'}",
            f"reaction seed:   {self.seed}",
        ]
        return lines

    def write_summary(self, runtime_dir: Path, n_save: int, steps_done: int) -> Path:
        path = Path(runtime_dir) / SUMMARY_FILENAME
        path.write_text("\n".join(self.summary_lines(n_save, steps_done)) + "\n")
        return path

    def write_outputs(self, runtime_dir: Path, n_save: int, steps_done: int) -> None:
        """Both output files. Called at every checkpoint, so a killed run keeps its log.

        The same two files are copied into the simulation's tracked ``results/``
        folder, so the conversion and the intra/inter split survive `sim clean`
        and reach GitHub with the run that produced them. The checkpoint
        (``crosslink_state.json``) is deliberately not copied: it is restart
        state, not a result.
        """
        self.write_events(runtime_dir, n_save)
        self.write_summary(runtime_dir, n_save, steps_done)
        publish(runtime_dir, SUMMARY_FILENAME, EVENTS_FILENAME)


def _chunk_sizes(remaining: int, settings: CrosslinkSettings, ramping: bool) -> int:
    """How far to step before stopping to update bond parameters.

    Normally one check interval. While a bond is mid-ramp, smaller steps, so the
    ramp is resolved over several parameter updates instead of jumping to full
    stiffness at the next check (see RAMP_UPDATES).
    """
    chunk = settings.check_every
    if ramping and settings.ramp_steps > 0:
        chunk = min(chunk, max(1, settings.ramp_steps // RAMP_UPDATES))
    return min(chunk, remaining)


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


def run_reactive(sim, settings: CrosslinkSettings, runtime_dir: Path | None = None) -> "Crosslinker":
    """Run `sim` with reactive crosslinking. Mirrors calvados Sim.simulate().

    Same integrator, platform, reporters, checkpointing and final PDBs as the
    stock path — the only difference is that stepping is broken into chunks so
    bond parameters can be updated between them. Chunking does not change the
    dynamics (stepping 500 + 500 is stepping 1000) and reporters fire on the
    simulation's own step counter, so the trajectory is unaffected by it.

    Only the configuration `sim new` actually generates is supported; the
    equilibration and clock-time modes are refused rather than half-replicated.
    """
    import os

    import openmm
    from openmm import app, unit
    from tqdm import tqdm

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

    print(settings.banner(), flush=True)
    sites = find_sites(sim.top, settings.selection)
    seed = _reaction_seed(settings, sim)
    xl = Crosslinker(settings, sites, seed)

    fcheck_in = f"{path}/{sim.frestart}"
    fcheck_out = f"{path}/restart.chk"
    restarting = os.path.isfile(fcheck_in) and sim.restart == "checkpoint"
    append = False

    if restarting:
        # Restoring a network is the whole reason crosslink_state.json exists: a
        # restart that rebuilt every bond dormant would silently make a
        # different network from the one the run had already formed.
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
    else:
        # A stale state file from a previous attempt must not resurrect bonds
        # into a run that is starting from scratch.
        stale = runtime_dir / STATE_FILENAME
        if stale.is_file():
            print(f"crosslinking: ignoring stale {STATE_FILENAME} (not restarting from a checkpoint)")
            stale.unlink()

    xl.add_force(sim.system)
    print(f"crosslinking: {xl.n_pairs} dormant bonds added, checks every "
          f"{settings.check_every} steps ({settings.check_every * DT_PS:.1f} ps)"
          + (f", ignored before step {settings.start_step}" if settings.start_step else ""),
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
        print("Minimizing energy.")
        simulation.minimizeEnergy()

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
    next_check = simulation.currentStep + settings.check_every

    print("STARTING SIMULATION", flush=True)
    done = 0
    with tqdm(total=total, mininterval=1) as bar:
        while done < total:
            chunk = _chunk_sizes(total - done, settings, bool(xl.ramping))
            chunk = max(1, min(chunk, next_checkpoint - done))
            simulation.step(chunk)
            done += chunk
            step = simulation.currentStep

            changed = xl.advance_ramps(step)
            if step >= next_check and step >= settings.start_step:
                changed |= xl.check_reactions(simulation.context, step)
                next_check = step + settings.check_every
            if changed:
                xl.force.updateParametersInContext(simulation.context)

            if done >= next_checkpoint or done >= total:
                simulation.saveCheckpoint(fcheck_out)
                # State first, outputs after: the checkpoint and the network it
                # belongs to have to land together.
                xl.save_state(runtime_dir, step)
                xl.write_outputs(runtime_dir, n_save, done)
                next_checkpoint = done + checkpoint_every
            bar.update(chunk)

    simulation.saveCheckpoint(fcheck_out)
    xl.save_state(runtime_dir, simulation.currentStep)
    xl.write_outputs(runtime_dir, n_save, done)

    stamp = datetime.now().strftime("%Y%d%m_%Hh%Mm%Ss")
    state_final = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    app.pdbreporter.PDBReporter(f"{path}/{sim.sysname}_{stamp}.pdb", 0).report(simulation, state_final)
    app.pdbreporter.PDBReporter(f"{path}/checkpoint.pdb", 0).report(simulation, state_final)

    n = len(xl.events)
    intra = sum(1 for r in xl.event_rows(n_save) if r["kind"] == "intra")
    print(f"crosslinking: {n} crosslinks formed ({intra} intra, {n - intra} inter); "
          f"see {EVENTS_FILENAME} and {SUMMARY_FILENAME}", flush=True)
    return xl


def post_hoc_events(traj_path: Path, top_path: Path, settings: CrosslinkSettings,
                    n_save: int, seed: int = 0, verbose: bool = True) -> "Crosslinker":
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
    """
    import mdtraj as md

    traj = md.load(str(traj_path), top=str(top_path))
    sites = find_sites(traj.topology, settings.selection, verbose=verbose)
    xl = Crosslinker(settings, sites, seed)     # force stays None: nothing to ramp

    if traj.unitcell_lengths is None:
        raise ValueError(f"{traj_path} has no box information, so the minimum image convention "
                         "cannot be applied.")

    for frame in range(traj.n_frames):
        step = frame * n_save
        if step < settings.start_step:
            continue
        xl.react(np.asarray(traj.xyz[frame], dtype=np.float64),
                 np.asarray(traj.unitcell_lengths[frame], dtype=np.float64), step)
    if verbose:
        rows = xl.event_rows(n_save)
        intra = sum(1 for r in rows if r["kind"] == "intra")
        print(f"post-hoc: {len(rows)} events over {traj.n_frames} frames "
              f"({intra} intra, {len(rows) - intra} inter)")
    return xl
