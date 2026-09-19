"""Attaching chains to the surface through their lysines — the ``free`` and
``preattached`` modes.

The default, ``brush`` mode, is untouched by this module: residue 0 of every
chain is tagged "Z", CALVADOS gives it zero mass and OpenMM holds it fixed on a
lattice point. In the two modes here residue 0 is an ordinary residue and the
*only* thing holding a chain to the surface is one or more surface-bonded
lysines:

``free``
    At initialisation every chain has exactly one lysine bonded to the surface,
    picked at random from that chain's lysines. During the run further lysines
    bond when they come within ``distance`` of the surface, with probability
    ``prob`` per check, until the global cap ``max_fraction`` is reached.

``preattached``
    At initialisation a random fraction ``preattached_fraction`` of *all* the
    lysines in the system is bonded, subject to every chain getting at least
    one. Whether more lysines may bind during the run is ``dynamic``; by
    default the set of surface bonds is fixed at t = 0, which is what makes the
    mode different from ``free``.

What a surface bond is
----------------------
A 3D pin: a harmonic tether ``0.5 k |r - r_pin|^2`` from the lysine bead to a
fixed point ``r_pin`` on the plane ``z = z_wall + tether``, where ``z_wall`` is
the onset of CALVADOS's repulsive wall and ``tether`` is how far above it a
bound lysine sits. It holds x, y and z, so a bound lysine does not slide along
the surface, and it never breaks. Like the lysine-lysine crosslinks in
``tools.crosslink``, every lysine is given a dormant (k = 0) tether up front —
OpenMM cannot add forces to a live Context — and binding means ramping that
lysine's stiffness up over ``ramp_steps``. Bonds present at initialisation are
at full stiffness from the first step, with the lysine already sitting exactly
on its pin, so they inject no energy.

A lysine is in exactly one of three states: free, surface-bonded, or bonded to
another lysine. A surface-bonded lysine never crosslinks and a crosslinked
lysine never binds the surface; the two reactors share one site list and mask
each other out.

The reaction rule during a run is the vertical gap: a free lysine binds when
``|z - (z_wall + tether)| < distance``, and is pinned where it is in x, y and at
the tether height in z. ``distance`` is *not* ``crosslink_distance``; the two
are separate knobs and both dominate their result. Sweep them.

Initial configurations
----------------------
The starting geometry is a construction, not an equilibrium state. Each chain
is laid out as a stack of vertical hairpins on a cubic grid of the CALVADOS bond
length, standing on the tether plane inside the chain's own lattice cell: every
pinned lysine sits exactly on its pin at the bottom of a column, no bead is
below the tether plane, and no two beads are closer than one bond length. That
is what keeps the first integration step from blowing up — the potential energy
is printed before and after minimisation, and the pins are checked to still be
on their targets, before a single step is taken. Chains are spread on the same
lattice brush mode uses, so the pins are spread with them; within one chain
the pins sit inside that chain's footprint, as they would for an adsorbed coil.
Equilibrate before sampling, and set ``start_step`` past that point: beads
left on the tether plane by the construction are within any reasonable reaction
distance from step 0.

Read before believing any of it
-------------------------------
* CALVADOS has no explicit crosslinker, no solvent and no activation barrier,
  so binding is a pure distance criterion on coarse-grained beads. The sequence
  and topology of binding events is meaningful; the *rate* is not.
* CALVADOS2 has no temperature-dependent hydrophobicity, so the LCST transition
  that drives real ELP adsorption is not reproduced.
* The wall is CALVADOS's purely repulsive ``step(z_wall-z)`` potential unless
  ``attraction`` is set, so how often a second lysine reaches the surface in
  ``free`` mode is set by diffusion and by ``distance``, not by any chemistry.
* ``distance`` and ``tether`` dominate every result; sweep them and report the
  sensitivity rather than one number.
* The initial configurations in both modes are constructions. Equilibrate
  before sampling.

Units are nanometres throughout, matching OpenMM and the box planning — not
the angstroms an MDAnalysis-based post-hoc script works in. The startup line
prints both.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from yaml import safe_dump, safe_load

from tools.crosslink import frame_of_step

# Written by prepare.py next to the trajectory; read by run.py. Its absence is
# what makes a run "brush", exactly as the absence of crosslink.yaml makes a run
# non-reactive.
SURFACE_FILENAME = "surface.yaml"
STATE_FILENAME = "surface_state.json"

MODES = ("brush", "free", "preattached")

# Same fixed Langevin timestep as tools.crosslink (calvados/sim.py).
DT_PS = 0.01

# CALVADOS's CA-CA bond length: the grid the initial construction is built on.
BOND_L = 0.38  # nm

# Geometry of the initial construction, in bond lengths. Columns stand two bond
# lengths apart so parallel legs are 0.76 nm apart — comfortably outside a
# bead's ~0.6 nm diameter — and intermediate hairpin bottoms turn 4 bond
# lengths (1.5 nm) above the tether plane so only the pins and their feet touch
# it. Columns are capped at BUILD_HEIGHT_NM so a long chain becomes a block of
# hairpins rather than one spike through the box.
COLUMN_PITCH = 2
FOLD_HEIGHT = 4
BUILD_HEIGHT_NM = 10.0

# The OpenMM force group the tethers are put in, so their energy can be read on
# its own and shown to be ~0 at the start.
FORCE_GROUP = 7


@dataclass(frozen=True)
class SurfaceSettings:
    """Everything that defines how chains attach through their lysines."""

    mode: str                                   # "free" | "preattached"
    distance: float = 0.8                       # nm — |z - tether plane| to bind
    prob: float = 1.0                           # P(bind | within distance at a check)
    max_fraction: float = 0.5                   # global cap, fraction of all lysines
    preattached_fraction: float | None = None   # preattached mode: fraction bound at t=0
    dynamic: bool | None = None                 # may lysines bind during the run? None = mode default
    k: float = 2000.0                           # kJ/mol/nm^2, tether stiffness
    tether: float = 0.6                         # nm above the wall onset a bound lysine sits
    ramp_steps: int = 500                       # steps over which k rises 0 -> k
    check_every: int = 1000                     # MD steps between binding checks
    start_step: int = 0                         # ignore binding before this step
    attraction: float = 0.0                     # kJ/mol depth of an optional lysine-wall well (0 = off)
    attraction_width: float = 0.5               # nm, width of that well
    seed: int | None = None                     # RNG for pin choice + binding; None = drawn at prepare
    z_wall: float | None = None                 # nm, wall onset; filled from the config at prepare

    def __post_init__(self) -> None:
        if self.mode not in ("free", "preattached"):
            raise ValueError(f"surface mode must be 'free' or 'preattached', got {self.mode!r} "
                             f"(brush mode has no surface settings).")
        if self.distance is None or self.distance <= 0:
            raise ValueError("surface_distance must be a positive number of nanometres.")
        if self.distance > 5.0:
            print(f"⚠  surface_distance is {self.distance} nm. That is very large; did you mean "
                  f"angstroms? ({self.distance} A = {self.distance / 10} nm)")
        if not 0.0 < self.prob <= 1.0:
            raise ValueError("surface_prob must be in (0, 1].")
        if not 0.0 < self.max_fraction <= 1.0:
            raise ValueError("surface_max_fraction must be in (0, 1].")
        if self.mode == "preattached":
            if self.preattached_fraction is None:
                raise ValueError("preattached mode needs surface_preattached_fraction.")
            if not 0.0 < self.preattached_fraction <= 1.0:
                raise ValueError("surface_preattached_fraction must be in (0, 1].")
            if self.preattached_fraction > self.max_fraction + 1e-12:
                raise ValueError(
                    f"surface_preattached_fraction ({self.preattached_fraction}) exceeds the global "
                    f"cap surface_max_fraction ({self.max_fraction}).")
        if self.k <= 0:
            raise ValueError("surface_k must be positive.")
        if self.tether <= 0:
            raise ValueError("surface_tether must be positive (nm above the wall onset).")
        if self.ramp_steps < 0:
            raise ValueError("surface_ramp_steps cannot be negative.")
        if self.check_every < 1:
            raise ValueError("surface_check_every must be at least 1 step.")
        if self.start_step < 0:
            raise ValueError("surface_start_step cannot be negative.")
        if self.attraction < 0:
            raise ValueError("surface_attraction is a well depth in kJ/mol and cannot be negative.")
        if self.attraction_width <= 0:
            raise ValueError("surface_attraction_width must be positive.")

    # ---- derived ----------------------------------------------------------

    @property
    def is_dynamic(self) -> bool:
        """Whether lysines may bind during the run. free: yes; preattached: no, unless asked."""
        if self.dynamic is not None:
            return bool(self.dynamic)
        return self.mode == "free"

    @property
    def distance_angstrom(self) -> float:
        return self.distance * 10.0

    @property
    def z_pin(self) -> float:
        """Height of the tether plane, nm: wall onset + tether."""
        if self.z_wall is None:
            raise ValueError("z_wall is not set; write_settings fills it from the config.")
        return self.z_wall + self.tether

    def banner(self) -> str:
        return (f"surface attachment: mode {self.mode}, reaction distance {self.distance:.2f} nm "
                f"({self.distance_angstrom:.1f} A), tether {self.tether:.2f} nm "
                f"({self.tether * 10:.1f} A) above the wall at {self.z_wall} nm, p={self.prob}, "
                f"cap {self.max_fraction:.3f} of all lysines, "
                f"binding during the run: {'on' if self.is_dynamic else 'off'}"
                + (f", lysine-wall well {self.attraction} kJ/mol" if self.attraction else ""))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SurfaceSettings":
        known = {f: data[f] for f in cls.__dataclass_fields__ if f in data}
        return cls(**known)


def write_settings(runtime_dir: Path, settings: SurfaceSettings | None, *,
                   z_wall: float | None = None, run_seed: int | None = None) -> Path | None:
    """Record the surface settings next to the trajectory, or remove a stale file.

    Nothing is written for brush mode: the *absence* of surface.yaml is what
    run.py reads as brush, so a regenerated brush run can never inherit a
    surface.yaml from an earlier attempt.

    The seed has to be fixed here, at prepare time, not at run time: it decides
    which lysines are pinned and where, and both a fresh start and a restart
    must agree on that. None becomes the run's OpenMM seed, else a recorded
    random draw.
    """
    path = Path(runtime_dir) / SURFACE_FILENAME
    if settings is None:
        path.unlink(missing_ok=True)
        return None
    data = settings.to_dict()
    if z_wall is not None:
        data["z_wall"] = float(z_wall)
    if data.get("z_wall") is None:
        raise ValueError("surface settings need z_wall (the wall onset from the config).")
    if data.get("seed") is None:
        data["seed"] = int(run_seed) if run_seed is not None \
            else int(np.random.SeedSequence().entropy % (2 ** 31))
    SurfaceSettings.from_dict(data)   # validate the final form
    with open(path, "w") as f:
        safe_dump(data, f, sort_keys=True)
    return path


def load_settings(runtime_dir: Path) -> SurfaceSettings | None:
    """The run's surface settings, or None when the run is brush mode."""
    path = Path(runtime_dir) / SURFACE_FILENAME
    if not path.is_file():
        return None
    with open(path) as f:
        data = safe_load(f) or {}
    if not data or data.get("mode") in (None, "brush"):
        return None
    return SurfaceSettings.from_dict(data)


def mode_of(runtime_dir: Path) -> str:
    """'brush', 'free' or 'preattached' for a prepared run."""
    settings = load_settings(runtime_dir)
    return settings.mode if settings else "brush"


# ---------------------------------------------------------------------------
# Choosing the initial pins
# ---------------------------------------------------------------------------

def minimum_fraction(n_chains: int, n_sites: int) -> float:
    """The smallest cap that still lets every chain have one bond."""
    return n_chains / n_sites if n_sites else math.inf


def check_capacity(settings: SurfaceSettings, n_chains: int, n_sites: int) -> int:
    """The cap in lysines, after checking it is compatible with one bond per chain.

    Raises with the minimum spelled out rather than leaving chains unattached.
    """
    if n_sites == 0:
        raise ValueError("No lysines in the system, so nothing can attach to the surface. "
                         "free/preattached modes need lysines in the sequence.")
    chains_without = n_chains - min(n_chains, n_sites)
    n_max = int(math.floor(settings.max_fraction * n_sites + 1e-9))
    needed = minimum_fraction(n_chains, n_sites)
    if n_max < n_chains:
        raise ValueError(
            f"surface_max_fraction = {settings.max_fraction} allows {n_max} of {n_sites} lysines to "
            f"bind, but every one of the {n_chains} chains needs at least one — set it to at least "
            f"{needed:.4f} ({n_chains}/{n_sites}).")
    if settings.mode == "preattached":
        n_pre = int(round(settings.preattached_fraction * n_sites))
        if n_pre < n_chains:
            raise ValueError(
                f"surface_preattached_fraction = {settings.preattached_fraction} would pin "
                f"{n_pre} of {n_sites} lysines, fewer than the {n_chains} chains that each need one "
                f"— set it to at least {needed:.4f} ({n_chains}/{n_sites}).")
    if chains_without:
        raise ValueError(f"{chains_without} chain(s) have no lysine at all and cannot attach.")
    return n_max


def choose_initial_pins(settings: SurfaceSettings, chain_of_site: Sequence[int],
                        n_chains: int, rng: np.random.Generator) -> np.ndarray:
    """Which sites (indices into the site list) are bonded at t = 0.

    Every chain gets one lysine chosen uniformly from its own; preattached mode
    then tops up to ``round(preattached_fraction x n_sites)`` uniformly from the
    lysines still free. Returned sorted.
    """
    chain_of_site = np.asarray(chain_of_site, dtype=int)
    n_sites = len(chain_of_site)
    n_max = check_capacity(settings, n_chains, n_sites)

    chosen: list[int] = []
    for chain in range(n_chains):
        own = np.flatnonzero(chain_of_site == chain)
        chosen.append(int(rng.choice(own)))

    if settings.mode == "preattached":
        n_pre = min(n_max, int(round(settings.preattached_fraction * n_sites)))
        extra = n_pre - len(chosen)
        if extra > 0:
            free = np.setdiff1d(np.arange(n_sites), np.array(chosen))
            chosen.extend(int(i) for i in rng.choice(free, size=extra, replace=False))
    return np.array(sorted(chosen), dtype=int)


# ---------------------------------------------------------------------------
# Building the initial conformations
# ---------------------------------------------------------------------------

def _serpentine(n_slots: int, width: int) -> list[tuple[int, int]]:
    """Grid (u, v) of column slot k, boustrophedon over rows of `width` slots."""
    out = []
    for k in range(n_slots):
        row, col = divmod(k, width)
        if row % 2:
            col = width - 1 - col
        out.append((col * COLUMN_PITCH, row * COLUMN_PITCH))
    return out


class _Path:
    """A chain path on the (t, h) grid: t = horizontal path coordinate along the
    serpentine of column slots (COLUMN_PITCH units per slot), h = height in bonds."""

    def __init__(self) -> None:
        self.t: list[int] = []
        self.h: list[int] = []

    def start(self, t: int, h: int) -> None:
        self.t.append(t)
        self.h.append(h)

    def move(self, dt: int, dh: int) -> None:
        """Walk |dt| horizontal then |dh| vertical unit steps (one bead per step)."""
        for _ in range(abs(dt)):
            self.t.append(self.t[-1] + (1 if dt > 0 else -1))
            self.h.append(self.h[-1])
        for _ in range(abs(dh)):
            self.t.append(self.t[-1])
            self.h.append(self.h[-1] + (1 if dh > 0 else -1))

    def __len__(self) -> int:
        return len(self.t)


def _hairpin_plan(n: int, o_a: int, o_b: int, hmax: int) -> tuple[list[int], int] | None:
    """Column tops for a segment of `n` bonds between pins at slot offsets o_a, o_b.

    Returns (tops, n_hairpins) or None when the segment is too short for any
    hairpin and must lie flat. The bond count of a segment with H hairpins is
        (p - o_a) + sum_i [(top_i - bot_i) + p + (top_i - bot'_i)] + (H - 1) p + o_b
    with bot_1 = bot'_H = 0 and every other bottom at FOLD_HEIGHT.
    """
    p, lo = COLUMN_PITCH, FOLD_HEIGHT
    for n_hp in range(1, n // 2 + 2):
        fixed = (p - o_a) + (2 * n_hp - 1) * p + o_b - (2 * n_hp - 2) * lo
        twice_sum = n - fixed
        if twice_sum < 0 or twice_sum % 2:
            continue
        total = twice_sum // 2
        min_top = 1 if n_hp == 1 else lo + 1
        if total < n_hp * min_top:
            return None          # too short to rise at all: lay it flat
        if total > n_hp * hmax:
            continue             # too tall for one more column: add a hairpin
        base, rem = divmod(total, n_hp)
        tops = [base + (1 if i < rem else 0) for i in range(n_hp)]
        return tops, n_hp
    return None


def _build_chain_grid(n_beads: int, pins: Sequence[int], hmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Lay one chain out as hairpin columns. Returns (t, h) integer arrays per bead.

    `pins` are bead indices (within the chain) that must sit at h = 0.
    """
    p = COLUMN_PITCH
    pins = sorted(set(int(i) for i in pins))
    if not pins:
        raise ValueError("a chain needs at least one pinned lysine")
    if pins[0] < 0 or pins[-1] >= n_beads:
        # Otherwise this surfaces as an AssertionError on the path length,
        # which says nothing about the site that caused it.
        raise ValueError(f"pinned bead index out of range for a {n_beads}-bead chain: {pins}")

    # Slot offset of each pin: 0 = on a column slot, 1 = one unit past it. Set
    # by the parity of the segment leading to it (a Manhattan path between two
    # points at h = 0 has the parity of their horizontal separation).
    path = _Path()

    # --- dangling start, built outward from the first pin and reversed -------
    first = pins[0]
    head = _Path()
    head.start(0, 0)
    remaining = first
    if remaining:
        # straight up the column above the pin, then hairpins outward (t decreasing)
        up = min(remaining, hmax)
        head.move(0, up)
        remaining -= up
        direction = -1
        level_top, level_bot = hmax, FOLD_HEIGHT
        going_down = True
        while remaining > 0:
            step = min(remaining, p)
            head.move(direction * step, 0)
            remaining -= step
            if remaining <= 0:
                break
            if going_down:
                dh = max(0, min(remaining, head.h[-1] - level_bot))
                head.move(0, -dh)
            else:
                dh = max(0, min(remaining, level_top - head.h[-1]))
                head.move(0, dh)
            remaining -= dh
            going_down = not going_down
    t0 = -min(head.t)                       # shift so the leftmost slot is t = 0
    ts = [t + t0 for t in reversed(head.t)]
    hs = list(reversed(head.h))
    path.t.extend(ts)
    path.h.extend(hs)
    t_pin = ts[-1]
    o_pin = 0                                # first pin sits on a slot

    # --- segments between pins ------------------------------------------------
    for a, b in zip(pins, pins[1:]):
        n = b - a
        o_b = (n + o_pin) % 2
        plan = _hairpin_plan(n, o_pin, o_b, hmax)
        if plan is None:
            # too short to rise: lie flat along the path
            path.move(n, 0)
            o_pin = o_b
        else:
            tops, n_hp = plan
            path.move(p - o_pin, 0)                    # foot to the next slot
            for i, top in enumerate(tops):
                bot = 0 if i == 0 else FOLD_HEIGHT
                bot_next = 0 if i == n_hp - 1 else FOLD_HEIGHT
                path.move(0, top - bot)                # up
                path.move(p, 0)                        # bridge at the top
                path.move(0, -(top - bot_next))        # down
                if i < n_hp - 1:
                    path.move(p, 0)                    # bridge at the fold height
            path.move(o_b, 0)                          # landing offset
            o_pin = o_b
        t_pin = path.t[-1]
        assert path.h[-1] == 0 and len(path) == b + 1, (len(path), b)

    # --- dangling end ----------------------------------------------------------
    remaining = n_beads - 1 - pins[-1]
    if remaining:
        step = min(remaining, p - o_pin)
        path.move(step, 0)
        remaining -= step
        going_up = True
        while remaining > 0:
            if going_up:
                dh = max(0, min(remaining, hmax - path.h[-1]))
                path.move(0, dh)
            else:
                dh = max(0, min(remaining, path.h[-1] - FOLD_HEIGHT))
                path.move(0, -dh)
            remaining -= dh
            if remaining <= 0:
                break
            step = min(remaining, p)
            path.move(step, 0)
            remaining -= step
            going_up = not going_up

    assert len(path) == n_beads, (len(path), n_beads)
    t = np.array(path.t, dtype=int)
    h = np.array(path.h, dtype=int)
    t -= t.min()
    return t, h


def _grid_to_xy(t: np.ndarray, width: int) -> np.ndarray:
    """Path coordinate -> (u, v) grid position on the serpentine of slots."""
    n_slots = int(t.max()) // COLUMN_PITCH + 2
    slots = _serpentine(n_slots, width)
    out = np.empty((len(t), 2), dtype=float)
    for i, ti in enumerate(t):
        k, r = divmod(int(ti), COLUMN_PITCH)
        a = np.array(slots[k], dtype=float)
        b = np.array(slots[k + 1], dtype=float)
        out[i] = a + (b - a) * (r / COLUMN_PITCH)
    return out


def build_chain(n_beads: int, pins: Sequence[int], hmax: int) -> np.ndarray:
    """Bead positions (nm) of one chain, in a local frame: pins at z = 0, centred in x/y.

    Hairpin columns on a cubic grid of BOND_L: every bond is exactly BOND_L
    long, every non-bonded pair is at least BOND_L apart, every pinned bead is
    at z = 0 and nothing is below it.
    """
    t, h = _build_chain_grid(n_beads, pins, hmax)
    n_slots = int(t.max()) // COLUMN_PITCH + 2
    width = max(1, int(math.ceil(math.sqrt(n_slots))))
    uv = _grid_to_xy(t, width)
    xyz = np.column_stack([uv[:, 0], uv[:, 1], h.astype(float)]) * BOND_L
    centre = 0.5 * (xyz[:, :2].min(axis=0) + xyz[:, :2].max(axis=0))
    xyz[:, :2] -= centre
    return xyz


def check_geometry(pos: np.ndarray, pin_beads: Sequence[int], pin_xyz: np.ndarray,
                   z_floor: float, box: Sequence[float], min_sep: float = BOND_L * 0.999) -> dict[str, float]:
    """Verify a constructed configuration before it goes anywhere near an integrator.

    Raises on: a bead below `z_floor`, a pinned bead off its pin, or two beads
    closer than `min_sep` (a bead-on-bead overlap is the NaN in waiting).
    Returns the numbers it checked for the startup report.
    """
    pos = np.asarray(pos, dtype=float)
    lowest = float(pos[:, 2].min())
    if lowest < z_floor - 1e-9:
        raise ValueError(f"initial configuration has a bead at z = {lowest:.3f} nm, below the "
                         f"floor at {z_floor:.3f} nm.")
    if len(pin_beads):
        off = np.linalg.norm(pos[np.asarray(pin_beads)] - np.asarray(pin_xyz), axis=1)
        if off.max() > 1e-6:
            raise ValueError(f"a pinned lysine is {off.max():.4f} nm off its pin in the initial "
                             "configuration.")
    highest = float(pos[:, 2].max())
    if highest > box[2] - 1.0:
        raise ValueError(f"initial configuration reaches z = {highest:.2f} nm, within 1 nm of the "
                         f"{box[2]} nm box top.")
    # Closest non-identical pair. cKDTree handles tens of thousands of beads in
    # well under a second; the box is periodic in x/y but the chains are laid
    # out inside their own cells, so the minimum image only matters across the
    # boundary — handled by boxsize on the lateral axes.
    from scipy.spatial import cKDTree

    # The z axis is periodic too, but no *false* overlap can hide across it:
    # the checks above guarantee every bead is >= z_floor (2.5 nm) above the
    # bottom and >= 1 nm below the top, so the closest pair across the z
    # boundary is > 3.5 nm apart, far outside min_sep. Passing 4 x the height
    # as the z period makes the tree treat z as effectively non-periodic while
    # still accepting the coordinates (cKDTree needs 0 <= x < boxsize on every
    # periodic axis; z is in [z_floor, box_z - 1] so it qualifies).
    lateral = np.array([box[0], box[1]], dtype=float)
    wrapped = pos.copy()
    wrapped[:, :2] = np.mod(wrapped[:, :2], lateral)
    # np.mod(-tiny, L) rounds to exactly L, which cKDTree rejects as outside
    # the box; fold that edge back to 0.
    wrapped[:, :2] = np.where(wrapped[:, :2] >= lateral, 0.0, wrapped[:, :2])
    tree = cKDTree(wrapped, boxsize=[box[0], box[1], box[2] * 4.0])
    pairs = tree.query_pairs(min_sep, output_type="ndarray")
    n_close = int(len(pairs))
    d_min = float(tree.query(wrapped, k=2)[0][:, 1].min())
    if n_close:
        i, j = pairs[0]
        raise ValueError(f"{n_close} bead pair(s) closer than {min_sep:.3f} nm in the initial "
                         f"configuration (e.g. beads {i} and {j}); chains overlap — lower the "
                         "concentration or the chain count.")
    return {"lowest_z": lowest, "highest_z": highest, "min_pair_distance": d_min}


@dataclass
class InitialPlan:
    """What build_sim decided: which lysines are pinned at t = 0 and where."""

    site_beads: list[int]         # every lysine bead index (the site list order)
    chain_of_site: list[int]
    pinned_sites: list[int]       # indices into site_beads
    pin_xyz: list[list[float]]    # pin coordinates, one per pinned site, nm
    seed: int
    footprint_nm: float           # lateral extent of one chain's construction
    height_nm: float              # tallest column

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def place_chains(sim, settings: SurfaceSettings, *, nx: int, ny: int, spacing: float,
                 margin: float, verbose: bool = True) -> InitialPlan:
    """Build every chain's starting conformation for free/preattached mode.

    Sets ``sim.pos`` (the positions CALVADOS writes to top.pdb), stores the
    plan on ``sim.surface_plan`` for the run loop, and returns it. Called from
    the generated prepare.py's build_sim in place of the brush layout.
    """
    from tools.crosslink import find_sites

    sites = find_sites(sim.top, None, verbose=False, drop_anchor=False)
    chain_of_site = np.array([s.chain for s in sites], dtype=int)
    n_chains = int(sim.top.n_chains)
    rng = np.random.default_rng(int(settings.seed))
    pinned = choose_initial_pins(settings, chain_of_site, n_chains, rng)

    z_pin = settings.z_pin
    box = [float(v) for v in sim.box]
    hmax = max(1, int(math.floor(min(BUILD_HEIGHT_NM, 0.5 * box[2] - z_pin) / BOND_L)))

    pos = np.asarray(sim.pos, dtype=float).copy()
    pin_beads: list[int] = []
    pin_xyz: list[list[float]] = []
    footprint = 0.0
    height = 0.0

    # Same lattice, same cell centres, as brush mode's build_sim.
    chain = 0
    ibead = 0
    for comp in sim.components:
        for _ in range(comp.nmol):
            n = int(comp.nbeads)
            own = [int(sites[s].res_in_chain) for s in pinned if chain_of_site[s] == chain]
            xyz = build_chain(n, own, hmax)
            x0 = margin + (chain % nx + 0.5) * spacing
            y0 = margin + (chain // nx + 0.5) * spacing
            xyz += np.array([x0, y0, z_pin])
            pos[ibead:ibead + n] = xyz
            footprint = max(footprint, float(np.ptp(xyz[:, 0])), float(np.ptp(xyz[:, 1])))
            height = max(height, float(xyz[:, 2].max() - z_pin))
            for s in pinned:
                if chain_of_site[s] == chain:
                    bead = sites[s].index
                    pin_beads.append(int(bead))
                    pin_xyz.append([float(v) for v in pos[bead]])
            ibead += n
            chain += 1

    numbers = check_geometry(pos, pin_beads, np.array(pin_xyz), z_pin, box)
    sim.pos = pos

    plan = InitialPlan(
        site_beads=[int(s.index) for s in sites],
        chain_of_site=[int(c) for c in chain_of_site],
        pinned_sites=[int(s) for s in pinned],
        pin_xyz=pin_xyz,
        seed=int(settings.seed),
        footprint_nm=footprint,
        height_nm=height,
    )
    sim.surface_plan = plan

    if verbose:
        per_chain = np.bincount(chain_of_site[pinned], minlength=n_chains)
        print(f"surface: mode {settings.mode}; {len(sites)} lysines in {n_chains} chains "
              f"(residue 0 is an ordinary residue, not an anchor)")
        print(f"surface: {len(pinned)} of {len(sites)} lysines pinned at t=0 "
              f"({len(pinned) / len(sites):.3f}); per chain min {per_chain.min()} / "
              f"max {per_chain.max()}; pins on z = {z_pin:.3f} nm ({z_pin * 10:.1f} A)")
        print(f"surface: initial construction {footprint:.2f} nm across, {height:.2f} nm tall, "
              f"lowest bead z = {numbers['lowest_z']:.3f} nm, closest pair "
              f"{numbers['min_pair_distance']:.3f} nm")
        if footprint > spacing:
            print(f"⚠  a chain's starting construction ({footprint:.2f} nm) is wider than the "
                  f"{spacing} nm lattice spacing — neighbours start interleaved. It minimises, "
                  f"but expect a long equilibration; a lower concentration is cleaner.")
    return plan


# ---------------------------------------------------------------------------
# The binder: the force, the state, the reaction
# ---------------------------------------------------------------------------

@dataclass
class SurfaceEvent:
    step: int
    site: int
    pin: tuple[float, float, float]
    distance: float          # gap to the tether plane when it bound (0 for initial)
    initial: bool

    def to_dict(self) -> dict[str, Any]:
        return {"step": self.step, "site": self.site, "pin": list(self.pin),
                "distance": self.distance, "initial": self.initial}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "SurfaceEvent":
        return cls(int(d["step"]), int(d["site"]), tuple(float(v) for v in d["pin"]),
                   float(d["distance"]), bool(d["initial"]))


class SurfaceBinder:
    """Every lysine's tether to the surface: dormant until it binds, then permanent."""

    def __init__(self, settings: SurfaceSettings, sites, seed: int):
        self.settings = settings
        self.sites = sites
        self.seed = int(seed)
        self.rng = np.random.default_rng(self.seed)
        n = len(sites)
        self.site_beads = np.array([s.index for s in sites], dtype=int)
        self.chain_of_site = np.array([s.chain for s in sites], dtype=int)
        self.pinned = np.zeros(n, dtype=bool)
        self.pin_xyz = np.zeros((n, 3), dtype=float)
        self.ramping: dict[int, int] = {}       # site -> step it bound
        self.events: list[SurfaceEvent] = []
        self.force = None
        self.resume_step = 0
        self.cap_reached_step: int | None = None
        self.n_max = 0

    # ---- setup ------------------------------------------------------------

    @property
    def n_sites(self) -> int:
        return len(self.site_beads)

    @property
    def n_pinned(self) -> int:
        return int(self.pinned.sum())

    def apply_plan(self, plan: InitialPlan, n_chains: int) -> None:
        """Adopt the pins build_sim constructed (a fresh start, not a restart)."""
        if [int(i) for i in plan.site_beads] != [int(i) for i in self.site_beads]:
            raise ValueError("the surface plan was built for a different set of lysines than "
                             "this system has.")
        self.n_max = check_capacity(self.settings, n_chains, self.n_sites)
        for s, xyz in zip(plan.pinned_sites, plan.pin_xyz):
            self.pinned[s] = True
            self.pin_xyz[s] = xyz
            self.events.append(SurfaceEvent(0, int(s), tuple(float(v) for v in xyz), 0.0, True))
        if self.n_pinned >= self.n_max:
            self.cap_reached_step = 0

    def add_force(self, system) -> Any:
        """One dormant tether per lysine, at full stiffness for the ones already bound.

        Called before the Context exists: pins present at t = 0 (or restored
        from a checkpoint) are at their stiffness from the first step, and their
        lysines are already on their pins, so nothing switches on abruptly.
        """
        from openmm import CustomExternalForce

        expr = "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"
        if self.settings.attraction > 0:
            # An optional lysine-wall well, centred on the tether plane, so a
            # lysine that wanders close is held near the surface long enough
            # to react. Off by default: the stock wall is purely repulsive.
            expr += (f" - {self.settings.attraction}*exp(-0.5*((z-{self.settings.z_pin})"
                     f"/{self.settings.attraction_width})^2)")
        force = CustomExternalForce(expr)
        force.setName("surface tethers")
        force.setForceGroup(FORCE_GROUP)
        for name in ("k", "x0", "y0", "z0"):
            force.addPerParticleParameter(name)
        for s in range(self.n_sites):
            k = self.settings.k * self._ramp_fraction(s, self.resume_step) if self.pinned[s] else 0.0
            force.addParticle(int(self.site_beads[s]), [k, *self.pin_xyz[s]])
        system.addForce(force)
        self.force = force
        return force

    # ---- the reaction -----------------------------------------------------

    def _ramp_fraction(self, site: int, step: int) -> float:
        formed_at = self.ramping.get(int(site))
        if formed_at is None or self.settings.ramp_steps <= 0:
            return 1.0
        return min(1.0, max(0.0, (step - formed_at) / self.settings.ramp_steps))

    def _set_particle(self, site: int, k: float) -> None:
        self.force.setParticleParameters(int(site), int(self.site_beads[site]),
                                         [k, *self.pin_xyz[site]])

    def advance_ramps(self, step: int) -> bool:
        if not self.ramping:
            return False
        done = []
        for s in list(self.ramping):
            frac = self._ramp_fraction(s, step)
            self._set_particle(s, frac * self.settings.k)
            if frac >= 1.0:
                done.append(s)
        for s in done:
            del self.ramping[s]
        return True

    def check_reactions(self, context, step: int, blocked: np.ndarray | None = None) -> bool:
        """Bind every free lysine within `distance` of the tether plane, nearest first."""
        from openmm.unit import nanometer

        state = context.getState(getPositions=True, enforcePeriodicBox=False)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        return self.react(pos, step, blocked)

    def react(self, pos: np.ndarray, step: int, blocked: np.ndarray | None = None) -> bool:
        settings = self.settings
        if not settings.is_dynamic or self.n_pinned >= self.n_max:
            return False
        free = ~self.pinned
        if blocked is not None:
            free &= ~np.asarray(blocked, dtype=bool)
        candidates = np.flatnonzero(free)
        if candidates.size == 0:
            return False
        z = pos[self.site_beads[candidates], 2]
        gap = np.abs(z - settings.z_pin)
        within = gap < settings.distance
        if not within.any():
            return False
        order = np.argsort(gap[within], kind="stable")
        sites = candidates[within][order]
        gaps = gap[within][order]

        changed = False
        for s, g in zip(sites, gaps):
            if self.n_pinned >= self.n_max:
                break
            if settings.prob < 1.0 and self.rng.random() > settings.prob:
                continue
            bead = self.site_beads[s]
            pin = (float(pos[bead, 0]), float(pos[bead, 1]), float(settings.z_pin))
            self.pinned[s] = True
            self.pin_xyz[s] = pin
            if self.force is not None:
                self.ramping[int(s)] = step
                self._set_particle(int(s), self.settings.k * self._ramp_fraction(int(s), step))
            self.events.append(SurfaceEvent(step, int(s), pin, float(g), False))
            changed = True
        if self.n_pinned >= self.n_max and self.cap_reached_step is None:
            self.cap_reached_step = int(step)
        return changed

    # ---- restart safety ---------------------------------------------------

    def save_state(self, runtime_dir: Path, step: int) -> Path:
        path = Path(runtime_dir) / STATE_FILENAME
        payload = {
            "step": int(step),
            "seed": self.seed,
            "settings": self.settings.to_dict(),
            "site_beads": [int(i) for i in self.site_beads],
            "pinned": np.flatnonzero(self.pinned).tolist(),
            "pin_xyz": {str(int(s)): [float(v) for v in self.pin_xyz[s]]
                        for s in np.flatnonzero(self.pinned)},
            "ramping": {str(int(s)): int(t) for s, t in self.ramping.items()},
            "events": [e.to_dict() for e in self.events],
            "cap_reached_step": self.cap_reached_step,
            "n_max": int(self.n_max),
        }
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=1))
        tmp.replace(path)
        return path

    def load_state(self, runtime_dir: Path) -> int | None:
        path = Path(runtime_dir) / STATE_FILENAME
        if not path.is_file():
            return None
        payload = json.loads(path.read_text())
        saved = [int(i) for i in payload.get("site_beads", [])]
        if saved != [int(i) for i in self.site_beads]:
            raise ValueError(
                f"{STATE_FILENAME} was written for a different set of lysines ({len(saved)}) than "
                f"this system has ({len(self.site_beads)}). Delete the runtime folder and start "
                "over rather than restarting into pins that mean something else.")
        self.pinned[:] = False
        self.pin_xyz[:] = 0.0
        for s in payload["pinned"]:
            self.pinned[int(s)] = True
            self.pin_xyz[int(s)] = payload["pin_xyz"][str(int(s))]
        self.events = [SurfaceEvent.from_dict(e) for e in payload.get("events", [])]
        if len(self.events) != self.n_pinned:
            raise ValueError(f"{STATE_FILENAME} is inconsistent: {len(self.events)} surface events "
                             f"but {self.n_pinned} pinned lysines. Do not restart from it.")
        step = int(payload.get("step", 0))
        self.ramping = {int(s): min(int(t), step) for s, t in payload.get("ramping", {}).items()}
        self.cap_reached_step = payload.get("cap_reached_step")
        if "n_max" in payload:
            self.n_max = int(payload["n_max"])
        else:
            # A state file from before n_max was saved. 0 would mean "cap
            # reached" and silently stop all binding, so recompute it; every
            # chain has a lysine (check_capacity enforced that at t = 0), so the
            # chain count is the highest chain index + 1.
            self.n_max = check_capacity(self.settings, int(self.chain_of_site.max()) + 1,
                                        self.n_sites)
        self.rng = np.random.default_rng([self.seed, step])
        self.resume_step = step
        return step

    # ---- outputs ----------------------------------------------------------

    def event_rows(self, n_save: int) -> list[dict[str, Any]]:
        """Surface events in the crosslink_events.csv schema (kind = 'surface')."""
        rows = []
        for e in self.events:
            a = self.sites[e.site]
            # 0-based nearest saved frame; an initial bond (step 0) predates the
            # first frame and goes on frame 0 — see CSV_COLUMNS in tools.crosslink.
            frame = frame_of_step(e.step, n_save)
            x, y, z = e.pin
            rows.append({
                "frame": frame,
                "time_ps": round(e.step * DT_PS, 6),
                "resid_i": a.resid,
                "chain_i": a.chain,
                "resid_j": "",
                "chain_j": "",
                "kind": "surface",
                "span": "",
                "x": round(x, 4), "y": round(y, 4), "z": round(z, 4),
                "pymol_i": a.index + 1,
                "pymol_j": "",
                "origin": "initial" if e.initial else "run",
            })
        return rows

    def per_chain_counts(self, n_chains: int) -> np.ndarray:
        return np.bincount(self.chain_of_site[self.pinned], minlength=n_chains)

    def summary_lines(self, n_chains: int) -> list[str]:
        counts = self.per_chain_counts(n_chains)
        hist = np.bincount(counts, minlength=4)
        n_run = sum(1 for e in self.events if not e.initial)
        n_init = len(self.events) - n_run
        lines = [
            "surface attachment",
            "------------------",
            f"mode:                  {self.settings.mode}",
            f"lysines total:         {self.n_sites}",
            f"surface-bonded:        {self.n_pinned}  ({self.n_pinned / self.n_sites:.4f} of all lysines)"
            if self.n_sites else "surface-bonded:        0",
            f"  at initialisation:   {n_init}",
            f"  during the run:      {n_run}",
            f"cap:                   {self.n_max} lysines ({self.settings.max_fraction} of all)"
            + (f", reached at step {self.cap_reached_step}" if self.cap_reached_step is not None
               else ", not reached"),
            f"chains with 0 bonds:   {int(hist[0])}",
            f"chains with 1 bond:    {int(hist[1])}",
            f"chains with 2 bonds:   {int(hist[2])}",
            f"chains with 3+ bonds:  {int(hist[3:].sum())}",
            f"bonds per chain:       min {int(counts.min())}, mean {counts.mean():.2f}, "
            f"max {int(counts.max())}" if n_chains else "",
            "",
            "The bonds-per-chain distribution is the number that says whether this is a",
            "coating (many pins per chain) or loosely tethered chains (one pin each).",
            "Event *times* are not kinetics — see the caveats in tools/surface.py.",
            "",
            f"distance:        {self.settings.distance} nm ({self.settings.distance_angstrom} A) "
            "to the tether plane",
            f"tether:          {self.settings.tether} nm above the wall at {self.settings.z_wall} nm "
            f"(pins on z = {self.settings.z_pin:.3f} nm)",
            f"prob:            {self.settings.prob}",
            f"binding in run:  {'on' if self.settings.is_dynamic else 'off'}",
            f"check_every:     {self.settings.check_every} steps",
            f"start_step:      {self.settings.start_step}",
            f"k:               {self.settings.k} kJ/mol/nm^2",
            f"ramp_steps:      {self.settings.ramp_steps}",
            f"attraction:      {self.settings.attraction} kJ/mol"
            + (f" (width {self.settings.attraction_width} nm)" if self.settings.attraction else ""),
            f"seed:            {self.seed}",
        ]
        return [line for line in lines if line is not None]


# ---------------------------------------------------------------------------
# Reading the record back (metadata, notebook)
# ---------------------------------------------------------------------------

def summarize_event_rows(rows: list[dict[str, Any]], n_chains: int | None = None) -> dict[str, Any]:
    """Surface-bond statistics from crosslink_events.csv rows (any source)."""
    surface = [r for r in rows if r.get("kind") == "surface"]
    initial = [r for r in surface if r.get("origin") == "initial"]
    per_chain: dict[int, int] = {}
    for r in surface:
        try:
            chain = int(r["chain_i"])
        except (KeyError, TypeError, ValueError):
            continue
        per_chain[chain] = per_chain.get(chain, 0) + 1
    counts = list(per_chain.values())
    if n_chains:
        counts += [0] * max(0, n_chains - len(per_chain))
    hist = np.bincount(np.array(counts, dtype=int), minlength=4) if counts else np.zeros(4, dtype=int)
    return {
        "surface_bonds": len(surface),
        "surface_bonds_initial": len(initial),
        "surface_bonds_run": len(surface) - len(initial),
        "chains_with_1": int(hist[1]),
        "chains_with_2": int(hist[2]),
        "chains_with_3plus": int(hist[3:].sum()),
        "chains_with_0": int(hist[0]),
        # time_ps is step x 0.01 rounded to 6 places; round back to the integer step.
        "last_run_bond_step": max((int(round(float(r["time_ps"]) / DT_PS)) for r in surface
                                   if r.get("origin") == "run"), default=None),
    }
