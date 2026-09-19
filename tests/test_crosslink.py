"""Unit tests for tools.crosslink — pure-geometry reactor, no OpenMM Context.

Run from the worktree root:
    /home/tnartey/ELP-Simulations/.venv/bin/python -m pytest tests -q
"""

# OpenMM's System.addForce() takes ownership of the force: a temporary
# System is garbage-collected straight away and deletes the force with it,
# leaving every later getParticleParameters/getBondParameters call on a
# dangling pointer (it raises std::bad_array_new_length). Production code
# passes the long-lived sim.system, so this only ever bites tests — keep a
# reference alive here.

from __future__ import annotations

import csv
import json
import math

import numpy as np
import pytest

from tools.crosslink import (
    CSV_COLUMNS,
    DT_PS,
    RAMP_UPDATES,
    STATE_FILENAME,
    Crosslinker,
    CrosslinkSettings,
    SiteInfo,
    _chunk_sizes,
    frame_of_step,
    min_image_distances,
    post_hoc_events,
    write_event_record,
)
from tools.surface import SurfaceBinder, SurfaceSettings

BOX = np.array([20.0, 20.0, 40.0])


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def make_sites(n_chains: int, lys_res: list[int], chain_len: int) -> list[SiteInfo]:
    """Sites for `n_chains` identical chains with lysines at residues `lys_res`
    (0-based within the chain). Bead index = chain * chain_len + res."""
    sites = []
    for c in range(n_chains):
        for r in lys_res:
            sites.append(SiteInfo(index=c * chain_len + r, chain=c, resid=r + 1,
                                  res_in_chain=r, resname="LYS", atomname="CA"))
    return sites


def far_positions(n_beads: int, seed: int = 0) -> np.ndarray:
    """Every bead > 3 nm from every other (a 4 nm grid), nothing reacts."""
    rng = np.random.default_rng(seed)
    side = math.ceil(n_beads ** (1 / 3))
    grid = np.array([[i, j, k] for i in range(side) for j in range(side) for k in range(side)],
                    dtype=float)[:n_beads] * 4.0 + 1.0
    return grid + rng.uniform(-0.2, 0.2, grid.shape)


def settings(**kw) -> CrosslinkSettings:
    base = dict(distance=0.8, valence=1, prob=1.0, check_every=1000, k=2000.0, r0=0.6,
                ramp_steps=500, start_step=0, seed=1, min_span=3)
    base.update(kw)
    return CrosslinkSettings(**base)


def surface_settings(**kw) -> SurfaceSettings:
    base = dict(mode="free", distance=0.8, prob=1.0, max_fraction=1.0, tether=0.6,
                ramp_steps=500, check_every=1000, start_step=0, seed=7, z_wall=1.9)
    base.update(kw)
    return SurfaceSettings(**base)


def bonded_slots(xl: Crosslinker) -> set[int]:
    return set(xl.pair_slot_a[xl.reacted].tolist()) | set(xl.pair_slot_b[xl.reacted].tolist())


# ---------------------------------------------------------------------------
# (e) min_span
# ---------------------------------------------------------------------------

def test_min_span_drops_intra_chain_neighbours_only():
    sites = make_sites(2, [5, 6, 7, 20], 30)     # 8 sites, C(8,2) = 28 all pairs
    xl = Crosslinker(settings(min_span=3), sites, seed=0)
    # per chain the pairs (5,6), (6,7), (5,7) are closer than 3 residues -> 3 x 2 chains
    assert xl.n_pairs == 28 - 6
    assert xl.n_pairs_dropped_by_span == 6
    for a, b in zip(xl.pair_slot_a, xl.pair_slot_b):
        sa, sb = sites[a], sites[b]
        if sa.chain == sb.chain:
            assert abs(sa.res_in_chain - sb.res_in_chain) >= 3
    # inter-chain pairs are all present, including (chain0 res5, chain1 res6)
    inter = sum(1 for a, b in zip(xl.pair_slot_a, xl.pair_slot_b) if sites[a].chain != sites[b].chain)
    assert inter == 16

    xl0 = Crosslinker(settings(min_span=0), sites, seed=0)
    assert xl0.n_pairs == 28 and xl0.n_pairs_dropped_by_span == 0


def test_adjacent_lysines_cannot_crosslink_by_backbone_alone():
    """...KK... on one chain: the CA beads are 0.38 nm apart forever, inside any
    reaction distance. With min_span they are not even candidates."""
    sites = make_sites(1, [10, 11, 12], 30)
    pos = far_positions(30)
    pos[10] = [5.0, 5.0, 5.0]
    pos[11] = pos[10] + [0.38, 0, 0]
    pos[12] = pos[11] + [0.38, 0, 0]
    xl = Crosslinker(settings(min_span=3), sites, seed=0)
    assert xl.n_pairs == 0
    assert xl.react(pos, BOX, step=1000) is False
    assert xl.events == []
    # and the pre-min_span behaviour, for contrast
    xl0 = Crosslinker(settings(min_span=0), sites, seed=0)
    assert xl0.react(pos, BOX, step=1000) is True
    assert len(xl0.events) == 1 and xl0.events[0].distance == pytest.approx(0.38)


def test_settings_from_dict_tolerates_missing_min_span_and_roundtrips():
    d = settings().to_dict()
    del d["min_span"]
    s = CrosslinkSettings.from_dict(d)
    assert s.min_span == 3
    assert CrosslinkSettings.from_dict(settings(min_span=5).to_dict()).min_span == 5
    with pytest.raises(ValueError):
        settings(min_span=-1)


def test_state_file_without_hash_loads_only_if_no_pairs_dropped(tmp_path):
    # no intra pairs within 3 residues: pair list identical to the all-pairs one
    sites = make_sites(2, [3, 48], 60)
    xl = Crosslinker(settings(), sites, seed=0)
    pos = far_positions(120)
    pos[3] = [5, 5, 5]
    pos[63] = [5.5, 5, 5]   # chain 1 res 3
    assert xl.react(pos, BOX, 1000)
    xl.save_state(tmp_path, 1000)
    payload = json.loads((tmp_path / STATE_FILENAME).read_text())
    del payload["pairs_hash"]
    del payload["n_pairs"]
    del payload["settings"]["min_span"]
    (tmp_path / STATE_FILENAME).write_text(json.dumps(payload))
    fresh = Crosslinker(settings(), sites, seed=0)
    assert fresh.load_state(tmp_path) == 1000
    assert fresh.reacted.sum() == 1

    # now a system where min_span DID drop pairs: an old (all-pairs) file must be refused
    sites2 = make_sites(2, [3, 4, 48], 60)
    xl2 = Crosslinker(settings(min_span=0), sites2, seed=0)
    xl2.save_state(tmp_path, 500)
    payload = json.loads((tmp_path / STATE_FILENAME).read_text())
    del payload["pairs_hash"]
    del payload["settings"]["min_span"]
    (tmp_path / STATE_FILENAME).write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="different bond list"):
        Crosslinker(settings(min_span=3), sites2, seed=0).load_state(tmp_path)


def test_state_file_with_other_min_span_is_refused(tmp_path):
    sites = make_sites(2, [3, 4, 48], 60)
    Crosslinker(settings(min_span=0), sites, seed=0).save_state(tmp_path, 100)
    with pytest.raises(ValueError, match="different bond list"):
        Crosslinker(settings(min_span=3), sites, seed=0).load_state(tmp_path)
    # same min_span: fine
    assert Crosslinker(settings(min_span=0), sites, seed=0).load_state(tmp_path) == 100


# ---------------------------------------------------------------------------
# (a) mutual exclusion surface <-> crosslink
# ---------------------------------------------------------------------------

def _live_step(xl: Crosslinker, sb: SurfaceBinder, pos: np.ndarray, step: int,
               do_sb: bool, do_xl: bool) -> None:
    """The reaction part of run_reactive's loop body, in its order: surface
    first with blocked = (xl.used > 0), then crosslinks with blocked = sb.pinned."""
    if do_sb:
        sb.react(pos, step, blocked=(xl.used > 0))
    if do_xl:
        xl.react(pos, BOX, step, blocked=sb.pinned)


def assert_exclusive(xl: Crosslinker, sb: SurfaceBinder) -> None:
    both = sb.pinned & (xl.used > 0)
    assert not both.any(), f"lysine slots both surface-pinned and K-K bonded: {np.flatnonzero(both)}"
    # bookkeeping consistency
    counts = np.bincount(xl.pair_slot_a[xl.reacted], minlength=len(xl.used)) + \
        np.bincount(xl.pair_slot_b[xl.reacted], minlength=len(xl.used))
    assert np.array_equal(counts, xl.used)
    assert (xl.used <= xl.settings.valence).all()
    assert len(sb.events) == sb.n_pinned


def test_surface_pinned_lysine_never_crosslinks_same_check():
    sites = make_sites(2, [3, 48], 60)
    z_pin = 1.9 + 0.6
    pos = far_positions(120)
    # chain0 res3 and chain1 res3 are 0.5 nm apart AND both on the tether plane
    pos[3] = [5.0, 5.0, z_pin]
    pos[63] = [5.5, 5.0, z_pin]
    xl = Crosslinker(settings(), sites, seed=0)
    sb = SurfaceBinder(surface_settings(), sites, seed=0)
    sb.n_max = 4
    _live_step(xl, sb, pos, 1000, do_sb=True, do_xl=True)
    assert sb.pinned[0] and sb.pinned[2]          # both bound the surface
    assert xl.reacted.sum() == 0                  # ...so neither could crosslink
    assert_exclusive(xl, sb)


def test_kk_bonded_lysine_never_pins_even_with_spare_valence():
    """valence 2: used = 1 is not saturated, but the lysine is still bonded and
    must be blocked for the surface (blocked = used > 0, not used >= valence)."""
    sites = make_sites(2, [3, 48], 60)
    z_pin = 2.5
    pos = far_positions(120)
    pos[3] = [5.0, 5.0, 10.0]
    pos[63] = [5.5, 5.0, 10.0]
    xl = Crosslinker(settings(valence=2), sites, seed=0)
    sb = SurfaceBinder(surface_settings(), sites, seed=0)
    sb.n_max = 4
    _live_step(xl, sb, pos, 1000, do_sb=True, do_xl=True)   # xl only fires (far from plane)
    assert xl.used[0] == 1 and xl.used[2] == 1 and xl.used[0] < 2
    # now drop them onto the tether plane: a surface check must not pin them
    pos[3, 2] = z_pin
    pos[63, 2] = z_pin
    _live_step(xl, sb, pos, 2000, do_sb=True, do_xl=False)
    assert not sb.pinned[0] and not sb.pinned[2]
    assert_exclusive(xl, sb)


def test_exclusion_holds_under_different_cadences_and_restart(tmp_path):
    """Randomised: sb checks every 500 steps, xl every 1000, positions
    re-drawn each check so lysines wander through both reaction zones; the
    state is saved / reloaded into fresh reactors half way. The invariant
    'pinned & used > 0 is empty' must hold after every check."""
    n_chains, lys = 4, [3, 8, 13, 18]     # 16 sites, spacing 5 >= min_span
    sites = make_sites(n_chains, lys, 20)
    beads = np.array([s.index for s in sites])
    rng = np.random.default_rng(3)
    z_pin = 2.5

    def frame() -> np.ndarray:
        pos = far_positions(80)
        # crowd the lysines into a 2 nm cube that straddles the tether plane
        pos[beads] = rng.uniform([5, 5, z_pin - 1.0], [7, 7, z_pin + 1.0], (len(beads), 3))
        return pos

    xl = Crosslinker(settings(valence=2, prob=0.7), sites, seed=11)
    sb = SurfaceBinder(surface_settings(prob=0.6, max_fraction=0.75), sites, seed=5)
    sb.n_max = 12
    step = 0
    for _ in range(12):
        step += 500
        _live_step(xl, sb, frame(), step, do_sb=True, do_xl=(step % 1000 == 0))
        assert_exclusive(xl, sb)
    assert sb.n_pinned > 0 and xl.reacted.sum() > 0, "test setup produced no reactions"

    xl.save_state(tmp_path, step)
    sb.save_state(tmp_path, step)
    xl2 = Crosslinker(settings(valence=2, prob=0.7), sites, seed=11)
    sb2 = SurfaceBinder(surface_settings(prob=0.6, max_fraction=0.75), sites, seed=5)
    assert xl2.load_state(tmp_path) == step and sb2.load_state(tmp_path) == step
    assert np.array_equal(xl2.used, xl.used) and np.array_equal(sb2.pinned, sb.pinned)
    assert_exclusive(xl2, sb2)
    for _ in range(12):
        step += 500
        _live_step(xl2, sb2, frame(), step, do_sb=True, do_xl=(step % 1000 == 0))
        assert_exclusive(xl2, sb2)


# ---------------------------------------------------------------------------
# (b) start_step and chunking
# ---------------------------------------------------------------------------

def test_chunk_sizes():
    assert _chunk_sizes(10_000, 1000, 500, ramping=False) == 1000
    assert _chunk_sizes(10_000, 1000, 500, ramping=True) == 500 // RAMP_UPDATES
    assert _chunk_sizes(300, 1000, 500, ramping=False) == 300        # check_every > remaining
    assert _chunk_sizes(10_000, 1000, 0, ramping=True) == 1000       # no ramp -> full chunk
    assert _chunk_sizes(10_000, 1000, 3, ramping=True) == 1          # tiny ramp -> 1-step chunks
    assert _chunk_sizes(10_000, 50, 500, ramping=True) == 50         # check_every already finer
    assert _chunk_sizes(10_000, 1000, 500, False, until=[700, 2000]) == 700   # stop at the next due check
    assert _chunk_sizes(10_000, 1000, 500, False, until=[0]) == 1             # never a zero-step chunk


def test_loop_schedule_never_skips_a_check():
    """Replays run_reactive's scheduling arithmetic (chunk = min(check_every,
    ramp sub-chunk, to-checkpoint, remaining); check when step >= next_check
    and step >= start_step). Every check interval is honoured: the gaps between
    checks are >= check_every and < check_every + sub-chunk, no check lands
    before start_step, and the run ends at exactly `total`."""
    total, check_every, ramp_steps, start_step = 20_000, 1000, 500, 2_500
    # two reactors at different cadences, as the loop supports
    check_xl, check_sb = 1000, 700
    checkpoint_every = max(1, total // 10)
    next_checkpoint = checkpoint_every
    next_xl, next_sb = check_xl, check_sb
    ramping_until = -1
    done = step = 0
    checks_xl, checks_sb, stops = [], [], 0
    while done < total:
        ramping = step < ramping_until
        until = [next_checkpoint - done, next_xl - step, next_sb - step]
        chunk = _chunk_sizes(total - done, min(check_xl, check_sb), ramp_steps, ramping, until)
        step += chunk
        done += chunk
        stops += 1
        if step >= next_sb:
            next_sb = step + check_sb
            if step >= start_step:
                checks_sb.append(step)
        if step >= next_xl:
            next_xl = step + check_xl
            if step >= start_step:
                checks_xl.append(step)
                ramping_until = step + ramp_steps      # pretend every xl check formed a bond
        if done >= next_checkpoint or done >= total:
            next_checkpoint = done + checkpoint_every
    assert done == total == step
    # every check lands exactly on its own multiple: none late, none skipped, none early
    assert checks_xl == [s for s in range(check_xl, total + 1, check_xl) if s >= start_step]
    assert checks_sb == [s for s in range(check_sb, total + 1, check_sb) if s >= start_step]
    # and the loop did not degenerate into single-step chunks before start_step
    assert stops < total // 50


def test_react_itself_has_no_start_step_guard_the_loop_does():
    """Documented: react() is pure geometry; start_step is enforced by the loop
    and by post_hoc_events. A caller of react() must do that itself."""
    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(start_step=10_000), sites, seed=0)
    assert xl.react(pos, BOX, step=1) is True


# ---------------------------------------------------------------------------
# (c) ramp
# ---------------------------------------------------------------------------

def _bond_k(force, b: int) -> float:
    from openmm.unit import kilojoule_per_mole, nanometer

    return force.getBondParameters(b)[3].value_in_unit(kilojoule_per_mole / nanometer ** 2)


def test_ramp_is_monotone_and_reaches_k_without_context():
    import openmm

    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(ramp_steps=500, k=2000.0), sites, seed=0)
    _sys1 = openmm.System()   # kept alive: addForce takes ownership
    xl.add_force(_sys1)
    assert all(_bond_k(xl.force, b) == 0.0 for b in range(xl.n_pairs))
    assert xl.react(pos, BOX, step=1000)
    b = int(np.flatnonzero(xl.reacted)[0])
    assert _bond_k(xl.force, b) == 0.0                 # formed: still soft
    ks = []
    for step in range(1100, 1600, 100):
        assert xl.advance_ramps(step) is True
        ks.append(_bond_k(xl.force, b))
    assert ks == pytest.approx([400, 800, 1200, 1600, 2000])
    assert np.all(np.diff(ks) > 0)
    assert xl.ramping == {}                            # finished -> dropped
    assert xl.advance_ramps(1700) is False
    assert _bond_k(xl.force, b) == 2000.0
    # the other bonds stayed dormant throughout
    assert all(_bond_k(xl.force, c) == 0.0 for c in range(xl.n_pairs) if c != b)


def test_ramp_steps_zero_switches_on_at_formation():
    import openmm

    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(ramp_steps=0, k=1500.0), sites, seed=0)
    _sys2 = openmm.System()   # kept alive: addForce takes ownership
    xl.add_force(_sys2)
    assert xl.react(pos, BOX, step=1000)
    b = int(np.flatnonzero(xl.reacted)[0])
    assert _bond_k(xl.force, b) == 1500.0    # not delayed to the next chunk


def test_restart_mid_ramp_resumes_at_reached_fraction(tmp_path):
    import openmm

    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(ramp_steps=500, k=2000.0), sites, seed=0)
    _sys3 = openmm.System()   # kept alive: addForce takes ownership
    xl.add_force(_sys3)
    xl.react(pos, BOX, step=1000)
    xl.advance_ramps(1200)
    xl.save_state(tmp_path, 1200)
    b = int(np.flatnonzero(xl.reacted)[0])

    fresh = Crosslinker(settings(ramp_steps=500, k=2000.0), sites, seed=0)
    assert fresh.load_state(tmp_path) == 1200
    assert fresh.ramping == {b: 1000}            # formed_at kept as is (the "step - elapsed" no-op)
    _sys4 = openmm.System()   # kept alive: addForce takes ownership
    fresh.add_force(_sys4)
    assert _bond_k(fresh.force, b) == pytest.approx(800.0)     # 200/500 of 2000, not 0, not 2000
    fresh.advance_ramps(1300)
    assert _bond_k(fresh.force, b) == pytest.approx(1200.0)
    fresh.advance_ramps(1500)
    assert _bond_k(fresh.force, b) == 2000.0 and fresh.ramping == {}

    # a fully ramped bond comes back at full k
    xl.advance_ramps(1500)
    xl.save_state(tmp_path, 1500)
    fresh2 = Crosslinker(settings(ramp_steps=500, k=2000.0), sites, seed=0)
    fresh2.load_state(tmp_path)
    _sys5 = openmm.System()   # kept alive: addForce takes ownership
    fresh2.add_force(_sys5)
    assert _bond_k(fresh2.force, b) == 2000.0


# ---------------------------------------------------------------------------
# (d) restart determinism and state consistency
# ---------------------------------------------------------------------------

def _random_frames(beads: np.ndarray, n: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        pos = far_positions(80)
        pos[beads] = rng.uniform([5, 5, 5], [6.5, 6.5, 6.5], (len(beads), 3))
        out.append(pos)
    return out


def test_restart_is_deterministic_and_reseeded(tmp_path):
    sites = make_sites(4, [3, 8, 13, 18], 20)
    beads = np.array([s.index for s in sites])
    frames = _random_frames(beads, 20, seed=42)
    xl = Crosslinker(settings(valence=2, prob=0.4), sites, seed=99)
    for i in range(10):
        xl.react(frames[i], BOX, (i + 1) * 1000)
    xl.save_state(tmp_path, 10_000)
    before = [e.to_dict() for e in xl.events]

    runs = []
    for _ in range(2):
        r = Crosslinker(settings(valence=2, prob=0.4), sites, seed=99)
        assert r.load_state(tmp_path) == 10_000
        assert [e.to_dict() for e in r.events] == before
        assert np.array_equal(r.used, xl.used) and np.array_equal(r.reacted, xl.reacted)
        for i in range(10, 20):
            r.react(frames[i], BOX, (i + 1) * 1000)
        runs.append([e.to_dict() for e in r.events])
    assert runs[0] == runs[1]
    # the restarted stream is not the fresh stream replayed from the start
    a = np.random.default_rng(99).random(5)
    b = np.random.default_rng([99, 10_000]).random(5)
    assert not np.allclose(a, b)


def test_load_state_rejects_inconsistent_used(tmp_path):
    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(), sites, seed=0)
    xl.react(pos, BOX, 1000)
    xl.save_state(tmp_path, 1000)
    payload = json.loads((tmp_path / STATE_FILENAME).read_text())
    payload["used"][0] = 0                      # bond says slot 0 is used; used[] says no
    (tmp_path / STATE_FILENAME).write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="inconsistent"):
        Crosslinker(settings(), sites, seed=0).load_state(tmp_path)


def test_load_state_rejects_other_sites(tmp_path):
    Crosslinker(settings(), make_sites(2, [3, 48], 60), seed=0).save_state(tmp_path, 0)
    with pytest.raises(ValueError, match="different set of reactive sites"):
        Crosslinker(settings(), make_sites(2, [3, 40], 60), seed=0).load_state(tmp_path)


# ---------------------------------------------------------------------------
# (f) frame convention, (g) summary, (i) misc
# ---------------------------------------------------------------------------

def test_frame_of_step_matches_dcdreporter_numbering():
    n_save = 1000
    assert frame_of_step(0, n_save) == 0          # initial bond: before any frame -> frame 0
    assert frame_of_step(1000, n_save) == 0       # first frame written at step n_save
    assert frame_of_step(2000, n_save) == 1
    assert frame_of_step(2400, n_save) == 1       # nearest
    assert frame_of_step(2600, n_save) == 2
    assert frame_of_step(500, 0) == 0


def test_event_rows_and_record_sorting(tmp_path):
    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3], pos[63] = [5, 5, 5], [5.5, 5, 5]
    xl = Crosslinker(settings(), sites, seed=0)
    xl.react(pos, BOX, step=3000)
    sb = SurfaceBinder(surface_settings(), sites, seed=0)
    sb.n_max = 4
    sb.pinned[1] = True
    sb.pin_xyz[1] = [1, 1, 2.5]
    from tools.surface import SurfaceEvent
    sb.events.append(SurfaceEvent(0, 1, (1.0, 1.0, 2.5), 0.0, True))
    pos2 = far_positions(120)
    pos2[108] = [3, 3, 2.5]                        # chain 1 res 48
    sb.react(pos2, step=2000)
    write_event_record(tmp_path, n_save=1000, steps_done=5000, xl=xl, sb=sb, n_chains=2)
    with open(tmp_path / "crosslink_events.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0].keys()) == CSV_COLUMNS
    assert [r["origin"] for r in rows] == ["initial", "run", "run"]
    assert [float(r["time_ps"]) for r in rows] == [0.0, 20.0, 30.0]
    assert [int(r["frame"]) for r in rows] == [0, 1, 2]
    kk = rows[2]
    assert kk["kind"] == "inter" and kk["span"] == "" and int(kk["pymol_i"]) == 4
    summary = (tmp_path / "crosslink_summary.txt").read_text()
    assert "surface-bonded sites:  2 of 4" in summary
    assert "conversion (available): 1.0000" in summary        # 2 available, both bonded
    assert "conversion:            0.5000" in summary


def test_midpoint_is_minimum_image():
    sites = make_sites(2, [3, 48], 60)
    pos = far_positions(120)
    pos[3] = [0.1, 5, 5]
    pos[63] = [BOX[0] - 0.2, 5, 5]                  # 0.3 nm apart across the boundary
    xl = Crosslinker(settings(), sites, seed=0)
    assert xl.react(pos, BOX, 1000)
    mx = xl.events[0].midpoint[0]
    assert mx == pytest.approx(-0.05) or mx == pytest.approx(BOX[0] - 0.05)
    assert xl.events[0].distance == pytest.approx(0.3)


def test_min_image_distances():
    pos = np.array([[0.1, 0, 0], [19.9, 0, 0], [10, 0, 0]])
    d = min_image_distances(pos, np.array([0, 0]), np.array([1, 2]), BOX)
    assert d == pytest.approx([0.2, 9.9])


def test_nearest_first_and_valence():
    sites = make_sites(3, [3], 10)                  # beads 3, 13, 23
    pos = far_positions(30)
    pos[3], pos[13], pos[23] = [5, 5, 5], [5.7, 5, 5], [5.3, 5, 5.0]
    xl = Crosslinker(settings(valence=1), sites, seed=0)
    xl.react(pos, BOX, 1000)
    assert len(xl.events) == 1
    e = xl.events[0]
    assert {e.site_a, e.site_b} == {0, 2}           # 0.3 nm beats 0.7 and 0.4
    assert xl.used.tolist() == [1, 0, 1]
    # the leftover pair (13, 23) is 0.4 apart but 23 is saturated; never reacts
    assert xl.react(pos, BOX, 2000) is False


def test_post_hoc_frame_step_mapping_and_surface_block(tmp_path):
    """Frame k of a DCD is step (k+1)*n_save. A pair in contact only in frame 0
    with start_step = n_save must react (step n_save is not before start_step),
    and a lysine whose surface bond is at time_ps = n_save*DT_PS is blocked in
    that same frame."""
    import mdtraj as md

    n_save = 1000
    top = md.Topology()
    for c in range(2):
        ch = top.add_chain()
        for r in range(6):
            res = top.add_residue("LYS" if r in (1, 5) else "GLY", ch)
            top.add_atom("CA", md.element.carbon, res)
    n = top.n_atoms                                  # 12 beads; lysines 1,5,7,11
    frames = []
    f0 = far_positions(n)
    f0[1], f0[7] = [5, 5, 5], [5.5, 5, 5]            # chain0 res1 <-> chain1 res1 in contact
    f0[5], f0[11] = [8, 8, 8], [8.5, 8, 8]           # chain0 res5 <-> chain1 res5 in contact
    frames.append(f0)
    frames.append(far_positions(n, seed=1))          # nothing in contact afterwards
    xyz = np.array(frames, dtype=np.float32)
    traj = md.Trajectory(xyz, top, unitcell_lengths=np.tile(BOX, (2, 1)),
                         unitcell_angles=np.full((2, 3), 90.0))
    traj[0].save_pdb(str(tmp_path / "top.pdb"))
    traj.save_dcd(str(tmp_path / "t.dcd"))

    s = settings(start_step=n_save, min_span=0)
    xl = post_hoc_events(tmp_path / "t.dcd", tmp_path / "top.pdb", s, n_save=n_save,
                         verbose=False, drop_anchor=False)
    assert len(xl.events) == 2 and all(e.step == n_save for e in xl.events)

    # the same, with chain 0 res 5 (resid 6) surface-bonded at exactly step n_save
    # mdtraj numbers resSeq from 0 across the whole topology when add_residue is
    # given no resSeq, so chain 0's sixth residue is resid 5 (not 6). Take the
    # resid from the site list rather than assuming the numbering.
    blocked_site = next(s for s in xl.sites if s.chain == 0 and s.res_in_chain == 5)
    rows = [{"kind": "surface", "chain_i": "0", "resid_i": str(blocked_site.resid),
             "time_ps": str(round(n_save * DT_PS, 6))}]
    xl2 = post_hoc_events(tmp_path / "t.dcd", tmp_path / "top.pdb", s, n_save=n_save,
                          verbose=False, drop_anchor=False, surface_rows=rows)
    assert len(xl2.events) == 1
    a, b = xl2.sites[xl2.events[0].site_a], xl2.sites[xl2.events[0].site_b]
    assert {a.index, b.index} == {1, 7}
    # start_step one past the only frame with contacts -> nothing
    xl3 = post_hoc_events(tmp_path / "t.dcd", tmp_path / "top.pdb",
                          settings(start_step=n_save + 1, min_span=0), n_save=n_save,
                          verbose=False, drop_anchor=False)
    assert xl3.events == []
