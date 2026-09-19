"""Unit tests for the SurfaceBinder side of tools.surface (reactor, state,
outputs). The chain-construction geometry (place_chains / build_chain) is
covered elsewhere.

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

import json

import numpy as np
import pytest

from tools.crosslink import DT_PS, SiteInfo
from tools.surface import (
    STATE_FILENAME,
    SurfaceBinder,
    SurfaceEvent,
    SurfaceSettings,
    check_capacity,
    choose_initial_pins,
    summarize_event_rows,
)

Z_WALL, TETHER = 1.9, 0.6
Z_PIN = Z_WALL + TETHER


def make_sites(n_chains: int, lys_res: list[int], chain_len: int) -> list[SiteInfo]:
    return [SiteInfo(index=c * chain_len + r, chain=c, resid=r + 1, res_in_chain=r,
                     resname="LYS", atomname="CA")
            for c in range(n_chains) for r in lys_res]


def settings(**kw) -> SurfaceSettings:
    base = dict(mode="free", distance=0.8, prob=1.0, max_fraction=1.0, tether=TETHER,
                ramp_steps=500, check_every=1000, start_step=0, seed=7, z_wall=Z_WALL)
    base.update(kw)
    return SurfaceSettings(**base)


def high_positions(n_beads: int) -> np.ndarray:
    """Everything 10 nm above the tether plane: nothing binds."""
    return np.column_stack([np.arange(n_beads) * 0.5, np.zeros(n_beads), np.full(n_beads, Z_PIN + 10.0)])


def make_binder(sites, n_max=None, **kw) -> SurfaceBinder:
    sb = SurfaceBinder(settings(**kw), sites, seed=3)
    sb.n_max = len(sites) if n_max is None else n_max
    return sb


# ---------------------------------------------------------------------------
# capacity (sb.n_max) with non-integer max_fraction x n_sites
# ---------------------------------------------------------------------------

def test_check_capacity_floors_the_cap():
    assert check_capacity(settings(max_fraction=0.5), n_chains=2, n_sites=5) == 2      # 2.5 -> 2
    assert check_capacity(settings(max_fraction=0.7), n_chains=2, n_sites=5) == 3      # 3.5 -> 3
    assert check_capacity(settings(max_fraction=0.6), n_chains=3, n_sites=5) == 3      # 3.0000000004 -> 3
    with pytest.raises(ValueError, match="every one of the 3 chains"):
        check_capacity(settings(max_fraction=0.5), n_chains=3, n_sites=5)             # cap 2 < 3 chains
    with pytest.raises(ValueError, match="No lysines"):
        check_capacity(settings(), n_chains=1, n_sites=0)


def test_choose_initial_pins_one_per_chain_and_preattached_topup():
    chain_of_site = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    rng = np.random.default_rng(0)
    pins = choose_initial_pins(settings(), chain_of_site, 3, rng)
    assert len(pins) == 3 and sorted(np.array(chain_of_site)[pins].tolist()) == [0, 1, 2]
    s = settings(mode="preattached", preattached_fraction=0.6, max_fraction=0.7)
    pins = choose_initial_pins(s, chain_of_site, 3, np.random.default_rng(0))
    assert len(pins) == round(0.6 * 9) == 5
    assert len(set(np.array(chain_of_site)[pins].tolist())) == 3
    # the preattached count is capped by the (floored) max
    s = settings(mode="preattached", preattached_fraction=0.5, max_fraction=0.5)
    pins = choose_initial_pins(s, chain_of_site, 3, np.random.default_rng(0))
    assert len(pins) == 4 == check_capacity(s, 3, 9)            # round(4.5)=4, floor(4.5)=4


# ---------------------------------------------------------------------------
# react()
# ---------------------------------------------------------------------------

def test_react_binds_within_distance_nearest_first_and_pins_in_plane():
    sites = make_sites(2, [3, 8], 10)               # beads 3, 8, 13, 18
    sb = make_binder(sites, n_max=1)
    pos = high_positions(20)
    pos[3] = [2.0, 3.0, Z_PIN + 0.5]                # 0.5 nm gap
    pos[13] = [4.0, 1.0, Z_PIN - 0.2]               # 0.2 nm gap -> nearest, wins the single slot
    pos[8] = [6.0, 6.0, Z_PIN + 0.9]                # outside 0.8
    assert sb.react(pos, step=1000) is True
    assert sb.pinned.tolist() == [False, False, True, False]
    assert sb.pin_xyz[2].tolist() == [4.0, 1.0, Z_PIN]      # x, y where it was; z on the plane
    assert sb.cap_reached_step == 1000
    e = sb.events[0]
    assert e.step == 1000 and e.site == 2 and e.distance == pytest.approx(0.2) and not e.initial
    # cap reached: nothing more binds, even though bead 3 still qualifies
    assert sb.react(pos, step=2000) is False


def test_react_respects_blocked_and_dynamic_flag():
    sites = make_sites(1, [3, 8], 10)
    pos = high_positions(10)
    pos[3, 2] = pos[8, 2] = Z_PIN
    sb = make_binder(sites)
    assert sb.react(pos, 1000, blocked=np.array([True, False])) is True
    assert sb.pinned.tolist() == [False, True]
    # preattached without dynamic: never binds during the run
    sb2 = make_binder(sites, mode="preattached", preattached_fraction=0.5, max_fraction=1.0)
    assert not sb2.settings.is_dynamic
    assert sb2.react(pos, 1000) is False and sb2.n_pinned == 0
    sb3 = make_binder(sites, mode="preattached", preattached_fraction=0.5, max_fraction=1.0, dynamic=True)
    assert sb3.react(pos, 1000) is True


def test_prob_draws_are_deterministic_and_reseeded_on_load(tmp_path):
    sites = make_sites(4, [1, 3, 5, 7], 8)
    pos = high_positions(32)
    pos[[s.index for s in sites], 2] = Z_PIN + 0.1

    def run(sb, steps):
        for st in steps:
            sb.react(pos, st)
        return [e.to_dict() for e in sb.events]

    a = run(make_binder(sites, prob=0.3), range(1000, 6000, 1000))
    b = run(make_binder(sites, prob=0.3), range(1000, 6000, 1000))
    assert a == b and 0 < len(a) < 16

    sb = make_binder(sites, prob=0.3)
    run(sb, range(1000, 3000, 1000))
    sb.save_state(tmp_path, 2000)
    tails = []
    for _ in range(2):
        r = make_binder(sites, prob=0.3)
        assert r.load_state(tmp_path) == 2000
        tails.append(run(r, range(3000, 8000, 1000)))
    assert tails[0] == tails[1]


# ---------------------------------------------------------------------------
# state file
# ---------------------------------------------------------------------------

def test_save_load_roundtrip_and_consistency_checks(tmp_path):
    sites = make_sites(2, [3, 8], 10)
    sb = make_binder(sites, n_max=3)
    sb.pinned[0] = True
    sb.pin_xyz[0] = [1.0, 2.0, Z_PIN]
    sb.events.append(SurfaceEvent(0, 0, (1.0, 2.0, Z_PIN), 0.0, True))
    pos = high_positions(20)
    pos[13] = [4.0, 1.0, Z_PIN - 0.2]
    sb.react(pos, 1000)
    sb.ramping[2] = 1000
    sb.save_state(tmp_path, 1200)

    r = make_binder(sites, n_max=0)
    assert r.load_state(tmp_path) == 1200
    assert np.array_equal(r.pinned, sb.pinned) and np.allclose(r.pin_xyz, sb.pin_xyz)
    assert r.n_max == 3 and r.ramping == {2: 1000} and r.resume_step == 1200
    assert [e.to_dict() for e in r.events] == [e.to_dict() for e in sb.events]
    assert r.cap_reached_step is None

    payload = json.loads((tmp_path / STATE_FILENAME).read_text())
    payload["events"].pop()
    (tmp_path / STATE_FILENAME).write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="inconsistent"):
        make_binder(sites).load_state(tmp_path)

    sb.save_state(tmp_path, 1200)
    with pytest.raises(ValueError, match="different set of lysines"):
        make_binder(make_sites(2, [3, 7], 10)).load_state(tmp_path)


def test_load_state_recomputes_missing_n_max(tmp_path):
    """An n_max of 0 would read as 'cap reached' and stop all binding."""
    sites = make_sites(2, [3, 8], 10)
    sb = make_binder(sites, n_max=2, max_fraction=0.5)
    sb.save_state(tmp_path, 0)
    payload = json.loads((tmp_path / STATE_FILENAME).read_text())
    del payload["n_max"]
    (tmp_path / STATE_FILENAME).write_text(json.dumps(payload))
    r = make_binder(sites, n_max=0, max_fraction=0.5)
    r.load_state(tmp_path)
    assert r.n_max == 2


# ---------------------------------------------------------------------------
# force and ramp (no Context)
# ---------------------------------------------------------------------------

def _k(force, s: int) -> float:
    return force.getParticleParameters(s)[1][0]


def test_add_force_and_ramp_resume(tmp_path):
    import openmm

    sites = make_sites(2, [3, 8], 10)
    sb = make_binder(sites, k=1000.0, ramp_steps=500)
    sb.pinned[0] = True
    sb.pin_xyz[0] = [1.0, 2.0, Z_PIN]
    sb.events.append(SurfaceEvent(0, 0, (1.0, 2.0, Z_PIN), 0.0, True))
    _sys1 = openmm.System()   # kept alive: addForce takes ownership
    force = sb.add_force(_sys1)
    assert force.getNumParticles() == 4
    assert _k(force, 0) == 1000.0 and all(_k(force, s) == 0.0 for s in (1, 2, 3))
    assert force.getParticleParameters(0)[0] == 3                  # bead index of site 0

    pos = high_positions(20)
    pos[13] = [4.0, 1.0, Z_PIN - 0.2]
    sb.react(pos, 1000)
    assert _k(force, 2) == 0.0 and sb.ramping == {2: 1000}
    ks = []
    for st in range(1100, 1600, 100):
        sb.advance_ramps(st)
        ks.append(_k(force, 2))
    assert ks == pytest.approx([200, 400, 600, 800, 1000]) and sb.ramping == {}

    # restart caught mid-ramp
    sb2 = make_binder(sites, k=1000.0, ramp_steps=500)
    pos[3] = [5.0, 5.0, Z_PIN]
    sb2.pinned[0] = True
    sb2.pin_xyz[0] = [1.0, 2.0, Z_PIN]
    sb2.events.append(SurfaceEvent(0, 0, (1.0, 2.0, Z_PIN), 0.0, True))
    # The ramp start step is only recorded when a force exists (post-hoc replay
    # has none and nothing to ramp), so a binder standing in for a live run has
    # to have one — otherwise nothing is mid-ramp and the restart below cannot
    # test resuming one.
    _sys_live = openmm.System()   # kept alive: addForce takes ownership
    sb2.add_force(_sys_live)
    sb2.react(pos, 1000)
    sb2.save_state(tmp_path, 1300)
    r = make_binder(sites, k=1000.0, ramp_steps=500)
    r.load_state(tmp_path)
    _sys2 = openmm.System()   # kept alive: addForce takes ownership
    f2 = r.add_force(_sys2)
    assert _k(f2, 2) == pytest.approx(600.0)                        # 300/500 of the way
    assert _k(f2, 0) == 1000.0


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------

def test_event_rows_frame_convention_and_summary():
    sites = make_sites(2, [3, 8], 10)
    sb = make_binder(sites)
    sb.pinned[0] = True
    sb.pin_xyz[0] = [1.0, 2.0, Z_PIN]
    sb.events.append(SurfaceEvent(0, 0, (1.0, 2.0, Z_PIN), 0.0, True))
    pos = high_positions(20)
    pos[13] = [4.0, 1.0, Z_PIN]
    sb.react(pos, 2000)
    rows = sb.event_rows(n_save=1000)
    assert [r["frame"] for r in rows] == [0, 1]           # step 0 -> frame 0; step 2000 -> frame 1
    assert [r["time_ps"] for r in rows] == [0.0, 20.0]
    assert [r["origin"] for r in rows] == ["initial", "run"]
    assert rows[1]["pymol_i"] == 14 and rows[1]["kind"] == "surface" and rows[1]["resid_j"] == ""
    assert rows[1]["z"] == round(Z_PIN, 4)
    lines = "\n".join(sb.summary_lines(2))
    assert "surface-bonded:        2" in lines and "chains with 1 bond:    2" in lines
    stats = summarize_event_rows(rows, n_chains=2)
    assert stats["surface_bonds"] == 2 and stats["surface_bonds_run"] == 1
    assert stats["last_run_bond_step"] == 2000 and isinstance(stats["last_run_bond_step"], int)
    assert stats["chains_with_1"] == 2 and stats["chains_with_0"] == 0


def test_settings_validation_and_defaults():
    with pytest.raises(ValueError, match="brush"):
        SurfaceSettings(mode="brush", z_wall=Z_WALL)
    with pytest.raises(ValueError, match="preattached_fraction"):
        SurfaceSettings(mode="preattached", z_wall=Z_WALL)
    with pytest.raises(ValueError, match="exceeds"):
        SurfaceSettings(mode="preattached", preattached_fraction=0.6, max_fraction=0.5, z_wall=Z_WALL)
    s = SurfaceSettings.from_dict({"mode": "free", "z_wall": Z_WALL, "unknown_key": 1})
    assert s.is_dynamic and s.z_pin == pytest.approx(Z_PIN)
    assert round(2000 * DT_PS, 6) == 20.0
