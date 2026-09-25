"""SASA on geometries with an analytic answer."""
import numpy as np, pytest
from tools.sasa import (_sphere_points, sasa_of_beads, find_motif,
                        sigmas_by_resname, PROBE_NM)


def test_isolated_bead_is_a_full_sphere():
    r = np.array([0.3])
    pos = np.array([[0., 0., 0.]])
    got = sasa_of_beads(pos, r)[0]
    assert got == pytest.approx(4 * np.pi * (0.3 + PROBE_NM) ** 2, rel=1e-12)


def test_far_apart_beads_do_not_occlude():
    r = np.array([0.3, 0.3])
    pos = np.array([[0., 0., 0.], [50., 0., 0.]])
    full = 4 * np.pi * (0.3 + PROBE_NM) ** 2
    assert sasa_of_beads(pos, r) == pytest.approx([full, full], rel=1e-12)


def test_fully_buried_bead_has_zero_area():
    """A small bead inside a much larger one is completely occluded."""
    r = np.array([0.05, 2.0])
    pos = np.array([[0., 0., 0.], [0., 0., 0.]])
    assert sasa_of_beads(pos, r, which=[0])[0] == pytest.approx(0.0)


def test_touching_pair_loses_area_symmetrically():
    r = np.array([0.3, 0.3])
    g = 0.3 + PROBE_NM
    pos = np.array([[0., 0., 0.], [1.2 * g, 0., 0.]])   # overlapping grown spheres
    a = sasa_of_beads(pos, r)
    full = 4 * np.pi * g ** 2
    assert a[0] == pytest.approx(a[1], rel=1e-6)        # symmetry
    assert 0 < a[0] < full                              # some, not all, is lost


def test_periodic_images_occlude():
    """Across a small box a neighbour's image must still block the surface."""
    r = np.array([0.3, 0.3])
    box = np.array([2.0, 2.0, 2.0])
    pos = np.array([[0.05, 1.0, 1.0], [1.95, 1.0, 1.0]])   # 0.1 nm apart through the wall
    free = sasa_of_beads(pos, r, box=None)[0]
    wrapped = sasa_of_beads(pos, r, box=box)[0]
    assert wrapped < free * 0.9


def test_sphere_points_are_on_the_unit_sphere():
    p = _sphere_points(500)
    assert np.allclose(np.linalg.norm(p, axis=1), 1.0)
    assert abs(p.mean(axis=0)).max() < 0.05          # roughly isotropic


def test_find_motif():
    assert find_motif("AAGRGDSAA", "GRGDS") == [(2, 7)]
    assert find_motif("XX", "GRGDS") == []
    assert len(find_motif("GRGDSGRGDS", "GRGDS")) == 2


def test_tagged_residues_share_their_parent_sigma():
    """Z (anchor valine) and X (tagged lysine) must not conflict with V and K."""
    s = sigmas_by_resname()          # raises if they disagree
    assert s["VAL"] > 0 and s["LYS"] > 0
