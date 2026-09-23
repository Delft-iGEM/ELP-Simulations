"""Kirkwood-Riseman on cases with a known answer."""
import numpy as np, pytest
from tools.hydrodynamic import rh_of_frame, rg_of_frame, unwrap_chain

def test_two_beads():
    # N=2: 1/Rh = (1/4)*2*(1/d) = 1/(2d)  ->  Rh = 2d
    p = np.array([[0.,0.,0.],[3.,0.,0.]])
    assert rh_of_frame(p) == pytest.approx(6.0)

def test_equilateral_triangle():
    # N=3, all pairs at d: 1/Rh = (1/9)*6*(1/d) = 2/(3d) -> Rh = 1.5d
    d = 2.0
    p = np.array([[0,0,0],[d,0,0],[d/2, d*np.sqrt(3)/2, 0]], dtype=float)
    assert rh_of_frame(p) == pytest.approx(1.5*d)

def test_rh_scales_with_size():
    rng = np.random.default_rng(0)
    p = rng.normal(size=(50,3))
    assert rh_of_frame(p*3.0) == pytest.approx(3.0*rh_of_frame(p))

def test_rg_known():
    # two beads distance d apart: Rg = d/2
    p = np.array([[0.,0.,0.],[4.,0.,0.]])
    assert rg_of_frame(p) == pytest.approx(2.0)

def test_unwrap_restores_a_split_chain():
    box = np.array([10.,10.,10.])
    whole = np.array([[4.6,5,5],[5.0,5,5],[5.4,5,5],[5.8,5,5]])
    wrapped = np.mod(whole + np.array([5.0,0,0]), box)   # shift so it straddles
    out = unwrap_chain(wrapped, box)
    # bond lengths must survive the wrap
    assert np.allclose(np.linalg.norm(np.diff(out,axis=0),axis=1), 0.4)
    assert rh_of_frame(out) == pytest.approx(rh_of_frame(whole))

def test_unwrapping_matters():
    """Without unwrapping a split chain gives a badly wrong Rh."""
    box = np.array([10.,10.,10.])
    whole = np.array([[4.6,5,5],[5.0,5,5],[5.4,5,5],[5.8,5,5]])
    wrapped = np.mod(whole + np.array([5.0,0,0]), box)
    assert rh_of_frame(wrapped) != pytest.approx(rh_of_frame(whole), rel=0.05)
