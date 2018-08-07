from shapely import geometry
from ..compactness import *
from numpy import testing

shape = geometry.Polygon([(0,0),
                          (.25,.25),
                          (0,.5),
                          (.25,.75),
                          (0,1),
                          (1.25,1),
                          (.75,.5),
                          (1.25,0)])

ATOL = .001

def test_boundary_amplitude():
    observed = boundary_amplitude(shape)
    testing.assert_allclose(observed, .844527, atol=ATOL)

def test_convex_hull():
    observed = convex_hull(shape)
    testing.assert_allclose(observed, .7, atol=ATOL)

def test_eig_seitzinger():
    observed = eig_seitzinger(shape)
    testing.assert_allclose(observed, .25, atol=ATOL)

def test_flaherty_crumplin():
    observed = flaherty_crumplin_lw(shape)
    testing.assert_allclose(observed, .220863, atol=ATOL)
    observed = flaherty_crumplin_radius(shape)
    testing.assert_allclose(observed, .659366, atol=ATOL)

def test_iaq():
    observed = iaq(shape)
    testing.assert_allclose(observed, .622314, atol=ATOL)
    observed2 = schwartzberg(shape)
    testing.assert_allclose(observed, observed2, atol=ATOL)

def test_ipq():
    observed = ipq(shape)
    testing.assert_allclose(observed, .387275, atol=ATOL)
    observed2 = polsby_popper(shape)
    testing.assert_allclose(observed, observed2, atol=ATOL)

def test_moa():
    observed = moa_ratio(shape)
    testing.assert_allclose(observed, 3.249799, atol=ATOL)

def test_moment_of_interia():
    observed = moment_of_inertia(shape)
    testing.assert_allclose(observed, .315715, atol=ATOL)

def test_nmi():
    observed = nmi(shape)
    testing.assert_allclose(observed, .487412, atol=ATOL)

def test_reock():
    observed = reock(shape)
    testing.assert_allclose(observed, .434764, atol=ATOL)

def test_taylor():
    observed = taylor_reflexive(shape)
    testing.assert_allclose(observed, .25, atol=ATOL)

