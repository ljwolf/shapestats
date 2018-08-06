from shapely import geometry
from numpy import testing
from ..minbc import minimum_bounding_circle
from ..maxbc import maximum_contained_circle
from ..compactness import _get_pointset
from .test_measures import shape, ATOL
from shapely.geometry import Polygon

pointset = _get_pointset(shape)

def test_minbc():
    radius, center = minimum_bounding_circle(pointset) 
    testing.assert_allclose(radius, .800390, atol=ATOL)
    testing.assert_array_equal((0.625, 0.5), center)
    circ = minimum_bounding_circle(shape)
    assert isinstance(circ, Polygon)


def test_maxbc():
    radius, center = maximum_contained_circle(pointset)
    testing.assert_allclose(radius, .309359, atol=ATOL)
    testing.assert_array_equal((0.4375, 0.5), center)
    circ = minimum_bounding_circle(shape)
    assert isinstance(circ, Polygon)

