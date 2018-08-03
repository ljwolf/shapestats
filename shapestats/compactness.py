from __future__ import division
from math import pi as _PI
from pysal import cg as _cg
from scipy.spatial import distance as _dst
import numpy as np

from . import _util as _u
from .minbc import minimum_bounding_circle as _mbc
from ._amoments import second_moa

__all__ = ['ipq','iaq','convex_hull','reock', 'nmi', 'moa_ratio',
           'moment_of_inertia', 'flaherty_crumplin_radius','taylor_reflexive',
           'flaherty_crumplin_lw','eig_seitzinger', 
           'polsby_popper', 'schwartzberg']

### ---- Altman's PA/A measures ---- ##

def ipq(chain):
    """
    The Isoperimetric quotient, defined as the ratio of a chain's area to the 
    area of the equi-perimeter circle. 

    Altman's PA_1 measure

    Construction:
    --------------
    let:
    p_d = perimeter of district
    a_d = area of district
    
    a_c = area of the constructed circle
    r = radius of constructed circle

    then the relationship between the constructed radius and the district
    perimeter is:
    p_d = 2 \pi r
    p_d / (2 \pi) = r
    
    meaning the area of the circle can be expressed as:
    a_c = \pi r^2
    a_c = \pi (p_d / (2\pi))^2
    
    implying finally that the IPQ is:

    pp = (a_d) / (a_c) = (a_d) / ((p_d / (2*\pi))^2 * \pi) = (a_d) / (p_d**2 / (4\PI))
    """
    return (4 * _PI * chain.area) / (chain.perimeter**2)

def convex_hull(chain):
    """
    ratio of the convex hull area to the area of the shape itself

    Altman's A_3 measure, from Neimi et al 1991. 
    """
    pointset = [point for part in chain.parts for point in part] 
    chull = _cg.Polygon(_cg.convex_hull(pointset))
    return chain.area / chull.area

def iaq(chain):
    """
    The Isoareal quotient, defined as the ratio of a chain's perimeter to the
    perimeter of the equi-areal circle

    Altman's PA_3 measure, and proportional to the PA_4 measure 
    """
    return (2 * _PI * np.sqrt(chain.area/_PI)) / chain.perimeter

def reock(chain):
    """
    The Reock compactness measure, defined by the ratio of areas between the
    minimum bounding/containing circle of a shape and the shape itself. 

    Measure A1 in Altman's thesis, cited for Frolov (1974), but earlier from Reock
    (1963)
    """
    pointset = [pt for part in chain.parts for pt in part]
    radius, (cx, cy) = _mbc(pointset)
    return chain.area / (_PI * radius ** 2)

def nmi(chain):
    """
    Computes the Normalized Moment of Inertia from Li et al (2013), recognizing
    that it is the relationship between the area of a shape squared divided by
    its second moment of area. 
    """
    return chain.area**2 / (2 * second_moa(chain) * _PI)

def moa_ratio(chain):
    """
    Computes the ratio of the second moment of area (like Li et al (2013)) to
    the moment of area of a circle with the same area. 
    """
    r = chain.perimeter / (2 * _PI)
    return (_PI * .5 * r**4) / second_moa(chain)


## ---- Altman's OS Measures ---- ##

def moment_of_inertia(chain, dmetric=_dst.euclidean):
    """
    Computes the moment of inertia of the chain. 

    This treats each boundary point as a point-mass of 1.

    Thus, for constant unit mass at each boundary point, 
    the MoI of this pointcloud is

    \sum_i d_{i,c}^2

    where c is the centroid of the chain
    
    Altman's OS_1 measure, cited in Boyce and Clark (1964), also used in Weaver
    and Hess (1963).
    """
    pointset = [pt for part in chain.parts for pt in part]
    dists = [dmetric(pt, chain.centroid)**2 for pt in pointset]
    return chain.area / np.sqrt(2 * np.sum(dists))


def flaherty_crumplin_radius(chain):
    """
    The Flaherty & Crumplin (1992) index, OS_3 in Altman's thesis. 
    
    The ratio of the radius of the equi-areal circle to the radius of the MBC
    """
    pointset = [point for part in chain.parts for point in part[:-1]]
    r_eac = np.sqrt(chain.area/_PI)
    r_mbc, _ = _mbc(pointset)
    return r_eac / r_mbc

def taylor_reflexive(chain):
    """
    The Taylor reflexive angle index, measure OS_4 in Altman's Thesis

    (N-R)/(N+R), the difference in number between non-reflexive angles and
    reflexive angles in a polygon, divided by the number of angles in the
    polygon in general.
    """
    angles = _u.get_all_angles(chain)
    R = np.sum((np.array(A) >= 0).sum() for A in angles) #again ^^^
    return (N - R)/(N+R)

## ---- Altman's Length-Width Measures ---- ##

def flaherty_crumplin_lw(chain):
    """
    The Flaherty & Crumplin (1992) length-width measure, stated as measure LW_7
    in Altman's thesis. 

    It is given as the ratio between the minimum and maximum shape diameter. 
    """
    _, minlen, _, maxlen = _u.unique_lw(chain)
    return minlen / maxlen

def eig_seitzinger(chain):
    """
    The Eig & Seitzinger (1981) shape measure, defined as:

    L - W

    Where L is the maximal east-west extent and W is the maximal north-south
    extent. 

    Defined as measure LW_5 in Altman's thesis
    """
    ptset = [pt for part in chain.parts for pt in part]
    xs, ys = [p[0] for p in ptset], [p[1] for p in ptset]
    l = np.max(xs) - np.min(xs)
    w = np.max(ys) - np.min(ys)
    return l - w

## ---- Alternative Names ---- ##

def polsby_popper(chain):
    """
    Alternative name for the Isoperimetric Quotient
    """
    return ipq(chain)

def schwartzberg(chain):
    """
    Alterantive name for the Isoareal Quotient
    """
    return iaq(chain)
