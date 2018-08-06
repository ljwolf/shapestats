import libpysal.cg as cg
import scipy.spatial.distance as d
import numpy as np
from libpysal.weights._contW_lists import _get_verts as _get_pointset

def all_angles(chain):
    """
    Construct all angles for all parts of a polygon 
    """
    chain = cg.shapes.asShape(chain)
    parts = []
    for part in chain.parts:
        angles = []
        for i in range(len(part)):
            R1 = cg.Ray(part[i-1], part[i-2])
            R2 = cg.Ray(part[i-1], part[i])
            angles.append(cg.get_angle_between(R1, R2))
        parts.append(angles)
    return parts

def pairwise_lw(chain):
    """
    Construct the diameter and width of a polygon, as defined as the longest and
    shortest pairwise distances between a polygon's vertices. 

    Returns
    --------

    Returns the indices of all minima/maxima pairs, as well as their values. 

    That is, if two points are minimal with distance d and three pairs are
    maximal with distance D, this would return:

    ((p1, p2), d, (p1, p2, p3), d)
    
    Thus, you need to be careful with automated unpacking.
    """
    ptset = _get_pointset(chain)
    pwds = d.pdist(ptset)
    sqf = d.squareform(pwds)
    minval = pwds[np.nonzero(pwds)].min() #have to account for zero selfdist
    amin, amax = np.where(sqf == minval), np.where(sqf == pwds.max())
    return amin, sqf[amin], amax, sqf[amax]

def unique_lw(chain):
    """
    Return a unique longest and shortest length for points on a chain's
    diameter. 
    
    This will filter out multiple results from pairwise_lw, but it may be the
    case that other points have at least the minimal and maximal lengths
    returned by thi function. 

    Returns
    -------

    (min_idxs), min_dist, (max_idxs), max_dist

    where min_idxs and max_idxs are the indices of a minimal and maximal
    distance point pair, and min_dist/max_dist are the smallest and largest
    distances between points in the polygon..
    """
    mins, mindists, maxes, maxdists = pairwise_lw(chain)
    assert np.all(mindists == mindists.min())
    mindists = mindists.min() #take smallest value, should be all the same
    mins = mins[0][0], mins[1][0] #grab the first from the index pairs
    assert np.all(maxdists == maxdists.max())
    maxdists = maxdists.max() #take largest value, should all be the same
    maxes = maxes[0][0], maxes[1][0] #grab the first from the index pairs
    return mins, mindists, maxes, maxdists
