from scipy.spatial import Voronoi
from scipy.spatial import distance as d
from pysal.cg import Polygon, LineSegment
from pysal.cg import get_segment_point_dist as dist_to_ls
import numpy as np

def maximum_contained_circle(points):
    boundary = Polygon(points)
    voronoi = Voronoi(points)
    within = [boundary.contains_point(pt) for pt in voronoi.vertices]
    ivoronoi = voronoi.vertices[within]
   
    # Unfortunately, pysal.cg.get_polygon_point_dist returns 0 for points
    # internal to the polygon. So, we have to consider each line segment forming
    # the boundary independently from the polygon itself. It'd be nice if we
    # could modify that behavior conditionally. 
    linesegs  = [LineSegment(boundary.parts[0][i], boundary.parts[0][i+1])
                 for i in range(len(boundary.parts[0])-1)]
    
    radius = -np.inf
    center = (np.nan, np.nan)

    # The maximal contained circle is centered on a vertex of the voronoi
    # partition that has the largest distance to nearest side. 
    for pt in ivoronoi:
        closest = np.min([dist_to_ls(ls, pt) for ls in linesegs])
        if closest > radius:
            radius = closest
            center = pt
    return center, radius
