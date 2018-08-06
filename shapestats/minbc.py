from math import pi as PI
from scipy.spatial import ConvexHull
from libpysal.cg import get_angle_between, Ray, is_clockwise
import copy 
import numpy as np
import scipy.spatial.distance as dist
from warnings import warn as Warn
from itertools import cycle
import decimal as dec

not_clockwise = lambda x: not is_clockwise(x)

def minimum_bounding_circle(points):
    """
    Implements Skyum (1990)'s algorithm for the minimum bounding circle in R^2. 

    0. Store points clockwise. 
    1. Find p in S that maximizes angle(prec(p), p, succ(p) THEN radius(prec(p),
    p, succ(p)). This is also called the lexicographic maximum, and is the last
    entry of a list of (radius, angle) in lexicographical order. 
    2a. If angle(prec(p), p, succ(p)) <= 90 degrees, then finish. 
    2b. If not, remove p from set. 
    """
    was_polygon = not isinstance(points, (np.ndarray,list))
    if was_polygon:
        from .compactness import _get_pointset
        points = _get_pointset(points)
    chull = ConvexHull(points)
    points = np.asarray(points)[chull.vertices]
    points = points[::-1] #shift from ccw to cw
    points = list(map(tuple, points))
    POINTS = copy.deepcopy(points)
    removed = []
    i=0
    while True:
        triple = _nples(points)
        angles = [_angle(*next(triple)) for p in points]
        circles = [_circle(*next(triple)) for p in points]
        radii = [np.max(c[0]) for c in circles]
        lexord = np.lexsort((angles, radii))#really weird defaults
        lexmax = lexord[-1] #recall, we're addressing p, previous_p, two_before_p

        #recall the triple generator is in leading-index order.
        candidate = points[lexmax], points[lexmax-1], points[lexmax-2]
        
        if angles[lexmax] > PI/2:
            removed = points.pop(lexmax-1)
        else:
            radius, center = circles[lexmax]
            if was_polygon:
                from shapely import geometry
                return geometry.Point(tuple(center)).buffer(radius)
            return circles[lexmax]
        i+=1

def _mbc_animation(points, plotname=False, buffer_=.2):
    """
    Implements Skyum (1990)'s algorithm for the minimum bounding circle in R^2. 

    0. Store points clockwise. 
    1. Find p in S that maximizes angle(prec(p), p, succ(p) THEN radius(prec(p),
    p, succ(p)). This is also called the lexicographic maximum, and is the last
    entry of a list of (radius, angle) in lexicographical order. 
    2a. If angle(prec(p), p, succ(p)) <= 90 degrees, then finish. 
    2b. If not, remove p from set. 
    """
    if plotname:
        import matplotlib.pyplot as plt
        import seaborn as sb
        import os
        if not os.path.isdir(plotname):
            os.makedirs('./'+str(plotname))
        ORIGINAL_POINTS = points
        
    chull = ConvexHull(points)
    points = np.asarray(points)[chull.vertices]
    points = points[::-1] #shift from ccw to cw
    minx,miny = points.min(axis=0)
    maxx,maxy = points.max(axis=0)
    width, height = np.abs(minx-maxx), np.abs(miny-maxy)
    left, bottom = minx - width*buffer_, miny - height*buffer_
    right, top = maxx + width*buffer_, maxy + height*buffer_
    points = list(map(tuple, points))
    POINTS = copy.deepcopy(points)
    removed = []
    i=0
    while True:
        triple = _nples(points)
        angles = [_angle(*next(triple)) for p in points]
        circles = [_circle(*next(triple)) for p in points]
        radii = [np.max(c[0]) for c in circles]
        lexord = np.lexsort((angles, radii))#really weird defaults
        lexmax = lexord[-1] #recall, we're addressing p, previous_p, two_before_p

        #recall the triple generator is in leading-index order.
        candidate = points[lexmax], points[lexmax-1], points[lexmax-2]
        
        if plotname:
            sb.set_context('talk')
            fig = plt.figure(figsize=(10,10))
            ax = plt.gca()
            ax.set_aspect('equal')
            ax.plot([p[0] for p in POINTS], [p[-1] for p in POINTS], 'r')
            ax.plot([p[0] for p in ORIGINAL_POINTS], 
                    [p[-1] for p in ORIGINAL_POINTS], 'r--')
            ax.plot([p[0] for p in points], [p[-1] for p in points], 'b--')
            ax.plot([c[0] for c in candidate], [c[1] for c in candidate], 'bo')
            circ = plt.Circle(circles[lexmax][-1], radii[lexmax], 
                              fc='w', ec='k')
            ax.add_artist(circ)
            ax.set_xlim(left, right)
            ax.set_ylim(bottom, top)
            if angles[lexmax] > PI/2:
                ax.set_title('FAILURE')
                ax.set_xlabel('iteration: {}'.format(i))
                removed = points.pop(lexmax-1)
                ax.plot(removed[0], removed[1], 'kx')
                plt.savefig(str(plotname)+'/'+str(i))
                plt.close()
            else:
                ax.set_title('SUCCESS!')
                ax.set_xlabel('iteration: {}'.format(i))
                plt.savefig(str(plotname)+'/'+str(i))
                plt.close()
                return circles[lexmax]
        else:
            if angles[lexmax] > PI/2:
                removed = points.pop(lexmax-1)
            else:
                return circles[lexmax]
        i+=1

def _angle(p,q,r):
    """
    compute the positive angle formed by PQR
    """
    p,q,r = list(map(tuple, (p,q,r)))
    return np.abs(get_angle_between(Ray(q,p),Ray(q,r)))

def _nples(l, n=3):
    """
    continuously yield a cycle of adjacent n from list, looping around the
    ending of l.

    Most "current" is on the left. 
    """
    cyc = cycle(l)
    previous = [l[-i] for i in range(1,n)] #need -1, -2, .., -(n-1)
    while cyc:
        current = next(cyc)
        yield [current] + previous
        previous = [current] + previous[:-1] 

def _circle(A,B,C, dmetric=dist.euclidean):
    """
    Returns (radius, (center_x, center_y)) of the circumscribed circle by the
    triangle pqr.

    note, this does not assume that p!=q!=r
    """
    Ax,Ay = dec.Decimal(A[0]), dec.Decimal(A[1])
    Bx,By = dec.Decimal(B[0]), dec.Decimal(B[1])
    Cx,Cy = dec.Decimal(C[0]), dec.Decimal(C[1])
    if np.array_equal([Ax,Ay], [Bx,By]) or np.array_equal([Bx,By],[Cx,Cy]):
        Warn('Duplicate neighboring point detected!')
    elif np.allclose(_angle(A,B,C), 0):
        #Warn('angle close to zero')
        radii = dist.euclidean(A,B)/2.
        center_x = float(Ax + Bx)/2.
        center_y = float(Ay + By)/2.
    else:
        try:
            D = 2*(Ax*(By - Cy) + Bx*(Cy - Ay) + Cx*(Ay - By))
            center_x = float(((Ax**2 + Ay**2)*(By-Cy) 
                            + (Bx**2 + By**2)*(Cy-Ay) 
                            + (Cx**2 + Cy**2)*(Ay-By)) / D)
            center_y = float(((Ax**2 + Ay**2)*(Cx-Bx) 
                            + (Bx**2 + By**2)*(Ax-Cx) 
                            + (Cx**2 + Cy**2)*(Bx-Ax)) / D)
            radii = np.max([dmetric((center_x, center_y), pt) for pt in [A,B,C]])
        except dec.InvalidOperation:
            center_x = center_y = radii = -np.inf
    return radii, (center_x, center_y)

if __name__ == '__main__':
    import libpysal 
    from libpysal.weights._contW_lists import _get_verts as _get_pointset
    import geopandas

    df = geopandas.read_file(libpysal.examples.get_path('columbus.shp'))
    ix = np.random.randint(0,len(df.geometry))
    #ix = 5
    ptset = _get_pointset(df.geometry[ix])

    _mbc_animation(ptset, plotname='test', buffer_=.66)
    print(ix)
