from math import pi as PI
from scipy.spatial import ConvexHull
from pysal.cg import get_angle_between, Ray, is_clockwise
import copy 
import numpy as np
import scipy.spatial.distance as dist
from warnings import warn as Warn
from itertools import cycle
import decimal as dec
from warnings import warn as Warn

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
    points = [points[i] for i in ConvexHull(points).vertices]
    points.reverse() #shift from ccw to cw
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
            return circles[lexmax]
        i+=1

def mbc_animation(points, plotname=False):
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
        
    points = [points[i] for i in ConvexHull(points).vertices]
    points.reverse() #shift from ccw to cw
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
            plt.plot([p[0] for p in POINTS], [p[-1] for p in POINTS], 'r')
            plt.plot([p[0] for p in ORIGINAL_POINTS], 
                     [p[-1] for p in ORIGINAL_POINTS], 'r--')
            plt.plot([p[0] for p in points], [p[-1] for p in points], 'b--')
            plt.plot([c[0] for c in candidate], [c[1] for c in candidate], 'bo')
            fig.gca().add_artist(plt.Circle(circles[lexmax][-1], radii[lexmax], fc='w'))
            if angles[lexmax] > PI/2:
                plt.title('FAILURE')
                plt.xlabel('iteration: {}'.format(i))
                removed = points.pop(lexmax-1)
                plt.plot(removed[0], removed[1], 'kx')
                plt.savefig(str(plotname)+'/'+str(i))
                plt.close()
            else:
                plt.title('SUCCESS!')
                plt.xlabel('iteration: {}'.format(i))
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
    import pysal as ps
    import time

    df = ps.pdio.read_files(ps.examples.get_path('columbus.shp'))

    chain = df['geometry'][2]

    ptset = [pt for part in chain.parts for pt in part]

    mbc_animation(ptset, plotname='test')
    st = time.time()
    minimum_bounding_circle(ptset)
    elapsed = time.time() - st
    print('took {} s'.format(elapsed))
