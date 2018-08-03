import numpy as np
from libpysal.cg.shapes import Polygon, asShape


def second_moa(chain):
    """
    Using equation listed on en.wikipedia.org/Second_Moment_of_area, the second
    moment of area is actually the cross-moment of area between the X and Y
    dimensions:

    I_xy = (1/24)\sum^{i=N}^{i=1} (x_iy_{i+1} + 2*x_iy_i + 2*x_iy_{i+1} +
    x_{i+1}y_i)(x_iy_i - x_{i+1}y_i)

    where x_i, y_i is the current point and x_{i+1}, y_{i+1} is the next point,
    and where x_{n+1} = x_1, y_{n+1} = 1.

    This relation is known as the:
    - second moment of area
    - moment of inertia of plane area
    - area moment of inertia
    - second area moment

    and is *not* the mass moment of inertia, a property of the distribution of
    mass around a shape.
    """
    chain = asShape(chain)
    outer_I = 0
    if chain.holes == [[]]:
        hole_I = 0
    else:
        hole_I = np.sum([second_moa(Polygon(hole)) for hole in chain.holes])
    for part in chain.parts:
        part = list(reversed(part))  # using equation for ccw MoA
        moment = 0
        for i, this in enumerate(part[:-1]):  # iterate pairwise, so drop last
            _next = part[i + 1]
            thisx, thisy = this
            nextx, nexty = _next
            first = (thisx * nexty + 2 * thisx * thisy + 2 * nextx * nexty + nextx * thisy)
            second = (thisx * nexty - nextx * thisy)
            moment += first * second
        outer_I += moment
    return (1 / float(24)) * (np.abs(outer_I) - np.abs(hole_I))
