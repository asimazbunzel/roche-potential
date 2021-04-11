'''Collection of classes to study Roche geometry
'''

import numpy as np
import matplotlib.pyplot as plt
try:
    from astropy import units,constans as u,c
    is_astropy = True
except:
    print("astropy not installed, using internal units instead")
    from utils import Constants as u
    is_astropy = False
from utils import a_to_P,P_to_a


class Roche(object):
    '''Roche potential
    
    Contains information on the Roche model for the three-body problem.
    The potential in a binary is determined by the gravitational attraction of the binary pair
    and by the motion of the two stars. With the typical assumptions, and in the binary frame
    (co-rotating with the binary), the Roche potential is:

    \phi = - \dfrac{G M_1}{r_1} - G \dfrac{G M_2}{r_2} - \dfrac{\omega^2 r^{2}_3}{2}

    where G is the gravitational constant, r_1 and r_2 the distances to the center of the stars
    having masses M_1 and M_2, respectively; \omega is the orbital angular velocity which can be
    derived from Kepler's law:

    \omega = \dfrac{2 \pi}{P} = \sqrt{\dfrac{G (M_1 + M_2)}{a^3}}

    and r_3 is the distance to the axis of rotation.

    Parameters
    ----------
    m1 : `float`
       Mass of the primary (more massive star) in Msun.

    m2 : `float`
       Mass of the secondary (less massive star) in Msun.

    a : `float`
       Orbital separation (in Rsun).

    N : `int`
       Number of points in the grid.

    Methods
    -------
    '''

    def __init__(self, m1: float, m2: float, a: float=None, P: float=None, N: int=1000, is_astropy: bool=is_astropy) -> None:

        if P is None and a is None: raise ValueError('either `P` or `a` must be specified')

        if a is None: a = P_to_a(P, m1, m2)
        if P is None: P = a_to_P(a, m1, m2)

        self.m1 = m1 * u.Msun
        self.m2 = m2 * u.Msun
        self.a = a * u.Rsun
        self.P = P * 24e0 * 3600e0
        if is_astropy:
            self.P = self.P * u.s
        self.N = N

        self.x = np.linspace(-2*a, 2*a, self.N, dtype=np.float64)
        self.y = np.linspace(-2*a, 2*a, self.N, dtype=np.float64)

        self.X, self.Y = np.meshgrid(self.x, self.y)

        # this is the dimensionless potential, in the z=0 plane
        self.V = self.Vf(self.X, self.Y)


    def _solve(self, a: float, b: float, x0: float) -> float:
        '''Solver for roots in the potential using first and second derivates
        
        Parameters
        ----------

        Returns
        -------
        '''

        eps = 1e-6

        iters = 10000
        for n in range(iters):
            dVX = self.dVXdx(a, b, x0)
            d2VX = self.d2VXdx2(a, b, x0)

            x1 = x0 - dVX/d2VX
            err = abs(x1 - x0)/abs(x0 + eps)

            if err < eps: break
            x0 = x1

        return x0

    def Vf(self, x: list, y: list) -> np.ndarray:
        '''Roche potential
        
        Provided a grid of points, returns the value of the potential in the mesh
        
        Parameters
        ----------
        x : `array`

        y : `array`

        Returns
        -------
        V : `array`
        '''

        V = 0

        return V

    def dVXdx(self, _V1, h2, x):
        return sp.diff(V,x)
        

    def d2VXdx2(self, h1, h2, x):

        return None

