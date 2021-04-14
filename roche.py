'''Collection of classes to study Roche geometry
'''

import numpy as np
import matplotlib.pyplot as plt

from utils import Constants, a_to_P, P_to_a



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

    def __init__(self, m1: float, m2: float, a: float=None, P: float=None, N: int=1000) -> None:

        if P is None and a is None: raise ValueError('either `P` or `a` must be specified')

        if a is None: a = P_to_a(P, m1, m2)
        if P is None: P = a_to_P(a, m1, m2)

        self.m1 = m1 * c.Msun
        self.m2 = m2 * c.Msun
        self.a = a * c.Rsun
        self.P = P * 24e0 * 3600e0

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

    def Vf(self, x: np.ndarray, y: np.ndarray, m1: float, m2: float, r_m1: np.array, r_m2: np.array) -> np.ndarray:
        '''Roche potential
        
        Provided a grid of points, returns the value of the potential in the mesh.
        
        Parameters
        ----------
        x : `np.ndarray` x coordinate of field point

        y : `np.ndarray` y coordinates of field points

        m1 : `float` mass #1

        m2 : `float` mass #2

        r_m1 : `np.array` position of masses
        
        r_m2 : `np.array` position of masses

        Returns
        -------
        V : `np.ndarray` 
        '''

        x_m1 = r_m1[0]*np.ones_like(x)           #X coordinate of mass 1
        y_m1 = r_m1[1]*np.ones_like(y)           #Y coordinate of mass 1
        x_m2 = r_m2[0]*np.ones_like(x)           #X coordinate of mass 2
        y_m2 = r_m2[1]*np.ones_like(y)           #Y coordinate of mass 2
        x_mc=(x_m1*m1+x_m2*m2)/(m1+m2)           #X coordinate of mass center 
        y_mc=(y_m1*m1+y_m2*m2)/(m1+m2)           #Y coordinate of mass center

        r_1 = np.sqrt((x-x_m1)**2 + (y-y_m1)**2) #Mesh of |r-r1|
        r_2 = np.sqrt((x-x_m2)**2 + (y-y_m2)**2) #Mesh of |r-r2|
        
        a = np.linalg.norm((r_m2-r_m1),2)        #Distance from r2 to r1 = |r2 - r1|
        
        V = -c.G.value * ( m1/r_1 + m2/r_2 + (m1+m2)/a**3 * ((x-x_mc)**2 + (y-y_mc)**2))

        return V

    def dVXdx(self, h1, h2, x):

        return None

    def d2VXdx2(self, h1, h2, x):

        return None

