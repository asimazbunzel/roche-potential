'''Collection of classes to study Roche geometry
'''

from typing import Union

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

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

        c = Constants()

        self.m1 = m1 * c.Msun
        self.m2 = m2 * c.Msun
        self.a = a * c.Rsun
        self.P = P * 24e0 * 3600e0

        self.N = N

        self.M = self.m1 + self.m2

        if self.m1 > self.m2:
            self.q = self.m2 / self.m1
        else:
            self.q = self.m1 / self.m2


    def _Vdl(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray],
            z: Union[float, np.ndarray]=None) -> Union[float, np.ndarray]:
        '''Dimensionless Roche potential
        
        Provided a grid of points, returns the value of the potential in the mesh
        
        Parameters
        ----------
        x : `float/np.ndarray`

        y : `float/np.ndarray`

        z : `float/np.ndarray`

        Returns
        -------
        V : `float/np.ndarray`
        '''

        if z is None:
            r1 = np.sqrt(x**2 + y**2)
            r2 = np.sqrt((x-1)**2 + y**2)
        else:
            r1 = np.sqrt(x**2 + y**2 + z**2)
            r2 = np.sqrt((x-1)**2 + y**2 + z**2)

        V = 2 / ((1+self.q) * r1) + 2 * self.q / ((1+self.q) * r2) + (x - self.q / (1+self.q))**2
        V += y**2

        return V


    def V(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray],
            z: Union[float, np.ndarray]=None) -> Union[float, np.ndarray]:
        '''Dimensionless Roche potential
        
        Provided a grid of points, returns the value of the potential in the mesh
        
        Parameters
        ----------
        x : `float/np.ndarray`

        y : `float/np.ndarray`

        z : `float/np.ndarray`

        Returns
        -------
        V : `float/np.ndarray`
        '''
        
        c = Constants()

        x = x / self.a
        y = y / self.a
        if z is not None: z = z / self.a

        V = - c.standard_cgrav * self.M / (2 * self.a) * self._Vdl(x, y, z)
        
        return V


    def _dVdldx(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray]=0,
            z: Union[float, np.ndarray]=0) -> Union[float, np.ndarray]:

        if z is None:
            r1 = np.sqrt(x**2 + y**2)
            r2 = np.sqrt((x-1)**2 + y**2)
        else:
            r1 = np.sqrt(x**2 + y**2 + z**2)
            r2 = np.sqrt((x-1)**2 + y**2 + z**2)

        dV = 2 / (1+self.q) * (-x / r1**3) + 2 * self.q / (1+self.q) * (-(x-1) / r2**3)
        dV += 2 * (x - self.q / (1+self.q))

        return dV


    def L1(self):

        eps = 1e-5

        a = eps  # avoid being in m1 where potential diverges
        b = 1-eps  # also avoid m2 position

        xL1 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL1 = self.V(xL1, 0, 0)

        return xL1, VL1
    
    
    def L2(self):

        eps = 1e-5

        a = 1 + eps  # just move a little bit from m2 position
        b = 30  # go as far as 30 times the position of m2

        xL2 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL2 = self.V(xL2, 0, 0)

        return xL2, VL2
    
    
    def L3(self):

        eps = 1e-5

        a = -30  # go 30 times left of m1
        b = -eps  # also avoid m1 position

        xL3 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL3 = self.V(xL3, 0, 0)

        return xL3, VL3
