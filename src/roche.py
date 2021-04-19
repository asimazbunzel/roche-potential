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

    .. math::

       \\phi = - \\dfrac{G M_1}{r_1} - G \\dfrac{G M_2}{r_2} - \\dfrac{\\omega^2 r^{2}_3}{2}

    where G is the gravitational constant, r_1 and r_2 the distances to the center of the stars
    having masses M_1 and M_2, respectively; \omega is the orbital angular velocity which can be
    derived from Kepler's law:

    .. math::
    
       \\omega = \\dfrac{2 \\pi}{P} = \\sqrt{\\dfrac{G (M_1 + M_2)}{a^3}}

    and r_3 is the distance to the axis of rotation.

    Instead, we will use a reference frame centered on the more massive star (M_1) such that the
    potential is:

    .. math::

       \\phi = - \\dfrac{G M_1}{r_1} - \\dfrac{G M_2}{r_2} - \\dfrac{\\omega^2}{2} \\left[ \\left( x - \\dfrac{M_2}{M_1 + M_2} \\right)^2 + y^2 \\right]

    where r_1 = \\sqrt{x^2 + y^2 + z^2} and r_2 = \\sqrt{(x-1)^2 + y^2 + z^2}

    Parameters
    ----------
    m1 : `float`
       Mass of the primary (more massive star) in Msun.

    m2 : `float`
       Mass of the secondary (less massive star) in Msun.

    a : `float`
       Orbital separation (in Rsun).
    '''

    def __init__(self, m1: float, m2: float, a: float=None, P: float=None) -> None:

        if P is None and a is None: raise ValueError('either `P` or `a` must be specified')

        if a is None: a = P_to_a(P, m1, m2)
        if P is None: P = a_to_P(a, m1, m2)

        c = Constants()

        self.m1 = m1 * c.Msun
        self.m2 = m2 * c.Msun
        self.a = a * c.Rsun
        self.P = P * 24e0 * 3600e0

        self.M = self.m1 + self.m2

        if self.m1 > self.m2:
            self.q = self.m2 / self.m1
        else:
            self.q = self.m1 / self.m2


    def _Vdl(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray],
            z: Union[float, np.ndarray]=None) -> Union[float, np.ndarray]:
        '''Dimensionless Roche potential
        
        Provided a grid of points, returns the value of the dimensionless potential in the mesh
        
        Parameters
        ----------
        x : `float/np.ndarray`
           Coordinates for the x-axis.

        y : `float/np.ndarray`
           Coordinates for the y-axis.

        z : `float/np.ndarray`
           Coordinates for the z-axis.

        Returns
        -------
        V : `float/np.ndarray`
           Dimensionless Roche potential on the (x, y, z) coordinates.
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
        '''Roche potential
        
        Provided a grid of points, returns the value of the potential in the mesh
        
        Parameters
        ----------
        x : `float/np.ndarray`
           Coordinates for the x-axis.

        y : `float/np.ndarray`
           Coordinates for the y-axis.

        z : `float/np.ndarray`
           Coordinates for the z-axis.

        Returns
        -------
        V : `float/np.ndarray`
           Roche potential on the (x, y, z) coordinates.
        '''
        
        c = Constants()

        x = x / self.a
        y = y / self.a
        if z is not None: z = z / self.a

        V = - c.standard_cgrav * self.M / (2 * self.a) * self._Vdl(x, y, z)
        
        return V


    def _dVdldx(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray]=0,
            z: Union[float, np.ndarray]=0) -> Union[float, np.ndarray]:
        '''Derivative of the dimensionless Roche potential along the x-axis

        Parameters
        ----------
        x : `float/np.ndarray`
           Coordinates for the x-axis.

        y : `float/np.ndarray`
           Coordinates for the y-axis.

        z : `float/np.ndarray`
           Coordinates for the z-axis.

        Returns
        -------
        dV : `float/np.ndarray`
           Derivative of dimensionless Roche potential on the (x, y, z) coordinates.
        '''

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
        '''Find L1 Lagrangian point
        
        Returns
        -------
        xL1 : `float`
           x coordiante of L1 in cm.
        VL1 : `float`
           Value of Roche potential on L1 (cgs units).
        '''

        eps = 1e-5

        a = eps  # avoid being in m1 where potential diverges
        b = 1-eps  # also avoid m2 position

        xL1 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL1 = self.V(xL1, 0, 0)

        return xL1, VL1
    
    
    def L2(self):
        '''Find L2 Lagrangian point
        
        Returns
        -------
        xL2 : `float`
           x coordiante of L2 in cm.
        VL2 : `float`
           Value of Roche potential on L2 (cgs units).
        '''

        eps = 1e-5

        a = 1 + eps  # just move a little bit from m2 position
        b = 30  # go as far as 30 times the position of m2

        xL2 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL2 = self.V(xL2, 0, 0)

        return xL2, VL2
    
    
    def L3(self):
        '''Find L3 Lagrangian point
        
        Returns
        -------
        xL3 : `float`
           x coordiante of L3 in cm.
        VL3 : `float`
           Value of Roche potential on L3 (cgs units).
        '''

        eps = 1e-5

        a = -30  # go 30 times left of m1
        b = -eps  # also avoid m1 position

        xL3 = optimize.brentq(self._dVdldx, a, b) * self.a
        VL3 = self.V(xL3, 0, 0)

        return xL3, VL3


    def eq_vol(self, k, eq='L1', N=100000):
        '''Compute volume of an equipotential around a certain star in the binary
        
        For the three Lagrangian points (L1, L2, L3), one can compute the value of the
        dimensionless Roche potential. If a plane is assumed to cross the L1 point, one can compute
        the volume of each equipotential for the Lagrangian points which we do with a MonteCarlo
        approach. Additionally, one can compute the associated radii for an sphere of the same
        volume, similar to the approach used by Eggleton (1983)

        Parameters
        ----------
        k : `integer`
           Id of the star around which the volume will be computed. Either `1` or `2` for the more
           (less) massive, respectively.

        eq : `string`
           Id of the equipotential to use. Options are: `L1`, `L2` or `L3`.

        N : `integer`
           Number of random draws for the MonteCarlo computation of the volume.

        Returns
        -------
        Vol : `float`
           Dimensionless volume inside the equipotential chosen around one of the stars.

        Reff : `float`
           Associated effective radii for the `vol` computed in units of separation.
        '''

        if k != 1 and k != 2:
            raise ValueError('`k` must be either 1 (m1) or 2 (m2).')

        if eq != 'L1' and eq != 'L2' and eq != 'L3':
            raise ValueError('`eq` must be either `L1`, `L2` or `L3`')

        # Lagrangian point positions and potentials in dimensionless units
        xL1, _ = self.L1()
        xL2, _ = self.L2()
        xL3, _ = self.L3()
        xL1 = xL1 / self.a
        xL2 = xL2 / self.a
        xL3 = xL3 / self.a
        dVl1 = self._Vdl(xL1, 0, 0)
        dVl2 = self._Vdl(xL2, 0, 0)
        dVl3 = self._Vdl(xL3, 0, 0)

        if eq == 'L1': eq = dVl1
        if eq == 'L2': eq = dVl2
        if eq == 'L3': eq = dVl3

        # we split space with a plane crossing the L1 point
        f = lambda x: self._Vdl(x, 0, 0) - dVl1

        # depending on the id of the star, we will get the limits of the lobes which will
        # give the radius of a sphere around the id of the star
        eps = 1e-5
        if k == 1:
            x1 = optimize.brentq(f, xL3, eps)
            x2 = xL1
            r = max(abs(x1), abs(x2))
        else:
            x1 = xL1
            x2 = optimize.brentq(f, 1+eps, xL2)
            r = max(1-x1, x2-1)

        # get random numbers with the right value in x coordinate. y and z will be used
        # to count which cases are inside the lobe
        x = r * (2 * np.random.random(N) - 1) + (k-1)
        y = r * (2 * np.random.random(N) - 1)
        z = r * (2 * np.random.random(N) - 1)

        # count points inside lobe which is the effective volume of the lobe
        ps = self._Vdl(x, y, z)
        kk = np.where(ps > eq)[0]
        ef = len(kk) / N
        Vol = (2 * r)**3 * len(kk) / N

        # compute the radius of a sphere with the same volume
        Reff = (Vol * (3 / 4) / np.pi)**(1/3)

        return Vol, Reff
