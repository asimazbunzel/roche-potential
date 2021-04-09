'''Collection of utilities
'''

import numpy as np



class Constants(object):
    '''Object with useful constants
    '''

    def __init__(self):

        self.ln10 = 2.3025850929940455e+00
        self.pi = 3.1415926535897932384626433832795028841971693993751e0
        self.standard_cgrav = 6.67428e-8 # gravitational constant (g^-1 cm^3 s^-2)
        self.Msun = 1.9892e33  # solar mass (g)  <<< gravitational mass, not baryonic
        self.Rsun = 6.9598e10 # solar radius (cm)
        self.secyer = 3.1558149984e7 # seconds per year
        self.dayyer = 365.25e0 # days per year
        self.au = 1.495978921e13 # astronomical unit (cm)


def P_to_a(period: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Binary separation from a known period

    Parameters
    ----------
    period : `float/array`
       Binary period in days

    m1 : `float/array`
       Mass of primary star in Msun

    m2 : `flota/array`
       Mass of secondary star in Msun

    Returns
    -------
    a : `float/array`
       Binary separation in Rsun
    '''

    c = Constants()

    period = period * 24e0 * 3600e0  # in sec
    m1 = m1 * c.Msun; m2 = m2 * c.Msun  # in g

    separation = np.power(c.standard_cgrav * (m1+m2) * np.square(period/(2*pi)), one_third)

    return separation / Rsun


def a_to_P(separation: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2:Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Orbital period from a known separation

    Parameters
    ----------
    a : `float/array`
       Binary separation in Rsun

    m1: `float/array`
       Mass of primary star in Msun

    m2: `float/array`
       Mass of secondary star in Msun

    Returns
    -------
    P : `float/array`
       Binary period in days
    '''
    
    c = Constants()

    separation = separation * c.Rsun  # in cm
    m1 = m1 * c.Msun; m2 = m2 * c.Msun   # in g

    period = np.power(separation*separation*separation / (c.standard_cgrav * (m1+m2)), 0.5e0)
    period = (2*c.pi) * period

    return period / (24e0 * 3600e0)
