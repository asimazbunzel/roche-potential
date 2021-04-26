
from typing import Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
try:
    import bpy
    import bmesh
except ModuleNotFoundError:
    print('not going to import blender modules')

from roche import Roche
from utils import Constants

__all__ = ['get_xy_coords']



def get_xy_coords(R: Roche, lagrangian_point: str='L1') -> Tuple[np.ndarray, np.ndarray]:
    '''Extract (X,Y) coordinates of an equipotential of the Roche model
    
    Parameters
    ----------
    R : `roche.Roche`
       Class with the Roche model.

    lagrangian_point : `string`
       Id of the equipotential to extrat (X,Y) values.

    Returns
    -------
    X : `np.ndarray`
       X-axis coordinates of the Lagrangian point.

    Y : `np.ndarray`
       Y-axis coordinates of the Lagrangian point.
    '''

    if R is None: raise ValueError('`R` must be a valid Roche object')

    c = Constants()

    N = 1000
    x = np.linspace(-3*R.a, 3*R.a, N)
    y = np.linspace(-3*R.a, 3*R.a, N)
    X, Y = np.meshgrid(x, y)
    V = R.V(X, Y)

    if lagrangian_point == 'L1':
        xL, VL = R.L1()
    elif lagrangian_point == 'L2':
        xL, VL = R.L2()
    elif lagrangian_point == 'L3':
        xL, VL = R.L3()
    else:
        raise ValueError('`lagrangian_point` not valid. Options are: `L1`, `L2` or `L3`')

    # get contour and extract points
    cs = plt.contour(X/c.Rsun, Y/c.Rsun, -V, [-VL])
    v = cs.collections[0].get_paths()[1].vertices
    xv = v[:,0]
    yv = v[:,1]

    plt.close()

    return xv, yv
