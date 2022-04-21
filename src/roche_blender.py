
from typing import Tuple, Union
import sys
sys.path.append('/home/asimazbunzel/Projects/mass-transfer/roche-potential/src/')

import numpy as np
#import matplotlib.pyplot as plt
try:
    import bpy
    import bmesh
except ModuleNotFoundError:
    print('not going to import blender modules')
import sys
import os

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )
    #print(sys.path)
from roche import Roche
from utils import Constants

__all__ = ['get_xy_coords']

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()

def get_xy_coords(R, lagrangian_point: str='L1') -> Tuple[np.ndarray, np.ndarray]:
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

    if not isinstance(R,Roche): raise ValueError('`R` must be a valid Roche object')

    if lagrangian_point != 'L1' and lagrangian_point != 'L2' and lagrangian_point != 'L3':
        raise ValueError('`lagrangian_point` not valid. Options are: `L1`, `L2` or `L3`')

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
#    cs = plt.contour(X/c.Rsun, Y/c.Rsun, -V, [-VL])
    v = cs.collections[0].get_paths()[1].vertices
    xv = v[:,0]
    yv = v[:,1]
    
    mask = yv > 0
    
    xv = xv[mask]
    yv = yv[mask]

#    plt.close()
    
    return xv, yv


R = Roche(m1=10, m2=5, a=5)

x, y = get_xy_coords(R, 'L1')

vertices = []
for k in range(len(x)):
    vertices.append([x[k], y[k], 0e0])

# make mesh
# vertices = [(0, 0, 0), (1, 0, 0)]
edges = []
faces = []

new_mesh = bpy.data.meshes.new('L1_mesh')
new_mesh.from_pydata(vertices, edges, faces)
new_mesh.update()

# make object from mesh
new_object = bpy.data.objects.new('L1_obj', new_mesh)

# make collection
new_collection = bpy.data.collections.new('L1_coll')
bpy.context.scene.collection.children.link(new_collection)

# add object to scene collection
new_collection.objects.link(new_object)
