"""
Routines to play with the wall
@Author: Lina Velarde - <lvelarde@us.es>
"""
import pyvista as pv
from a5py import Ascot
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from a5py.ascot5io import wall_2D as wall_2D


def wall_to_stl(fn, qidwall, fnout):
    """
    Convert a 3D wall to .stl File

    Inputs:
    fn. Ascot filename
    qidwall. qid of the wall to convert. Can only be the 10 digit ID number.
    fnout. Filename of the stl file. Must end with .stl
    """

    h5 = Ascot(fn)
    dum = h5.data.wall['q'+qidwall].tomesh()
    pv.save_meshio(fnout, dum)


def move_wall(path, rext, device = 'MU'):
    """
    Script to move the wall close to separatrix. NOT ADAPTED TO NEW ASCOT.
    Also, only works for MAST-U at the moment.

    Inputs:
        path. Path to ascot file. Including ascot.h5.
        rext. R at which we want the outer wall. Smth like 1.485.
        device. MU or any other
    """
    h5 = Ascot(path)
    if device == 'MU':
        n = h5.wall.q1734340493.read()['nelements']
        r = h5.wall.q1734340493.read()['r']
        z = h5.wall.q1734340493.read()['z']
        indr = np.where(r == 1.521)
        rnew = r
        rnew[indr] = rext
    else:
        print('Other devices not implemented yet')
    # Write new input
    wall2D.write_hdf5(path, n, rnew, z, desc='wall2D_v1_1p'+str(rext)[2:])