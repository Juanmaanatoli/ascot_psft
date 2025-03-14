"""
Auxiliary utilities to work with ASCOT

NOT CHECKED

by Javier Hidalgo Salaverri (jhsalaverri@us.es)
"""

import numpy as np
import matplotlib.pyplot as plot
from a5py import Ascot
import subprocess
from datetime import datetime
import h5py
import matplotlib.path as mpltPath

def read_old_options(file, run_id = None):
    """
    Read the options from an old ASCOT run that's not readable by the up-to-date routines

    :param file: path to the h5 file
    :param run_id: run id. If None get active one
    """
    if run_id is None:
        options = Ascot(file).data.options.active
    else:
        options = Ascot(file).data.options.__getattribute__(run_id)

    fn   = options._root._ascot.file_getpath()
    path = options._path

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            val = f[path][key][:]
            out[key] = val
    return out

# def compare_options(file1, qid1 = None, file2, qid2=None):
#     f1 = Ascot(file1)
#     f2 = Ascot(file2)
#     if qid1 is None:
#         opt1 = f1.data.options.active
#     else:
#         opt1 = f1.data.options.__getattribute__(qid1)
#     if qid2 is None:
#         opt2 = f2.data.options.active
#     else:
#         opt2 = f2.data.options.__getattribute__(qid2)

    


def time_in_seconds():
    "Returns current time in seconds. Ref 01/01/1970"
    time = datetime.now()
    time_ref = datetime(1970, 1, 1)
    return (time-time_ref).total_seconds()

def date2sec(date:str):
    "Translate a h5 date to seconds to allow comparisons. Ref 01/01/1970"
    year = int(date[0:4])
    month = int(date[5:7])
    day = int(date[8:10])
    hour = int(date[11:13])
    min = int(date[14:16])
    sec = int(date[17:19])
    time = datetime(year, month, day, hour, min, sec)
    time_ref = datetime(1970, 1, 1)
    return (time-time_ref).total_seconds()


def wait4sim(user: str = 'jhidalg1', check_time: float = 1.0,
    return_time: bool = False):
    """
    Waits until the queue for a certain user is clean.

    :param user: user name
    :param check_time: how often it is checked (in seconds)
    :param return_time: return waited time

    :return time waited
    """
    t0 = time_in_seconds()

    cmd = 'squeue -u '+str(user)
    flag = True
    t1 = time_in_seconds()
    while flag == True:
        if time_in_seconds() > t1+check_time:
            check_queue = subprocess.check_output(cmd, shell = True)
            check_queue = check_queue.decode('utf-8')
            if check_queue.find(user) < 0:
                flag = False
            t1 = time_in_seconds()
    return time_in_seconds()-t0


def get_last_run_id(file, attribute: str = ''):
    """
    Return the id of the last added element.

    :params file_id
    :params attribute. If empty, return the last ASCOT run, otherwise, the given
        attribute (i.e., marker)
    """
    a5 = Ascot(file)
    if attribute != '':
        a5 = a5.__getattribute__(attribute)

    times = [date2sec(a5.__getattribute__(q).get_date()) for q in a5._qids]
    index = np.argmax(times)
    return a5._qids [index]

def copy_inputs(file_in, file_out):
    """
    Create a copy of a h5 file without the run data

    """
    inputs = ['options','bfield', 'efield', 'marker', 'plasma', 'neutral',
        'wall', 'boozer', 'mhd', 'nbi']
    a5 = Ascot(file_in)

    for inp in inputs:
        for qid in a5.__getattribute__(inp)._qids:
            a5.__getattribute__(inp).__getattribute__(qid).copy_to_hdf5file(file_out)


def pol2cart(rho, phi):
    """
    Convert from poloidal to cartesian coords.
    Phi should be in radians here.
    """
    x = rho*np.cos(phi)
    y = rho*np.sin(phi)
    return(x, y)

def check_point_inside_vessel(points, vessel = 'AUG'):
    """
    Check if a point is inside a vessel

    :param points: array with the points to check
    :param vessel: vessel to check. Options: AUG
    """
    # TODO: Dont hardcode the vessel location
    vessel_file = '/marconi_work/FUA38_ASCOT_3D/jhidalg1/ascot_psft/preprocessing/vessel_mesh/vessel_mesh.h5'
    v = Ascot(vessel_file)
    if vessel == 'AUG': 
        wall = v.data.wall.AUG2D.read()
    else:
        raise ValueError('Vessel not recognized')

    rw = wall['r'].ravel()
    zw = wall['z'].ravel()
    w = np.array([rw, zw]).T
    path = mpltPath.Path(w)
    inside = path.contains_points(points)

    return inside