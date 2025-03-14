"""
Save a certain orbit. NOT CHECKED.

by Javier Hidalgo Salaverri (jhsalaverri@us.es)
"""

import numpy as np
from a5py import Ascot


def save_orbit(file, orbit_id, file_out, run_qid: str = 'active'):
    """
    Save an orbit to a .txt file

    :param file: path to the .h5 file
    :orbit id: 
    :run_qid: run qid
    @TODO: add a header
    """

    # Get the orbits
    a5 = Ascot(file)
    orbits = a5.__getattribute__(run_qid).orbit.read()
    index = orbits['ids'] == orbit_id
    r = orbits['r'][index]
    z = orbits['z'][index]
    r = r
    z = z

    np.savetxt(file_out, [r,z], delimiter = ';')
