"""
Tools to calculate the power/markers that reach the wall and work with/export the data

Javier Hidalgo Salaverri (jhsalaverri@us.es)
Lina Velarde (lvelarde@us.es)
"""

import numpy as np
from a5py import Ascot


def wall2vtp(file_in, qid_wall, file_out, get_markers: bool =  False,
            qid_run: str = 'active'):
    """
    Save a file in vtp format to allow visualization

    @params file_in: path to the .h5 file
    @params qid_wall. QID of the wall to convert. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    @params file_out: path to the .vtp file to be saved
    @params get_markers: include numbers of markers that hit the wall_id
    @params run_id: run for the markers if get_markers is True. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    """
    r = Ascot(file_in)
    try:
        int(qid_wall) # If it is a number, will not raise an error
        qid_wall = 'q'+qid_wall
    except ValueError:
        pass
    wall = r.wall.__getattribute__(qid_wall)
    vtk = wall.toVtk()

    if get_markers:
        try:
            int(qid_run) # If it is a number, will not raise an error
            qid_run = 'q'+qid_run
        except ValueError:
            pass
        data = r.__getattribute__(qid_run).endstate.read()
        markers_in_triang = np.zeros((int(wall.read()['nelements']),))
        unique, unique_counts = np.unique(data['walltile'],
            return_counts = True)
        for i in range(len(markers_in_triang)):
            if i+1 in unique:
                markers_in_triang[i] = unique_counts[unique == i+1]
            else:
                markers_in_triang[i] = 0
        vtk.addScalarField('nmarkers', markers_in_triang)

    vtk.writeVtp(file_out)


def write_vtkwloads(fn, qid = 'active', fn_vtk = 'test_wallloads.vtk', coord = 'cyl',
                    w_indices=None, p_ids=None, flag_return = False):
    """
    Write a vtk file containing:
            # - "pload" particle load in units of prt/m^2 or prt/m^2s,
            # - "eload" power/energy load in units of W/m^2 or J/m^2
            # - "mload" marker load in units of markers
            # - "iangle" angle of incidence (the angle between power flux and the surface normal) in deg
            # - "area" area of the wetted triangles
            # - "barycenter" barycenter of the wetted triangles

    Inputs:
    fn : str
        Path to the ascot file.
    qid : str, optional 
        qid of the run. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    fn_vtk : str, optional 
        Path to the file that will be written.
    coord : str, optional
        Coordinates for the barycenter. Default is cylindrical ('cyl'), but can also be cartesians ('car')
    w_indices : array_like, optional
        List of triangle indecies for which the 3D mesh is made.
    p_ids : array_like, optional
        List of particle ids for which the wall loads are calculated.
    flag_return: bool, optional.
        If True, the mesh data will be returned.


    Note on the units:
    Ascot is usually run assuming steady-state, meaning that you'll get 
    a constant flux of losses (W/m^2). The code can be used also to 
    simulate existing  particle population, not source as in steady-state, 
    but then the losses are in units of J/m^2.
    The steady-state assumption is currently hard-coded (e.g. weight is in 
    particles/s and not ever just particles) so it's user's responsibility 
    to keep track of the units until the non-steady state support is fully developed
    """

    print('~~~~~~~~~~~~~~~~~~')
    print(fn)
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)
    
    polywall = sim.getwall_3dmesh(w_indices=w_indices, p_ids=p_ids)
    ids, _, _, _, _ = sim.getwall_loads(p_ids=p_ids)
    ids = ids - 1 # Convert IDs to indices
    if coord == 'cyl':
        barycenters = sim.wall.barycenters(cartesian=False)
    elif coord == 'car':
        barycenters = sim.wall.barycenters(cartesian=True)
    else:
        print('Coordinates not recognised. Exiting...')
        return

    n_tri = sim.wall.read()['nelements']

    wall_f = np.arange(n_tri, dtype=int)
    if w_indices is not None:
        wall_f = wall_f[w_indices]

    if coord == 'cyl':
        baryc_r = np.zeros((n_tri, )) + np.nan
        baryc_r[ids] = barycenters[ids][:,0]
        polywall.cell_data["barycenter(r)"]  = baryc_r[wall_f]
        baryc_phi = np.zeros((n_tri, )) + np.nan
        baryc_phi[ids] = barycenters[ids][:,1]
        polywall.cell_data["barycenter(phi)"]  = baryc_phi[wall_f]
    elif coord == 'car':
        baryc_x = np.zeros((n_tri, )) + np.nan
        baryc_x[ids] = barycenters[ids][:,0]
        polywall.cell_data["barycenter(x)"]  = baryc_x[wall_f]
        baryc_y = np.zeros((n_tri, )) + np.nan
        baryc_y[ids] = barycenters[ids][:,1]
        polywall.cell_data["barycenter(y)"]  = baryc_y[wall_f]
    baryc_z = np.zeros((n_tri, )) + np.nan
    baryc_z[ids] = barycenters[ids][:,2]
    polywall.cell_data["barycenter(z)"]  = baryc_z[wall_f]

    polywall.save(fn_vtk)

    if flag_return:
        return polywall


