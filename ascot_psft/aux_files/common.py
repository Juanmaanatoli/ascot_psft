"""
Script with some common routines that will be used in other scripts.

@author: L. Velarde - lvelarde@us.es


List of functions:
    read_bf. To return the most usual bf variables
"""
import numpy as np
from a5py import Ascot
import matplotlib.pyplot as plt
from a5py.ascot5io import nbi as a5nbi
from a5py.ascot5io import wall as a5wall
from a5py.ascot5io import bfield as a5bf
from a5py.ascot5io import marker as a5marker
from a5py.ascot5io import plasma as a5plasma
from a5py.ascot5io import options as a5options


def copy_input(fn_from: str, fn_to: str, type_inp: str, qid_from: str, desc_to = None):
    """
    Copy input from one ascot file to another.
    
    Inputs:
        fn_from. Name of the file to copy from.
        fn_to. Name of the file to copy to.
        type_inp. Type of input to copy. Can be:
            - 'options'
            - 'bfield'
            - 'marker'
            - 'plasma'
            - 'wall'
            - 'nbi'
        qid_from. QID of the input to copy.
        desc_to. Description of the input to copy to. If None, will keep the same description.
    """
    h5 = Ascot(fn_from)
    
    if type_inp == 'options':
        opt = h5.data.options['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.options['q'+qid_from].get_desc()
        a5options.Opt.write_hdf5(fn = foutput, desc = desc_to, **opt)
        
    elif type_inp == 'bfield':
        bf = h5.data.bfield['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.bfield['q'+qid_from].get_desc()
        formatinp = h5.data.bfield['q'+qid_from].get_type()
        if formatinp == 'B_2DS':
            a5bf.B_2DS.write_hdf5(fn = fn_to, desc = desc_to, **bf)
        elif formatinp == 'B_3DS':
            a5bf.B_3DS.write_hdf5(fn = fn_to, desc = desc_to, **bf)
        else:
            print('Format not implemented yet')
            return
            
    elif type_inp == 'marker':
        print('Reading original marker input...')
        marker = h5.data.marker['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.marker['q'+qid_from].get_desc()
        formatinp = h5.data.marker['q'+qid_from].get_type()
        if formatinp == 'prt':
            print('Writing new prt input...')
            a5marker.Prt.write_hdf5(fn = fn_to, desc = desc_to, **marker)
        else:
            print('Format not implemented yet')
            return    
        
    elif type_inp == 'plasma':
        plasma = h5.data.plasma['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.plasma['q'+qid_from].get_desc()
        formatinp = h5.data.plasma['q'+qid_from].get_type()
        if formatinp == 'plasma_1D':
            a5plasma.plasma_1D.write_hdf5(fn = fn_to, desc = desc_to, **plasma)
        else:
            print('Format not implemented yet')
            return
        
    elif type_inp == 'wall':
        wall = h5.data.wall['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.wall['q'+qid_from].get_desc()
        formatinp = h5.data.wall['q'+qid_from].get_type()
        if formatinp == 'wall_2D':
            a5wall.wall_2D.write_hdf5(fn = fn_to, desc = desc_to, **wall)
        elif formatinp == 'wall_3D':
            a5wall.wall_3D.write_hdf5(fn = fn_to, desc = desc_to, **wall)
        else:
            print('Format not implemented yet')
            return
    
    elif type_inp == 'nbi':
        nbi = h5.data.nbi['q'+qid_from].read()
        if desc_to is None:
            desc_to = h5.data.nbi['q'+qid_from].get_desc()
        formatinp = h5.data.nbi['q'+qid_from].get_type()
        if formatinp == 'nbi':
            a5nbi.NBI.write_hdf5(fn = fn_to, desc = desc_to, **nbi)
        else:
            print('Format not implemented yet')
            return
        
    # Give some info about the input copied: format, qid, description
    print('Input copied from {} to {}'.format(fn_from, fn_to))
    print('Input type: {}'.format(type_inp))
    print('Format: {}'.format(formatinp))
    print('QID: {}'.format(qid_from))
    print('Description: {}'.format(desc_to))
    

def plot_vessel_tor(device = 'MU', cwall = 'k', ax = None):
    """
    Plot the toroidal projection of the vessel. Mostly used as an aux for other plotting routines.
    @author: Lina Velarde <lvelarde@us.es>


    Inputs:
        device. Device to plot. Can be:
            - 'AUG'
            - 'MU'
        ax. axis where to plot, if None a new one wil be created
    """
    if ax is None:
        fig, ax = plt.subplots()

    if device == 'MU':
        win = 0.25 # position of the inner wall
        wout = 1.56 # position of the outter wall
        sepin = 0.33 # position of the inner separatrix
        sepout = 1.4 # position of the outter separatrix
        range_xy = [[-1.65,1.65],[-1.65,1.65]]
        xlim = [-1.65,1.65]
        ylim = [-1.65,1.65]
    
    elif device == 'AUG':
        print('Not fully implemented. Please correct position of the wall and separatrix.')
        win = 0.1 # position of the inner wall
        wout = 2 # position of the outter wall
        sepin = 0.33 # position of the inner separatrix
        sepout = 1.4 # position of the outter separatrix
        range_xy = [[-2,2],[-2,2]]
        xlim = range_xy[0]
        ylim = range_xy[1]
        
    
    circle1 = plt.Circle((0, 0), win, color=cwall, fill=False)
    circle2 = plt.Circle((0, 0), wout, color=cwall, fill=False)
    circle3 = plt.Circle((0, 0), sepin, color='r', fill=False)
    circle4 = plt.Circle((0, 0), sepout, color='r', fill=False)
    
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.add_patch(circle4)

    return xlim, ylim

def read_bf(bf):
    """
    To return the most usual bf variables

    Inputs:
        bf. Magnetic field dictionary read from the ascot input
    """
    psi = (bf["psi"])
    nr = len(bf["psi"]); nz = len(bf["psi"][0])
    try:
        rmin = bf["b_rmin"]; rmax = bf["b_rmax"]
        zmin = bf["b_zmin"]; zmax = bf["b_zmax"]
        phimin = bf["b_phimin"][0]; phimax = bf["b_phimax"][0]
        nphi = bf['b_nphi'][0]; 
        phivec = np.linspace(phimin,phimax,nphi)
        bf2D = False
    except KeyError:
        rmin = bf["rmin"]; rmax = bf["rmax"]
        zmin = bf["zmin"]; zmax = bf["zmax"]
        nphi = None; phivec = None
        bf2D = True
    zvec = np.linspace(zmin,zmax,nz); rvec = np.linspace(rmin,rmax,nr)
    zvec = zvec.transpose()[0]; rvec = rvec.transpose()[0]
    psi1 = bf["psi1"]; psi0 = bf["psi0"]

    return nr, rvec, nz, zvec, nphi, phivec, psi, psi1, psi0, bf2D


def eV2J(quant):
    """
    Convert from ev to Joules

    Inputs:
        quant. Quantity in eV to be converted.
    """
    quant_J = quant * 1,60218e-19

    return quant_J
