"""
Generate the SMART neutral beam injector geometry for BBNBI5.

J.F. Rivero-Rodriguez <jfrivero@us.es> || L. Velarde <linvelgal@alum.us.es>
"""
import numpy as np
from scipy.constants import physical_constants as const

import a5py.ascot5io.nbi as nbi
import a5py.marker.interpret as a5interpret


def generate(fn, desc, Rt = 0.44, inj = [1], specie = "H",
            E = 25e3, P = 1):
    """
    Generate SMART injector.

    Args:
        fn : str
            Full path to the HDF5 file.
        Rt : float
            Tangency radius. This is the parameter that we will scan.
        inj : list of int, default = all
            List of SMART injectors to be generated (there's only one,
            but there are multiple configurations).
            Config: Radialdist = 2.0, thetasteer = 0, phisteer = 14.036º
        desc : str, optional
            Input description.
        specie : {"D","H"}
            Specie injected by the NBI. Deuterium "D" (default) or Protons "H".
    """

    ## Convert injector numbers into index
    inj = np.array(inj)-1

    ## Define injectors geometry
    Rtan = np.array([Rt]) #defining the tangencial radius
                            #both phi steer and Rcenter can be calculated
    Vtan = np.array([0.00])   #vertical distance where the beam hits
                              #measured from the magn axis
    RadialDist = 2 #this distance is measured from the grid to the centre
                   #of the tokamak
    Rcenter = np.sqrt(RadialDist**2-Rtan**2)

    phicenter = np.array([-180])*np.pi/180
    zcenter = np.array([0.00])
    phisteer = -np.arctan(Rtan/RadialDist) #radians
    thetasteer = np.arctan(Vtan/Rcenter) #radians

    focal_length_h = np.array([5])
    focal_length_v = focal_length_h * 6 / 14

    tilt = 0.5*np.pi/180.

    div_halo_frac = np.array([0, #unknown data
                              ])
    div_halo_h= np.array([1e-10, #unknown data
                         ])
    div_halo_v= np.array([1e-10, #unknown data
                         ])
    if specie == "D":
        anum = 2
        znum = 1
        mass = a5interpret.mass_kg("deuterium")
        energy = np.array([1.5e4])*const["elementary charge"][0]
        power = np.array([2.88])*1e6
        efrac = np.array([[0.62,0.27,0.11]])
        div_h = np.array([0.6])*np.pi/180.
        div_v = np.array([0.6])*np.pi/180.

    elif specie == "H":
        anum = 1
        znum = 1
        mass = a5interpret.mass_kg("hydrogen")
        energy = np.array([E])*const["elementary charge"][0]
        power = np.array([P])*1e6
        efrac = np.array([[0.62,0.27,0.11]])
        div_h = np.array([0.6])*np.pi/180.
        div_v = np.array([0.6])*np.pi/180.
    else:
        raise ValueError("{} is not a recognised injection specie".format(specie))

    # Make nbi input of the selected injectors
    out = []
    for i in inj:
        ## Define beamlets grid
        grid = define_grid(Rcenter[i],phicenter[i],zcenter[i],
                           thetasteer[i],phisteer[i],
                           tilt,focal_length_h[i],focal_length_v[i])

        out.append({})
        out[-1]["id"] = i+1
        out[-1]["nbeamlet"] = grid["xyz"][:,0].size
        out[-1]["beamletx"] = grid["xyz"][:,0]
        out[-1]["beamlety"] = grid["xyz"][:,1]
        out[-1]["beamletz"] = grid["xyz"][:,2]
        out[-1]["beamletdx"] = grid["dxdydz"][:,0]
        out[-1]["beamletdy"] = grid["dxdydz"][:,1]
        out[-1]["beamletdz"] = grid["dxdydz"][:,2]
        out[-1]["div_h"] = div_h[i]
        out[-1]["div_v"] = div_v[i]
        out[-1]["div_halo_frac"] = div_halo_frac[i]
        out[-1]["div_halo_h"] = div_halo_h[i]
        out[-1]["div_halo_v"] = div_halo_v[i]
        out[-1]["anum"] = anum
        out[-1]["znum"] = znum
        out[-1]["mass"] = mass
        out[-1]["energy"] = energy[i]
        out[-1]["efrac"] = efrac[i,:]
        out[-1]["power"] = power[i]

    return nbi.write_hdf5(fn, out, desc)


def define_grid(Rcenter,phicenter,zcenter,
                thetasteer,phisteer,
                tilt,fl_h,fl_v):

    """
    Define beamlets grid in local coordinates  (centre of
    grid at x=0, y=0, z=0) and then allocate and orientate
    the grid according to the given location and orientation
    of the beam.

    Args:
        Rcenter, phicenter, zcenter : floats [m,rad,m].
            Location of the grid center in cilindrical coordinates.
        thetasteer: float [rad]
            Angle around the horizontal axis by which the beam center
            line is oriented.
        phisteer: float [rad]
            Angle around the vertical axis by which the beam center
            line is oriented.
        tilt: float [rad]
            Angle around the horizontal axis by which half of the beamline
            grid is tilted. This angle makes grid to be tilted in halves,
            producing a shape similar to the one shown below:
                \
                 \
                  \
                  /
                 /
                /
        fl_h, fl_v: floats [m,m]
            Horizontal and vertical focal lengths.
    """

    grid = {}

#     # Define beamlets position in local coordinates
#
#     Divide each gridhalf to 3 blocks: 2 of which create the main diamond
#     pattern, and one that add the top row.
#     Generate a rows-by-columns matrix of beamlets where
#     the extreme points are given as (front view):
#
#     TopLeft ---------- LastCol
#      \                       \
#       \                       \
#        \                       \
#       LastRow ------------------
#
#
# MAST-U GRID
#     Rows = np.array([6,6,1])
#     Cols = np.array([10,11,5])
#
#     TopLeft = np.array([[-7.42e-2, 20.04e-2],
#                         [-8.25e-2, 18.34e-2],
#                         [-3.3e-2, 21.74e-2]])
#
#     LastRow = np.array([[-7.42e-2, 3.04e-2],
#                         [-8.25e-2, 1.34e-2],
#                         [-3.3e-2, 0]])
#
#
#     LastCol = np.array([[7.42e-2, 20.04e-2],
#                         [8.25e-2, 18.34e-2],
#                         [3.3e-2, 21.74e-2]])


    Rows = np.array([6,6,1])
    Cols = np.array([7,8,5])
   # Cols = np.array([12,8,5])


    TopLeft = np.array([[-7.42e-2, 20.04e-2],
                        [-8.25e-2, 18.34e-2],
                        [-3.3e-2, 21.74e-2]])*0.75

    LastRow = np.array([[-7.42e-2, 3.04e-2],
                        [-8.25e-2, 1.34e-2],
                        [-3.3e-2, 0]])*0.75


    LastCol = np.array([[7.42e-2, 20.04e-2],
                        [8.25e-2, 18.34e-2],
                        [3.3e-2, 21.74e-2]])*0.75

    RowStep = ((LastRow - TopLeft).T/np.where(Rows>1,Rows-1,1)).T
    ColumnStep = ((LastCol - TopLeft).T/np.where(Cols>1,Cols-1,1)).T

    ## Create grid in 2D

    n = 0
    nbeamlets = np.sum(Rows*Cols)
    grid["xyz"] = np.zeros((nbeamlets,2))
    for k in range(Rows.size):
        for i in range(Rows[k]):
            for j in range(Cols[k]):
                grid["xyz"][i+j*Rows[k]+n,:] = TopLeft[k] + i*RowStep[k] + j*ColumnStep[k]
        n = n + Rows[k]*Cols[k]

    ## Add a third dimension to the grid

    grid["xyz"] = np.array([np.zeros(nbeamlets),grid["xyz"][:,0],grid["xyz"][:,1]]).T

    ## Tilt the gridhalf

    A_tilt = np.array([[np.cos(tilt), 0, -np.sin(tilt)],
                       [0,1,0],
                       [np.sin(tilt), 0, np.cos(tilt)]])

    grid["xyz"] = np.matmul(A_tilt,grid["xyz"].T).T

    ## Second gridhalf as a mirror of the first

    grid["xyz"] = np.append(grid["xyz"],
                           np.array([grid["xyz"][:,0],grid["xyz"][:,1],-grid["xyz"][:,2]]).T,
                           axis=0)
    nbeamlets = 2*nbeamlets

    ## Define the beamlets direction in local coordinates

    grid["dxdydz"] = np.array([-np.ones(nbeamlets),
                               -grid["xyz"][:,1]/(fl_h - grid["xyz"][:,0]),
                               -grid["xyz"][:,2]/(fl_v - grid["xyz"][:,0])])

    grid["dxdydz"] = (grid["dxdydz"]/np.linalg.norm(grid["dxdydz"],axis=0)).T

    ## Reorientate and allocate the grid to its specific orientation and location

    A_thetasteer = np.array([[np.cos(thetasteer), 0, -np.sin(thetasteer)],
                       [0,1,0],
                       [np.sin(thetasteer), 0, np.cos(thetasteer)]])

    A_phisteer = np.array([[np.cos(phisteer),-np.sin(phisteer),0],
                           [np.sin(phisteer),np.cos(phisteer),0],
                           [0,0,1]])

    A_phicenter = np.array([[np.cos(phicenter),-np.sin(phicenter),0],
                      [np.sin(phicenter),np.cos(phicenter),0],
                      [0,0,1]])

    grid["xyz"] = np.matmul(A_thetasteer,grid["xyz"].T).T
    grid["xyz"] = np.matmul(A_phisteer,grid["xyz"].T).T

    grid["xyz"][:,0] = grid["xyz"][:,0] + Rcenter
    grid["xyz"][:,2] = grid["xyz"][:,2] + zcenter

    grid["xyz"] = np.matmul(A_phicenter,grid["xyz"].T).T

    grid["dxdydz"] = np.matmul(A_thetasteer,grid["dxdydz"].T).T
    grid["dxdydz"] = np.matmul(A_phisteer,grid["dxdydz"].T).T
    grid["dxdydz"] = np.matmul(A_phicenter,grid["dxdydz"].T).T

    return grid
