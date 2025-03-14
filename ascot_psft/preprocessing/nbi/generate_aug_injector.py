"""
Generate the AUG neutral beam injector geometry for BBNBI5.

IMPORTANT DISCLAIMER: The NBI geometry is in the "new" AUG
coordinate system (origin between sector 16 and 1). Make 
sure you don't mix them with inputs in the old AUG coordinate
system. Only the toroidal angle changes:

AUGnew_phi = AUGold_phi - 67.5 [deg]

Two beam geometry versions are available: 
    "old": original BBNBI version.
    "2017": beam geometry as defined in the  2017 geometry
    verification.

DISCLAIMER: "old" version does not stand for an old configuration
of the AUG NBIs. It is merely the version that it was originally
implemented in BBNBI. Therefore, it is recommended to use version
"2017" by all means, as it updates the NBI configuration with the
information available since 2017. Both versions are in the new
coordinate system.

J.F. Rivero-Rodriguez <jfrivero@us.es>
"""
import numpy as np
from scipy.constants import physical_constants as const

import a5py.ascot5io.nbi as a5nbi
import ascot_psft.marker.interpret as a5interpret


def generate(fn, inj=[1,2,3,4,5,6,7,8], version="2017", specie = "D", desc="AUG_NBI"):
    """
    Generate AUG injector.

    Args:
        fn : str
            Full path to the HDF5 file.
        inj : list of int, default = all
            List of AUG injectors from 1 to 8 to be generated.
        desc : str, optional
            Input description.
        version : {"old","2017"}
            NBI geometry version. "2017" is the default and recommended.
        specie : {"D","H"}
            Specie injected by the NBI. Deuterium "D" (default) or Protons "H".
    """

    ## Convert injector numbers into index

    inj = np.array(inj)-1 

    ## Define injectors geometry

    if version == "old":
        Rcenter = np.array([9.3234,
                            9.2490,
                            9.2490,
                            9.3234,
                            9.7407,
                            9.6371,
                            9.6371,
                            9.7407])
        phicenter = np.array([0.4569,
                              0.3560,
                              0.3560,
                              0.4569,
                              3.4764,
                              3.3799,
                              3.3799,
                              3.4764])-67.5*np.pi/180
        zcenter = np.array([0.6,
                            0.6,
                            -0.6,
                            -0.6,
                            0.6,
                            0.7,
                            -0.7,
                            -0.6])
        thetasteer = np.array([0.085505,
                               0.085505,
                               -0.085505,
                               -0.085505,
                               0.085505,
                               0.116141,
                               -0.116141,
                               -0.085505])
        phisteer = np.array([-0.057486,
                             -0.100898,
                             -0.100898,
                             -0.057486,
                             -0.086347,
                             -0.134245,
                             -0.134245,
                             -0.086347])
        focal_length_v = np.array([8.5,
                                   8.5,
                                   8.5,
                                   8.5,
                                   8.5,
                                   8.5,
                                   8.5,
                                   8.5])
        focal_length_h = np.array([6.5,
                                   6.5,
                                   6.5,
                                   6.5,
                                   6.5,
                                   6.5,
                                   6.5,
                                   6.5])
        tilt = 0.87*np.pi/180.
        div_h = np.array([0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099])*np.sqrt(2)
        div_v = np.array([0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099])*np.sqrt(2)
        div_halo_frac = np.array([0, #unknown data
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0])
        div_halo_h= np.array([1e-10, #unknown data
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10])
        div_halo_v= np.array([1e-10, #unknown data
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10])

    elif version == "2017":
        Rcenter = np.array([9.3569,
                            9.2821,
                            9.2822,
                            9.3569,
                            9.7759,
                            9.7020,
                            9.7032,
                            9.7759])
        phicenter = np.array([5.5620,
                              5.4611,
                              5.4610,
                              5.5630,
                              2.2981,
                              2.2021,
                              2.2020,
                              2.2981])
        zcenter = np.array([0.6,
                            0.6,
                            -0.6,
                            -0.6,
                            0.60093,
                            0.70136,
                            -0.70272,
                            -0.60155])
        thetasteer = np.array([0.08498,
                               0.08498,
                               -0.08498,
                               -0.08498,
                               0.08601,
                               0.11678,
                               -0.11849,
                               -0.08670])
        phisteer = np.array([-0.05727,
                             -0.10054,
                             -0.10148,
                             -0.05727,
                             -0.08760,
                             -0.13013,
                             -0.13171,
                             -0.08760])
        focal_length_v = np.array([8.5,
                                   8.5,
                                   8.5,
                                   8.5,
                                   11.94,
                                   11.94,
                                   11.94,
                                   11.94])
        focal_length_h = np.array([6.5,
                                   6.5,
                                   6.5,
                                   6.5,
                                   7.23,
                                   7.23,
                                   7.23,
                                   7.23])
        tilt = 0.87*np.pi/180.
        div_h = np.array([0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099])*np.sqrt(2)
        div_v = np.array([0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099,
                          0.0099])*np.sqrt(2)
        div_halo_frac = np.array([0, #unknown data
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0])
        div_halo_h= np.array([1e-10, #unknown data
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10])
        div_halo_v= np.array([1e-10, #unknown data
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10,
                              1e-10])

    else:
        raise ValueError("{} is not a recognised version of NBI geometry".format(version))
                          
    if specie == "D":
        anum = 2
        znum = 1
        mass = a5interpret.mass_amu("deuterium")
        energy = np.array([6.0e4,
                           6.0e4,
                           6.0e4,
                           6.0e4,
                           9.3e4,
                           9.3e4,
                           9.3e4,
                           9.3e4])
        power = np.array([2.5,
                          2.5,
                          2.5,
                          2.5,
                          2.5,
                          2.5,
                          2.5,
                          2.5])*1e6
        efrac = np.array([[0.65,0.25,0.1],
                          [0.65,0.25,0.1],
                          [0.65,0.25,0.1],
                          [0.65,0.25,0.1],
                          [0.62,0.29,0.09],
                          [0.62,0.29,0.09],
                          [0.62,0.29,0.09],
                          [0.62,0.29,0.09]])

    elif specie == "H":
        anum = 1
        znum = 1
        mass = a5interpret.mass_amu("hydrogen")
        energy = np.array([5.5e4,
                           5.5e4,
                           5.5e4,
                           5.5e4,
                           7.2e4,
                           7.2e4,
                           7.2e4,
                           7.2e4])
        power = np.array([1.8,
                          1.8,
                          1.8,
                          1.8,
                          1.43,
                          1.43,
                          1.43,
                          1.43])*1e6
        efrac = np.array([[0.51,0.3,0.19],
                          [0.51,0.3,0.19],
                          [0.51,0.3,0.19],
                          [0.51,0.3,0.19],
                          [0.38,0.35,0.27],
                          [0.38,0.35,0.27],
                          [0.38,0.35,0.27],
                          [0.38,0.35,0.27]])
                          
    else:
        raise ValueError("{} is not a recognised injection specie".format(specie))

    # Make nbi input of the selected injectors
    out = []
    for i in inj:
        ## Define beamlets grid
        grid = define_grid(Rcenter[i],phicenter[i],zcenter[i],
                           thetasteer[i],phisteer[i],
                           tilt,focal_length_h[i],focal_length_v[i])

        # out.append({})
        # out[-1]["id"] = i+1
        # out[-1]["nbeamlet"] = grid["xyz"][:,0].size
        # out[-1]["beamletx"] = grid["xyz"][:,0]
        # out[-1]["beamlety"] = grid["xyz"][:,1]
        # out[-1]["beamletz"] = grid["xyz"][:,2]
        # out[-1]["beamletdx"] = grid["dxdydz"][:,0]
        # out[-1]["beamletdy"] = grid["dxdydz"][:,1]
        # out[-1]["beamletdz"] = grid["dxdydz"][:,2]
        # out[-1]["div_h"] = div_h[i]
        # out[-1]["div_v"] = div_v[i]
        # out[-1]["div_halo_frac"] = div_halo_frac[i]
        # out[-1]["div_halo_h"] = div_halo_h[i]
        # out[-1]["div_halo_v"] = div_halo_v[i]
        # out[-1]["anum"] = anum
        # out[-1]["znum"] = znum
        # out[-1]["mass"] = mass
        # out[-1]["energy"] = energy[i]
        # out[-1]["efrac"] = efrac[i,:]
        # out[-1]["power"] = power[i]

        out.append(a5nbi.Injector(i+1, anum, znum, mass, energy[i], efrac[i,:],
                                power[i], div_h[i], div_v[i], div_halo_frac[i],
                                div_halo_h[i], div_halo_v[i], grid["xyz"][:,0].size,
                                grid["xyz"][:,0], grid["xyz"][:,1], grid["xyz"][:,2],
                                grid["dxdydz"][:,0], grid["dxdydz"][:,1], grid["dxdydz"][:,2]))
            
    return a5nbi.NBI.write_hdf5(fn, len(out), out, desc)


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

    ## Define beamlets position in local coordinates

    # Divide each gridhalf to 6 pieces: 2 of which create the main
    # diamond pattern, and 4 that add the four top rows.
    # Generate a rows-by-columns matrix of beamlets where
    # the extreme points are given as (front view):
    #
    # TopLeft ---------- LastCol
    #  \                       \
    #   \                       \
    #    \                       \
    #   LastRow ------------------   

    Rows = np.array([9,8,1,1,1,1])
    Cols = np.array([19,20,18,17,14,7])

    TopLeft = np.array([[-104.4e-3, 201.1e-3],
                        [-110.2e-3, 189.2e-3],
                        [-98.6e-3, 213.0e-3],
                        [-92.8e-3, 224.9e-3],
                        [-75.4e-3, 236.8e-3],
                        [-34.8e-3, 248.7e-3]])

    LastRow = np.array([[-104.4e-3, 10.7e-3],
                        [-110.2e-3, 22.6e-3],
                        [-98.6e-3, 0],
                        [-92.8e-3, 0],
                        [- 75.4e-3, 0],
                        [- 34.8e-3, 0]])


    LastCol = np.array([[104.4e-3, 201.1e-3],
                        [110.2e-3, 189.2e-3],
                        [98.6e-3, 213.0e-3],
                        [92.8e-3, 224.9e-3],
                        [75.4e-3, 236.8e-3],
                        [34.8e-3, 248.7e-3]])

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
