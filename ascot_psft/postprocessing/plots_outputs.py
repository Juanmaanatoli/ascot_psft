"""
Sets of scripts for ASCOT plotting. The scripts for any plots related to FILD are in a separate file.

@authors: Javier Hidalgo-Salaverri <jhsalaverri@us.es>; Lina Velarde <lvelarde@us.es>

List of functions:
    plot_vessel_Rz. Plot the current vessel. Mostly used as an aux for other plotting routines.
    plot_orbit. Plot a series of simulated orbits.
    plot_markers_input_points. Plot initial markers distribution with points.
    plot_markers. Plot the markers from either the marker input or from the simulation. 
    compare_markers_1D. Compare the markers (from either the marker input or from the simulation)
                            between several scenarios, using 1D histograms.
    compare_markers_2D. Compare the markers (from either the marker input or from the simulation)
                            between two different scenarios, using 2D histograms.
    plot_hists. Plot histograms. Options so far: R vs pitch, pitch, R vs phi.
    plot_markers_velspace. Plot the VS of the markers.
        
"""

import numpy as np
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    print('Could not import tqdm')
from a5py import Ascot
import matplotlib
import matplotlib.cm as cm
import matplotlib
# # matplotlib.use('TkAgg')
import matplotlib.pyplot as plt; plt.ion()
import matplotlib.colors as colors
from scipy import interpolate as si
from ascot_psft.code_utilities import pol2cart
from ascot_psft.aux_files.common import read_bf
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap
from ascot_psft.aux_files.common import plot_vessel_tor
from ascot_psft.aux_files.vessel import vessel_pol as plot_vessel_pol

fs = 18 #fontsize for all figures and labels
plt.rc('xtick', labelsize=fs)
plt.rc('ytick', labelsize=fs)
cmapidl = LinearSegmentedColormap.from_list(
    'mycmap', ['black', 'blue', 'red', 'yellow', 'white'], N=256)
cmapinf = plt.cm.inferno
cmapbwr = plt.cm.bwr
cmapseismic = plt.cm.seismic


colorlist = ['blue', 'red', 'green', 'black', 'darkorange','cyan', 'purple', 'grey',
             'darkgreen', 'pink', 'peru', 'magenta', 'lime', 'gold', 'turquoise', 
             'indianred','slategray','brown']


def plot_orbit(file, orbits_id, qid_run: str = 'active', qid_marker: str = 'active',
               plot_vv: bool = True, qid_wall: str = 'active', ax = None, view = 'pol',
               device = 'AUG'):
    """
    Plot a series of simulated orbits.
    @author: Javier Hidalgo-Salaverri <jhsalaverri@us.es>

    :param file: path to the .h5 file
    :param orbits_id: List of the orbit ids to plot.
    :param qid_run: run qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    :param qid_marker: markers qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    :param plot_vv: plot the vacuum vessel 
    :param qid_wall: wall qid, only if plot_vv == True
    :param view: Plot in the RZ ('pol') or XY ('tor') plane
    :param ax: axis where to plot, if None a new one wil be created
    :param device: 'AUG' or 'MU' implemented so far
    """
    a5 = Ascot(file)
    try:
        int(qid_run) # If it is a number, will not raise an error
        qid_run = 'q'+qid_run
    except ValueError:
        pass
    try:
        int(qid_marker) # If it is a number, will not raise an error
        qid_marker = 'q'+qid_marker
    except ValueError:
        pass
    
    all_orbits = a5.data.__getattribute__(qid_run).orbit.read()
    all_pos0 = a5.data.marker.__getattribute__(qid_marker).read()

    if ax is None:
        fig, ax = plt.subplots()

    for orbit_id in orbits_id:
        index = all_orbits['ids'] == orbit_id
        r0 =  all_pos0['r'][orbit_id]
        z0 =  all_pos0['z'][orbit_id]
        phi0 = all_pos0['phi'][orbit_id]
        print('Careful! Is phi in radians here? Need to be checked for toroidal plotting!')
        r =  all_orbits['r'][index]
        z =  all_orbits['z'][index]
        phi = all_orbits['phi'][index]

        if view == 'pol':
            line, = ax.plot(r,z, label = str(orbit_id))
            ax.plot(r0,z0, 'x', color = line.get_color(), markersize = 15)
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')
            if plot_vv:
                if device == 'MU':
                    plot_vessel_pol.MASTu(ax = ax)
                elif device == 'AUG':
                    plot_vessel_pol.AUG(ax = ax)
                else:
                    raise NameError('Device can only be "AUG" or "MU" at the moment')
        
        elif view == 'tor':
            x0, y0 = pol2cart(r0,phi0)
            x, y = pol2cart(r,phi)
            line, = ax.plot(x,y, label = str(orbit_id))
            ax.plot(x0,y0, 'x', color = line.get_color(), markersize = 15)
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            if plot_vv:
                if device == 'MU':
                    plot_vessel_tor(device = 'MU', ax = ax, cwall = 'k')
                elif device == 'AUG':
                    plot_vessel_tor(device = 'AUG', ax = ax, cwall = 'k')
                else:
                    raise NameError('Device can only be "AUG" or "MU" at the moment')
        else:
            raise NameError('View can only be "tor" or "pol"')

    ax.legend()
    ax.set_aspect(1)
    plt.tight_layout()
    plt.show(block = False)


def plot_markers_input_points(file, qid_marker: str = 'active', plot_vv: bool = True, qid_wall: str = 'active', 
                              ax = None, view: str = 'pol'):
    """
    Plot initial markers distribution. Each marker will be a point.
    Can only plot a poloidal view.

    @author: Javier Hidalgo-Salaverri <jhsalaverri@us.es>

    :param file: path to the .h5 file
    :param qid_marker: markers qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    :param plot_vv: plot the vacuum vessel 
    :param qid_wall: wall qid, only if plot_vv == True
    :param ax: axis where to plot, if None a new one wil be created
    :param view: Plot in the RZ ('pol') or XY ('tor') plane
    """
    a5 = Ascot(file)
    try:
        int(qid_marker) # If it is a number, will not raise an error
        qid_marker = 'q'+qid_marker
    except ValueError:
        pass
    marker = a5.data.marker.__getattribute__(qid_marker).read()
    r = marker['r']
    z = marker['z']
    # Plot
    if ax is None:
        fig, ax = plt.subplots()
    # if plot_vv: 
    #     plot_vessel_Rz(file = file, qid_wall = qid_wall, view = view, ax = ax)

    ax.plot(r,z, 'o')
    ax.set_xlabel('r [m]')
    ax.set_ylabel('z [m]')
    ax.set_aspect(1)
    plt.tight_layout()
    plt.show(block = False)


def plot_markers(path, markertype: str = 'input', qid: str = 'active', qid_bf: str = 'active', state: str = 'end', ec: str = 'wall', 
                 mode: str = 'prt', norm: str = 'log', vmax: int = None, cmap: str = cmapidl, view: str = 'both', 
                 flag_phitheta: bool = False, flag_rho: bool = False, title: str = None, device: str = 'MU', 
                 rfild: float = 1.5, plotwall2D: bool = False, plotRMPs: bool = False):
    """
    Routine to plot the markers either from the simulation or from the input. 
    Can be inistate or endstate. Endcondition can be specified.
    
    @author: Lina Velarde <lvelarde@us.es>

    Inputs:
    path. Path to the ascot file.
    markertype. 'input' for input markers, 'sim' for markers read from the specified simulation.
    qid. If markertype=='input', this is the marker's qid. If markertype=='sim', this is the simulation's qid.
        By default, the active marker input/sim will be used, but it can also be 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    
        INPUTS USED ONLY IF markertype=='input':
    qid_bf. Qid of the bfield. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
        
        INPUTS USED ONLY IF markertype=='sim':
    state. Default is 'end' for endstate, can be 'ini'.
    ec. It is the end condition to plot. If None, will plot all.
        Default is wall, but it can be "none" for aborted particles, or tmax for confined ones.
    mode. 'prt' for particle position, 'gc' for gyrocenter position.
    
        COMMON INPUTS:
    norm. 'log' for logarithmic scale, None for linear scale.
    vmax. vmax for plotting.
    cmap. cmap for plotting.
    view. Poloidal 'pol' or toroidal 'tor' view. Can be 'both'.
    flag_phitheta. If True, will also plot the phi vs theta view.
    flag_rho. To plot an additional figure with the rho coordinate on top of the r-z plot.
    title. If included, will be displayed.
    device. 'MU' or 'AUG'.
        If MU, you can choose:
            rfild. FILD position, for plotting in poloidal view.
            plotwall2D. If the ASCOT option for 2D plotting should be used. This 
                view has some more details than the imported from vessel one.
            plotRMPs. To include RMPs in the plot.
    """
	
    h5 = Ascot(path)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
        
    if markertype == 'input':
        print('The markers from the input will be plotted!!')
        marker = h5.data.marker.__getattribute__(qid)
        print("The marker's description is:",marker.get_desc())
        
        try:
            int(qid_bf) # If it is a number, will not raise an error
            qid_bf = 'q'+qid_bf
        except ValueError:
            pass
        bf = h5.data.bfield.__getattribute__(qid_bf).read()
        r = marker.read()['r']
        z = marker.read()['z']
        phi = marker.read()['phi']*np.pi/180
    
    elif markertype == 'sim':
        print('The markers from the simulation will be plotted!!')
        sim = h5.data.__getattribute__(qid)
        print("The run's description is:",sim.get_desc())
        print("The state is:", state)
        print("The endcondtion is:", ec)
        
        bf = sim.bfield.read()
        r, z, phi = sim.getstate('r', 'z', 'phi', mode=mode, endcond=ec, state=state)
    
    # from here, everything is common for both sources of markers
    # (python knows phi is in degrees for the sim markers)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    # read bfield
    nr, rvec, nz, zvec, nphi, phivec, psi, psi1, psi0, bf2D = read_bf(bf)

    if norm == 'log':
        norm = colors.LogNorm()
        cwall = 'k'
    elif norm == None:
        cwall = 'w'
    # Marker: r vs z
    if view == 'pol' or view == 'both':
        # if view both, will plot the poloidal and toroidal views in the same figure
        if view == 'both':
            fig, axes = plt.subplots(1,2, figsize=(14,6))
            ax = axes[0]
        else:
            fig, ax = plt.subplots()
        if device == 'MU':
            if plotwall2D:
                sim.wall.active.plotRz()
            plot_vessel_pol.MASTu(ax=ax,color=cwall,FILD_R=rfild)
            range_rz = [[0.,2.3],[-2,2]] # range for histogram calculation
            xlim = [0,2] # limits for plotting
            ylim = [-2.1,2.1] # limits for plotting
        elif device == 'AUG':
            if plotwall2D:
                sim.wall.active.plotRz()
            else:
                plot_vessel_pol.AUG(ax=ax,color=cwall)
            range_rz = [[0.1,2.3],[-2,2]]
            xlim = range_rz[0]
            ylim = range_rz[1]
            print('Not fully implemented. Please correct ranges and wall plotting')
        else:
            raise NameError('Only "AUG" or "MU" devices implemented at the moment')
        h1 = ax.hist2d(r, z, bins=400, range=range_rz,
                        cmap=cmap, norm=norm, vmax=vmax)
        # ax.set_facecolor('k')
        cbar = plt.colorbar(h1[3])
        cbar.set_label('Counts',labelpad=fs-3,fontsize=fs)
        cont = ax.contour(rvec, zvec, psi.T, levels=psi1, colors="red", linestyles='solid')
        # h5.input_plotrhocontour(rho=1,axes=ax)
        if device == 'MU' and plotRMPs is True:
            Ru = 1.43
            Rl = 1.5
            zu = 0.77
            zl = 0.5
            ax.plot([Ru,Rl],[zu,zl],color='lime',linewidth=2)
            ax.plot([Ru,Rl],[-zu,-zl],color='lime',linewidth=2)    
        ax.set_xlabel('r [m]',fontsize=fs)
        ax.set_ylabel('z [m]',fontsize=fs)
        # plt.axis('scaled')
        ax.set_aspect(aspect='equal', adjustable='box')
        plt.tight_layout()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        if flag_rho == True:
            plt.figure()
            rho = np.sqrt((psi - psi0)/(psi1-psi0))
            if device == 'MU':
                vessel.MASTu(FILD_R=rfild,color='k')
            elif device == 'AUG':
                vessel.AUG(color='k')
            else:
                raise NameError('Only "AUG" or "MU" devices implemented at the moment')  
            cont = plt.contour(rvec, zvec, rho.T, levels = np.linspace(0,1.25,60))
            # cont = plt.contour(rvec, zvec, rho.T, levels = np.linspace(0.98,1.2,20), 
                    # cmap=plt.cm.jet)
            cbar = plt.colorbar()
            cbar.set_label('rhopol',labelpad=10,fontsize=fs)
            h1 = plt.hist2d(r, z, bins=400, range=range_rz, norm=norm, 
                            vmax=vmax, cmap=plt.cm.Oranges)
            cbar = plt.colorbar(h1[3])
            cbar.set_label('Counts',labelpad=10,fontsize=fs)
            cont = plt.contour(rvec, zvec, psi.T, levels=psi1, colors="red", linestyles='solid')
            plt.xlabel('r [m]',fontsize=fs)
            plt.ylabel('z [m]',fontsize=fs)
            plt.axis('scaled')
            plt.tight_layout()
        
    # Marker: x vs y
    if view == 'tor' or view == 'both':
        # if view both, will plot the poloidal and toroidal views in the same figure
        if view == 'both':
            ax = axes[1]
        else:
            fig, ax = plt.subplots()
        xlim, ylim = plot_vessel_tor(device = device, ax = ax, cwall = 'k')
        range_xy = [xlim, ylim] # will use same limits for plotting and hist calculation
        h1 = ax.hist2d(x, y, bins=400, range=range_xy,
                        cmap=cmap, norm=colors.LogNorm())
        # ax.set_facecolor('k')
        cbar = plt.colorbar(h1[3])
        cbar.set_label('Counts',labelpad=fs-3,fontsize=fs)
        if device == 'MU' and plotRMPs is True:
            rRMP = 1.44
            width = 0.2
            phi_vals = np.linspace(24*np.pi/180, 2*np.pi+24*np.pi/180, 9)[:-1] # 4 rectangles equispaced in phi
            for phi in phi_vals:
                x1 = rRMP*np.cos(phi-width)
                x2 = rRMP*np.cos(phi+width)
                y1 = rRMP*np.sin(phi-width)
                y2 = rRMP*np.sin(phi+width)
                plt.plot([x1,x2],[y1,y2],color='lime',linewidth=2)
        xlim, ylim = plot_vessel_tor(device = device, ax = ax, cwall = 'k')
        ax.set_xlabel('x [m]',fontsize=fs)
        ax.set_ylabel('y [m]',fontsize=fs)
        # plt.axis('scaled')
        ax.set_aspect(aspect='equal', adjustable='box')
        plt.tight_layout()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)


    if flag_phitheta:
        thetamod, phimod = sim.getstate('thetamod', 'phimod', mode=mode, state=state,endcond=ec)
        # convert to [-180 to 180]
        thetamod2 = thetamod.to_value()-180
        phimod2 = phimod.to_value()-180
        phimod2[phimod2<0] +=360.0
        phimod2 -= 180
        fig, ax = plt.subplots()
        h1 = ax.hist2d(phimod.to_value(), thetamod.to_value(),range=([0,360],[0,360]), bins=400,
                    cmap=cmap, norm=colors.LogNorm(), vmax=vmax)
        plt.title('Before conversion, [0,360]deg')
        plt.tight_layout()
        fig, ax = plt.subplots()
        h1 = ax.hist2d(phimod2, thetamod2, bins=400,range=([-180,180],[-180,180]),
                    cmap=cmapidl, norm=colors.LogNorm(), vmax=vmax)
        # ax.set_facecolor('k')
        cbar = plt.colorbar(h1[3])
        cbar.set_label('Counts',labelpad=10,fontsize=fs)
        plt.axvline(x = 15, color='r', linewidth=1, alpha=0.5, linestyle='--')
        plt.xlabel('phi [deg]',fontsize=fs)
        plt.ylabel('theta [deg]',fontsize=fs)
        plt.axis('scaled')
        if title:
            plt.title(title + '(unfiltered)', fontsize=fs+2)
        if device == 'MU':
            plt.text(20, -175, 'FILD', c='r', fontsize=18)
        plt.tight_layout()


def compare_markers_1D(pathlist, qidlist,  markertype = 'sim', variables = ['r','z','phi'], state = 'end', 
                       ec = 'wall', mode = 'prt', qid_bf: str = None, title: str = None, labellist = None,
                       cbarlabel: str = 'Counts', porcrz = 0.55, porcxy = 0.55, porcpt =  0.65,
                       device = 'MU', rfild = 1.5, plotwall2D = False):
    """
    Routine to compare the markers from several simulations. 1D  plots.
    Will just plot 1D histograms for each simulation, in the same figure.
    If the simulation's markers are being used, can be inistate or endstate.
    This can also be used to compare the birth distribution by setting the markertype to 'input'.

    @author: Lina Velarde <lvelarde@us.es>
    
    Inputs:
    pathlist. List of paths to the files.
    markertype. 'input' for input markers, 'sim' for markers read from the specified simulation.
    qidlist. List of qids to compare. Can only be the 10 digit ID number, without the Q. 
            If markertype=='input', these are the marker's qids. If markertype=='sim', 
             these are the simulation's qids.
    
        INPUTS USED ONLY IF markertype=='sim':
    state. Default is 'end' for endstate, can be 'ini'.
    ec. It is the end condition to plot. 
        Default is wall, but it can be "none" for aborted particles, or tmax for confined ones.
    mode. 'prt' for particle position, 'gc' for gyrocenter position.
    
        INPUTS USED ONLY IF markertype=='input':
    qid_bf. Qid of the bfield. If not specified, the active one will be used.
        
        COMMON INPUTS:
    variables. List of variables to plot. Can be 'r', 'z', 'phi', 'thetamod', 'phimod', 'ekin', 'pitch', 'rho'.
    title. If included, will be displayed.
    labellist. If None, the descript of the input/run will be used.
    device. 'MU' or 'AUG'.
        If MU, you can choose:
            rfild. FILD position, for plotting in poloidal view.
            plotwall2D. True if the ASCOT option for 2D plotting should be used. This 
                view has some more details than the imported from vessel one (but it's not great).
    """
    if type(pathlist)==str:
        print('Only one file specified. Will assume all runs are in the same file')
        pathlist = [pathlist]*len(qidlist)
    elif len(pathlist) == 1:
        print('Only one file specified. Will assume all runs are in the same file')
        pathlist = [pathlist[0]]*len(qidlist)
    elif len(pathlist) != len(qidlist):
        print('The number of paths and qids must be the same')
        return
    
    # Create as many subplots as variables. If up to 3, plot in the same row.
    if len(variables) <= 3:
        fig, axes = plt.subplots(1,len(variables), figsize=(6*len(variables),8))
    else:
        fig, axes = plt.subplots(2, int(np.ceil(len(variables) / 2)), figsize=(6*int(np.ceil(len(variables) / 2)), 16))
        
    # flatten axes list
    if len(variables) > 1:
        axes = axes.flatten()
    else:
        axes = [axes] #avoid crash later
    
    if markertype == 'input':
        print('The markers from the input will be compared!!')
    elif markertype == 'sim':
        print('The markers from the simulation will be plotted!!')
    else:
        print('Markertype not recognised. Exiting.')
        return

    # Go through the list of simulations plotting the histograms
    for ii, path in enumerate(pathlist):
        h5 = Ascot(path)
        print('--------------------------------------------')
        print("The file's name is:",path)
        # create dictionary to store the data
        data = {}
        
        if markertype == 'input':
            print("The marker's qid is:",qidlist[ii])

            sim = h5.data.marker["q"+qidlist[ii]]

            print("The marker's description is:",sim.get_desc())
            
            bf = h5.data.bfield.active.read()
            for var in variables:
                if var == 'thetamod' or var == 'phimod' or var == 'pitch' or var == 'ekin':
                    print('The ekin, pitch, phimod or theta histograms cannot be plotted for the input markers')
                    continue
                elif var == 'phi':
                    data[var] = sim.read()[var]
                # if the variable is rho, will be calculated from r and z later
                # for now, create an empty array
                elif var == 'rho':
                    print('The active bfield will be used for the rho calculation!')
                    data[var] = np.array([])
                    if 'z' not in variables:
                        data['z'] = sim.read()['z']
                else:
                    data[var] = sim.read()[var]
    
        elif markertype == 'sim':
            print("The run's qid is:",qidlist[ii])
            try:
                sim = h5.data["run_"+qidlist[ii]]
            except AttributeError:
                print('CAREFUL! THE QID DOES NOT CORRESPOND TO AN ASCOT RUN')
                print('Please run again setting the markertype to input')
                return

            print("The run's description is:",sim.get_desc())
            print("The endcondtion is:", ec)
            print("The state is:", state)
            
            bf = sim.bfield.read()
            # qid_bf = sim.bfield.get_qid()
            # print('Initialising bfield')
            # h5.input_init(bfield = qid_bf)
            for var in variables:
                # if the variable is rho, will be calculated from r and z later
                # for now, create an empty array
                if var == 'rho':
                    data[var] = np.array([])
                    if 'z' not in variables:
                        data['z'] = sim.getstate('z', mode=mode, state=state,endcond=ec)
                else:
                    data[var] = sim.getstate(var, mode=mode, state=state,endcond=ec)
        
            if 'thetamod' in variables:
            # calculate the phi,theta histogram
            # convert to [-180 to 180]
                data['thetamod'] = data['thetamod'].to_value()-180
                
            if 'phimod' in variables:
                data['phimod'] = data['phimod'].to_value()-180
                data['phimod'][data['phimod']<0] +=360.0
                data['phimod'] -= 180
        
        if 'rho' in variables:
            nr, rvec, nz, zvec, nphi, phivec, psi, psi1, psi0, bf2D = read_bf(bf)
            # psi is in [rvec, zvec] coordinates
            # need to calculate the psi for each marker knowing r and z
                        # interpolate psi to get the value for each marker
            f = si.RegularGridInterpolator((rvec,zvec),psi)
            points = np.vstack((data['r'],data['z'])).T
            psi_interp = f(points)
            data['rho'] = np.sqrt((psi_interp - psi0)/(psi1 - psi0))
        
        
        # Get label for the legend for each marker
        if labellist is None:
            lab = sim.get_desc()
        else:
            lab = labellist[ii]
                
        # with the variables read, plot the histograms
        for jj, var in enumerate(variables):
            # calculate the histogram
            hist = np.histogram(data[var], bins=100)
            # plot the histogram
            axes[jj].plot(hist[1][:-1], hist[0], label=lab, color=colorlist[jj])
            # in the last iteration, set the labels
            if ii == len(pathlist)-1:
                axes[jj].set_xlabel(var, fontsize=fs)
                axes[jj].set_ylabel('Counts', fontsize=fs)
                axes[jj].legend(loc='best')
    
    plt.tight_layout()


def compare_markers_2D(pathlist, qidlist,  markertype = 'sim', state = 'end', ec = 'wall', mode = 'prt', 
                    view: str = 'both', flag_phitheta = False, cmap=cmapbwr, cwall='gray', title: str = None, 
                    cbarlabel: str = 'Counts', porcrz = 0.55, porcxy = 0.55, porcpt =  0.65,
                    device = 'MU', rfild = 1.5, plotwall2D = False):
    """
    Routine to compare the markers from ONLY TWO different simulations. 2D  plots.
    Will calculate the histograms and densities for each simulation and plot the difference between the two.
    The difference is calculated as the second - the first sim.
    If the simulation's markers are being used, can be inistate or endstate.
    But this can also be used to compare the birth distribution by setting the markertype to 'input'.

    @author: Lina Velarde <lvelarde@us.es>
    
    Inputs:
    pathlist. List of paths to the files.
    markertype. 'input' for input markers, 'sim' for markers read from the specified simulation.
    qidlist. List of qids to compare. Can only be the 10 digit ID number, without the Q. 
            If markertype=='input', these are the marker's qids. If markertype=='sim', 
            these are the simulation's qids.
    
        INPUTS USED ONLY IF markertype=='sim':
    state. Default is 'end' for endstate, can be 'ini'.
    ec. It is the end condition to plot. 
        Default is wall, but it can be "none" for aborted particles, or tmax for confined ones.
    mode. 'prt' for particle position, 'gc' for gyrocenter position.
    
        INPUTS USED ONLY IF markertype=='input':
    qid_bf. Qid of the bfield. If not specified, the active one will be used.
        
        COMMON INPUTS:
    view. Poloidal 'pol' or toroidal 'tor' view. Can be 'both'.
    flag_phitheta. If True, will also plot the phi vs theta view.
    cmap. cmap for plotting.
    cwall. color for wall plotting.
    title. If included, will be displayed.
    cbarlabeel. Label for the colorbar. Default is 'Counts'.
    porc (rz, xy, pt). Percentages to adjust the colorbar limits.
    device. 'MU' or 'AUG'.
        If MU, you can choose:
            rfild. FILD position, for plotting in poloidal view.
            plotwall2D. True if the ASCOT option for 2D plotting should be used. This 
                view has some more details than the imported from vessel one (but it's not great).
    """
    if type(pathlist)==str:
        print('We will assume both runs are in the same file')
        pathlist = [pathlist, pathlist]     
    elif len(pathlist) == 1:
        print('We will assume both runs are in the same file')
        pathlist = [pathlist[0], pathlist[0]]    
    elif len(pathlist) != len(qidlist):
        print('The number of paths and qids must be the same')
        return
    hist_values = []
    dens_values = []
    # Go through the list of simulations calculating the histograms and densities
    for ii, path in enumerate(pathlist):
        h5 = Ascot(path)
        print("The file's name is:",path)
        
        if markertype == 'input':
            print('The markers from the input will be compared!!')
            print("The marker's qid is:",qidlist[ii])

            sim = h5.data.marker["q"+qidlist[ii]]

            print("The marker's description is:",sim.get_desc())
            
            bf = h5.data.bfield.active.read()
            r = sim.read()['r']
            z = sim.read()['z']
            phi = sim.read()['phi']*np.pi/180
            
            flag_phitheta = False
            print('The phi vs theta view cannot be plotted for the input markers')
    
        elif markertype == 'sim':
            print('The markers from the simulation will be plotted!!')
            print("The run's qid is:",qidlist[ii])
            sim = h5.data["run_"+qidlist[ii]]

            print("The run's description is:",sim.get_desc())
            print("The endcondtion is:", ec)
            print("The state is:", state)
            
            bf = sim.bfield.read()
            # qid_bf = sim.bfield.get_qid()
            # print('Initialising bfield')
            # h5.input_init(bfield = qid_bf)
            r, z, phi = sim.getstate('r', 'z', 'phi', mode=mode, endcond=ec, state=state)
            thetamod, phimod = sim.getstate('thetamod', 'phimod', mode=mode, state=state,endcond=ec)
        
            # calculate the phi,theta histogram
            # convert to [-180 to 180]
            thetamod2 = thetamod.to_value()-180
            phimod2 = phimod.to_value()-180
            phimod2[phimod2<0] +=360.0
            phimod2 -= 180
            histpt = np.histogram2d(phimod2, thetamod2, bins=[49,50], range=([-180,180],[-180,180]))
            denspt = histpt[0]

        # calculate the r,z histogram
        histrz = np.histogram2d(r, z, bins=[49,50], range=[[0.2,1.7],[-1.5,1.5]])
        densrz = histrz[0]

        # calculate the x,y histogram
        # python knows phi is in degrees for the sim markers
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        histxy = np.histogram2d(x, y, bins=[49,50], range=[[-1.5,1.5],[-1.5,1.5]])
        densxy = histxy[0]

        if markertype == 'input':
            dens_values.append([densrz, densxy])
            hist_values.append([histrz, histxy])
        else:
            dens_values.append([densrz, densxy, denspt])
            hist_values.append([histrz, histxy, histpt])


    # Marker: r vs z
    if view == 'pol' or view == 'both':
        # if view both, will plot the poloidal and toroidal views in the same figure
        if view == 'both':
            fig, axes = plt.subplots(1,2, figsize=(14,6))
            ax = axes[0]
        else:
            fig, ax = plt.subplots()
        # get the r and z centers from the rz hist of the first run (they're the same)
        rcenter = hist_values[0][0][1][:-1] + (hist_values[0][0][1][1:] - hist_values[0][0][1][:-1])/2
        zcenter = hist_values[0][0][2][:-1] + (hist_values[0][0][2][1:] - hist_values[0][0][2][:-1])/2
        # get the density for the rz hist for each run
        dens1 = dens_values[0][0]
        dens2 = dens_values[1][0]
        # prepare the norm for the plot
        print('vmax = ', (dens2.T-dens1.T).max())
        print('vmin = ', (dens2.T-dens1.T).min())
        vmax = (dens2.T-dens1.T).max()
        vmin = (dens2.T-dens1.T).min()
        vmaxporc = vmax*porcrz # 10%
        vminporc = vmin*porcrz # 10%
        norm = colors.TwoSlopeNorm(vmin=(dens2.T-dens1.T).min()-vminporc, 
                vcenter=0, vmax=(dens2.T-dens1.T).max()-vmaxporc)
        hidfrz = ax.pcolormesh(rcenter, zcenter, dens2.T-dens1.T,
                        cmap=cmap, norm = norm,zorder=1)
        cbar = plt.colorbar(hidfrz)
        cbar.set_label(cbarlabel,labelpad=fs-3,fontsize=fs)

        # plot separatrix
        nr, rvec, nz, zvec, nphi, phivec, psi, psi1, psi0, bf2D = read_bf(bf)
        print('A reference separatrix will be plotted. Please check it if this is important for you.')
        cont = ax.contour(rvec, zvec, psi.T, levels=psi1, colors="red", linestyles='solid')
        # h5.input_plotrhocontour(rho=1,axes=ax)
        if device == 'MU':
            if plotwall2D:
                sim.wall.active.plotRz()
            plot_vessel_pol.MASTu(ax=ax,color=cwall,FILD_R=rfild)
        elif device == 'AUG':
            plot_vessel_pol.AUG(ax=ax,color=cwall)
        else:
            raise NameError('Only "AUG" or "MU" devices implemented at the moment')
        plt.xlabel('r [m]',fontsize=fs)
        plt.ylabel('z [m]',fontsize=fs)
        plt.axis('scaled')
        plt.title(title,fontsize=fs+2)
        plt.xlim([rcenter.min(),rcenter.max()])
        plt.ylim([zcenter.min(),zcenter.max()])
        plt.tight_layout()
        
    # Marker: x vs y
    if view == 'tor' or view == 'both':
        # if view both, will plot the poloidal and toroidal views in the same figure
        if view == 'both':
            ax = axes[1]
        else:
            fig, ax = plt.subplots()
        # get the x and y centers from the xy hist of the first run (they're the same)
        xcenter = hist_values[0][1][1][:-1] + (hist_values[0][1][1][1:] - hist_values[0][1][1][:-1])/2
        ycenter = hist_values[0][1][2][:-1] + (hist_values[0][1][2][1:] - hist_values[0][1][2][:-1])/2
        dens1 = dens_values[0][1]
        dens2 = dens_values[1][1]
        print('vmax = ', (dens2.T-dens1.T).max())
        print('vmin = ', (dens2.T-dens1.T).min())
        vmax = (dens2.T-dens1.T).max()
        vmin = (dens2.T-dens1.T).min()
        vmaxporc = vmax*porcxy # 10%
        vminporc = vmin*porcxy # 10%
        norm = colors.TwoSlopeNorm(vmin=(dens2.T-dens1.T).min(), vcenter=0, vmax=(dens2.T-dens1.T).max())
        norm = colors.TwoSlopeNorm(vmin=(dens2.T-dens1.T).min()-vminporc, vcenter=0, vmax=(dens2.T-dens1.T).max()-vmaxporc)
        hdif = ax.pcolormesh(xcenter, ycenter, dens2.T-dens1.T,
                        cmap=cmap, norm=norm)
        cbar = plt.colorbar(hdif)
        cbar.set_label(cbarlabel,labelpad=fs-3,fontsize=fs)
        xlim, ylim = plot_vessel_tor(device = device, ax = ax, cwall = cwall)
        if device == 'MU':
            plt.plot([0, 1.5], [0,0.41], color='k', linewidth=1, alpha=0.7, linestyle='--')
            plt.text(0.5, 0.4, 'FILD', c='k', fontsize=18)            
        plt.xlabel('x [m]',fontsize=fs)
        plt.ylabel('y [m]',fontsize=fs)
        plt.axis('scaled')
        plt.title(title,fontsize=fs)
        plt.tight_layout()
        plt.xlim(xlim)
        plt.ylim(ylim)

    # Marker: phi vs theta
    if flag_phitheta:
        phicenter = hist_values[0][2][1][:-1] + (hist_values[0][2][1][1:] - hist_values[0][2][1][:-1])/2
        thetacenter = hist_values[0][2][2][:-1] + (hist_values[0][2][2][1:] - hist_values[0][2][2][:-1])/2
        dens1 = dens_values[0][2]
        dens2 = dens_values[1][2]
        print('vmax = ', (dens2.T-dens1.T).max())
        print('vmin = ', (dens2.T-dens1.T).min())  
        vmax = (dens2.T-dens1.T).max()
        vmin = (dens2.T-dens1.T).min()
        vmaxporc = vmax*porcpt # 10%
        vminporc = vmin*porcpt # 10%
        fig, ax = plt.subplots()
        norm = colors.TwoSlopeNorm(vmin=(dens2.T-dens1.T).min()-vminporc, vcenter=0, vmax=(dens2.T-dens1.T).max()-vmaxporc)
        # norm = colors.TwoSlopeNorm(vmin=-200, vcenter=0, vmax=600)
        # norm = colors.SymLogNorm(linthresh=0.1, linscale=0.5,vmin=-700, vmax=700, base=10)
        # norm = colors.AsinhNorm(linear_width=0.1, vmin=None, vmax=None)
        hdif = plt.pcolormesh(phicenter, thetacenter, dens2.T-dens1.T,
                        cmap=cmap, norm=norm)
        cbar = plt.colorbar(hdif)
        cbar.set_label(cbarlabel,labelpad=fs-2,fontsize=fs)
        plt.xlabel('phi [deg]',fontsize=fs)
        plt.ylabel('theta [deg]',fontsize=fs)
        plt.axis('scaled')
        if title:
            plt.title(title + '(unfiltered)', fontsize=fs+2)
        if device == 'MU':
            plt.axvline(x = 15, color='k', linewidth=1, alpha=0.7, linestyle='--')
            plt.text(20, -175, 'FILD', c='k', fontsize=18)
        plt.tight_layout()


def plot_hists(path, qid: str = 'active', plotxi = True, plotRxi = False, plotRphi = False,
               state='end', ec=None, mode='prt', ax=None, title=None, cmap = cmapidl, norm = None,
               flag_getval = False, pitchindeg = True, phi360 = True):
    """
    Routine to plot histograms of the final or initial state. 
    Options so far:
        R vs pitch
        pitch
        R vs phi

    @author: Lina Velarde <lvelarde@us.es>
    
    Inputs:
    path. Path to the ascot file.
    qid. Run's qid. . Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    flag_plotxi. Plots a histogram in xi
    flag_plotRxi. Plots a histogram in R vs xi
    flag_getval. If True, it returns the R and pitch ranges from the bin that
        has the maximum amount of particles.
    state. Can be either 'ini' or 'end'.
    ec. End condition of the markers. Typicall 'wall' or 'tlim'. If 'none', 
        it will return aborted markers. But, if None, it will return all markers.
    pitchindeg. Will plot the pitch in degrees. If false, will not convert.
    phi360. Will convert phi to a [0,360] range. If false, will not convert ([-180,180] range).
    """

    h5 = Ascot(path)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)
        
    # initialise bfield
    qid_bf = sim.bfield.get_qid()
    print('Initialising bfield')
    h5.input_init(bfield = qid_bf)

    # --- Get variables
    pitch, R, phi = sim.getstate('pitch', 'r', 'phi', mode=mode, state=state, endcond=ec)
    if pitchindeg:
        pitch = np.arccos(-pitch)*180/np.pi
        label = 'Pitch [degree]'
    else:
        label = 'Pitch'
    if phi360:
        phi = np.where(phi < 0, phi.to_value() + 360, phi)
        yrange = [0, 360]
    else:
        yrange = [-180, 180]

    # --- Plots
    if plotRphi:
        if ax is None:
            fig, ax = plt.subplots()
        H, Redges, phiedges, im = ax.hist2d(R, phi, bins=150, cmap=cmap, norm=norm)
        if flag_getval:
            pos = np.argwhere(H==H.max())
            if pos.shape == (1,2): # if there is only one bin with a max value
                R_max = [Redges[pos[0][0]], Redges[pos[0][0]+1]] # center is between x and x+1
                phi_max = [phiedges[pos[0][1]], phiedges[pos[0][1]+1]]
            else:
                raise Exception('More than one maximum. Reduce number of bins or code it!')
        plt.xlabel('R [m]', fontsize=fs)
        plt.ylabel('phi [deg]', fontsize=fs)
        plt.xlim(0.8,1.6)
        plt.ylim(yrange)
        if title is None:
            if ec is None:
                ec='all'
            title = state+' state, endcondition '+ec
        plt.title(title,fontsize=fs)
        plt.tight_layout()

    if plotRxi:
        if ax is None:
            fig, ax = plt.subplots()
        H, Redges, pedges, im = ax.hist2d(R, pitch, bins=50)
        if flag_getval:
            pos = np.argwhere(H==H.max())
            if pos.shape == (1,2): # if there is only one bin with a max value
                Rini = [Redges[pos[0][0]], Redges[pos[0][0]+1]] # center is between x and x+1
                pini = [pedges[pos[0][1]], pedges[pos[0][1]+1]]
            else:
                raise Exception('More than one maximum. Reduce number of bins or code it!')
        plt.xlabel('R [m]', fontsize=fs)
        plt.ylabel(label, fontsize=fs)
        plt.xlim(0.8,1.7)
        plt.ylim(-0.2,-0.9)
        if title is None:
            if ec is None:
                ec='all'
            title = state+' state, endcondition '+ec
        plt.title(title,fontsize=fs)
        plt.tight_layout()

    if plotxi:
        if ax is None:
            fig, ax = plt.subplots()
        n, binedg, c = ax.hist(pitch, bins=150, alpha=0.5, color='tomato')
        bin_center = binedg[:-1]+np.diff(binedg/2)
        ax.plot(bin_center,n,':r')
        ax.set_xlabel(label,fontsize=fs)
        ax.set_ylabel('Counts',fontsize=fs)
        ax.set_title(title,fontsize=fs+2)
        if title is None:
            if ec is None:
                ec='all'
            title = state+' state, endcondition '+ec
        plt.title(title,fontsize=fs)
        plt.tight_layout()

    if flag_getval:
        return Rini, pini


def plot_markers_velspace(file, qid_run: str = 'active', ax = None):
    """
    Plot markers initial and final distribution in the velocity space
    @author:    Javier Hidalgo-Salaverri <jhsalaverri@us.es>

    :param file: path to the .h5 file
    :qid_run: run qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    :param ax: axis where to plot, if None a new one wil be created
    """
    a5 = Ascot(file)
    try:
        int(qid_run) # If it is a number, will not raise an error
        qid_run = 'q'+qid_run
    except ValueError:
        pass
    data = a5.__getattribute__(qid_run)
    energy_beg = data.inistate.get('energy')/1e3
    pitch_beg = np.arccos(data.inistate.get('pitch'))/np.pi*180

    energy_end = data.endstate.get('energy')/1e3
    pitch_end = np.arccos(data.endstate.get('pitch'))/np.pi*180
    if ax is None:
            fig, ax = plt.subplots(2,2)
    elif ax.shape != (2,2):
        raise NameError('ax shape must be (2,2)')
    ax[0,0].hist(energy_beg.value)
    ax[0,0].set_xlabel('Energy [keV]')
    ax[0,0].set_title('Initial distribution')
    ax[1,0].hist(pitch_beg)
    ax[1,0].set_xlabel('Pitch [deg]')

    ax[0,1].hist(energy_end.value)
    ax[0,1].set_xlabel('Energy [keV]')
    ax[0,1].set_title('End distribution')
    ax[1,1].hist(pitch_end)
    ax[1,1].set_xlabel('Pitch [deg]')

    plt.tight_layout()
    plt.show(block = False)



