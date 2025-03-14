"""
Script to manage FILD data from MAST-U, should be easy to adapt to AUG (probably no changes needed).

@author: L. Velarde - lvelarde@us.es


List of functions:
    read_FILD. Routine to read the losses and return the variables for the losses that get to FILD.
    FILD_XIhist_varioussims. Routine to plot in a same plot the XI-histograms of the FILD losses from different simulations.
                             Can also plot experimental data.
    FILD_countlosses. Routine to count the number of markers that reach FILD.
    FILD_get_power. Routine to calculate the power at each triangle in the FILD head.
    FILD_get_hists. Routine to plot a histogram with the losses that reach FILD.
    FILD_comparestate_pitchrange. Routine to compare the initial/end aspects of populations 
        within a certain pitch range
    plot_FILD_FI_dist. Routine to plot the FILD FI distr.
    read_TPboundary. Read pickle with trapped-passing boundary info.
"""

import numpy as np
try:
    import xarray as xr
except ModuleNotFoundError:
    print('Could not import xarray')
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    print('Could not import tqdm')
from a5py import Ascot
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt; plt.ion()
import matplotlib.colors as colors
from ascot_psft.aux_files.common import read_bf, eV2J
import warnings
warnings.filterwarnings('ignore')
import IPython


colorlist = ['black', 'blue', 'red', 'green', 'darkorange','cyan', 'purple', 'grey',
             'darkgreen', 'pink', 'peru', 'magenta', 'lime', 'gold', 'turquoise', 
             'indianred','slategray','brown']


def read_FILD(h5, sim, variables = ['energy','pitch','r','z','ids'], state = 'end', xi_deg = True):
    """
    Read the losses and return the variables for the losses that get to FILD 
    Will also print some statistics about the losses getting to FILD
    Is usually called from other routines.

    Inputs:
        h5. Ascot file. So the bfield can be initialised to read the pitch.
        sim. Simulation read from the ASCOT output.
        variables. The variables that will be read from FILD and returned.
                   CAREFUL! If 'pitch' is in the list, it will be converted to degrees if xi_deg is True.
        state. State of the markers to be read. Default is 'end'.
        xi_deg. Flag to convert the pitch to degrees. If False, will be kept in [0,1] range.
        
    Returns:
        FILD_var. Dictionary with the variables read from the markers that hit FILD.
    """
    # initialise bfield
    bfqid = sim.bfield.get_qid()
    print('Initialising bfield')
    h5.input_init(bfield = bfqid)
    
    # get the tiles where markers hit
    wall_tile = sim.getstate('walltile',mode='prt',state='end',endcond='wall')
    nlost = len(wall_tile)
    # read the wall
    wall = sim.wall.read()
    # get their flags in a boolean type. True will be FILD flags
    # CAREFUL! The wall_tile values are literally the number of the wall tiles
    # Need to rest 1 to use them as index
    FILDflag = wall['flag'][wall_tile-1].astype(bool).reshape(nlost)
    # Read the info of the markers that get to FILD
    FILD_time = sim.getstate('time',mode='prt',state='end',endcond='wall')[FILDflag]
    if len(FILD_time) == 0:
        print('CAREFUL. You are using an old wall without FILD tags.')
        print('Will try to proceed the old way.')
        nwall_total = wall['nelements']
        nFILD = 309
        FILDflag = (wall_tile>nwall_total-nFILD)&(wall_tile<=nwall_total) #particles that hit wall tiles covered by FILD

    # create the dictionary to store the variables
    FILD_var = {}
    for var in variables:
        # if the variable is the pitch, we will convert it to degrees and store it as FILD_Pitch
        # the rest of the variables will be stored in a dictionary (FILD_var) with the name of the variable
        if var == 'pitch':
            FILD_pitch = sim.getstate(var,mode='prt',state=state,endcond='wall')[FILDflag]
            if xi_deg:
                FILD_pitch = np.arccos(-FILD_pitch)*180/np.pi
        else:
            FILD_var[var] = sim.getstate(var,mode='prt',state=state,endcond='wall')[FILDflag]

    if 'pitch' in variables:
        FILD_var['pitch'] = FILD_pitch
    
    # free the bfield
    h5.input_free(bfield=True)
    
    return FILD_var


def FILD_XIhist_varioussims(fnlist: list, qidlist: list, labellist: list = None, xi_deg = True, 
                normalise: bool = False, lineplot: bool = False, exptal_data: str = None, factor = 20, 
                histtype = 'bar', alpha = 0.5, range = [0,90], title = None, filter_R_range: list = None, 
                filter_z_range: list = None, filter_E_range: list = None, filter_time_range: bool = False):
    """
    Routine to plot in a same plot the XI-histograms of the FILD losses from different simulations.
    They can either be in different files or in the same file.
    In principle this only works for endstate, but can easily be modified so it also accets inistate.

    Inputs:
    fn. List of paths to the ascot files.
    qidlist. List of qids (str, can only be the 10 digit ID number). 
    labellist. List of labels for the legend (str). If none, the description of the simulation will be used.
    xi_deg. If True, the pitch will be plotted in degrees.
    normalise. If True, will normalise the histograms.
    lineplot. If True, a contour line of the histogram will be plotted.
    exptal_data. String pointing to a netcdf file with experimental data, to compare the pitch distribution.
    factor. Factor by which the experimental data will be divided for an easier comparison.
    histtype. Type of histogram. Can be 'bar' or 'step'.
    alpha. Transparency of the histogram.
    range. Range of pitch in the histogram.
    title. Title for the plot.
    filter_R_range. If a list is given with a range of R values (e.g. [1.5,1.52]),
        a second histogram on top of the other one with only the particles that end within
        that R range will be plotted. For this, rule of thumb: if FILD @ 1.505, use [1.5,1.51].
        Or if FILD @ 1.49, then [1.485,1.495]
    filter_z_range. If a list is given with a range of z values (e.g. [0.15,0.18]),
        a second histogram on top of the other one with only the particles that end within
        that z range will be plotted.
    filter_E_range. If a list is given with a range of E values in eV (e.g. [50e3,75e3]),
        a second histogram on top of the other one with only the particles that end within
        that E range will be plotted.
    filter_time_range. This is a boolean flag. If True, will first plot a 1D histogram of the 
        times it takes for the particles to reach FILD. Thenk, you will be asked the t_promptlosses to
        only plot the losses that get lost before that time.

    Returns:
        counts. Number of losses reaching FILD for each run.
    """    
    if type(fnlist)==str:
        print('Only one file is given. Will use the same file for all the runs')
        fnlist = [fnlist]*len(qidlist)
    elif len(fnlist) == 1:
        print('Only one file is given. Will use the same file for all the runs')
        fnlist = [fnlist[0]]*len(qidlist)
    
    if type(qidlist)==str:
        print('Only one qid is given. Will convert to list and use the same qid for all the runs')
        qidlist = [qidlist]*len(fnlist)
    elif len(qidlist) == 1:
        print('Only one qid is given. Will convert to list and use the same qid for all the runs')
        qidlist = [qidlist[0]]*len(fnlist)
        
    if len(fnlist) != len(qidlist) and len(qidlist) != 1 and len(fnlist) != 1:
        print('The number of files and qids do not match. Please check.')
        return 
    
    fig, axs = plt.subplots()
    counts = []
    for ii, fn in enumerate(fnlist):
        c = colorlist[ii]
        print('~~~~~~~~~~~~~~~~~~')
        print(fn + ' is being analysed')
        h5 = Ascot(fn)
        qid = qidlist[ii]
        try:
            int(qid) # If it is a number, will not raise an error
            qid = 'q'+qid
        except ValueError:
            pass
        sim = h5.data.__getattribute__(qid)
        print("The run's qid is:",qid)
        if labellist is None:
            lab = sim.get_desc()
        else:
            lab = labellist[ii]
        print('Run '+lab+' is being processed')

        fild_data = read_FILD(h5, sim, variables=['ekin', 'pitch', 'r', 'z', 'time'], state='end', xi_deg=True)
        FILD_energy = fild_data['ekin']
        FILD_pitch = fild_data['pitch']
        FILD_r = fild_data['r']
        FILD_z = fild_data['z']
        FILD_time = fild_data['time']
        
        
        if filter_R_range is not None and filter_z_range is not None and filter_E_range is not None:
            mask_R_z_E = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])&(FILD_z>filter_z_range[0])&(FILD_z<filter_z_range[1])&(FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])
            FILD_pitch = FILD_pitch[mask_R_z_E]
            title += ', R, z and E filters'
        elif filter_R_range is not None and filter_E_range is not None:
            mask_R_E = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])&(FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])
            FILD_pitch = FILD_pitch[mask_R_E]
            title += ', R and E filters'
        elif filter_R_range:
            mask_R = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])
            FILD_pitch = FILD_pitch[mask_R]
            title += ', R filter'
        elif filter_z_range:
            mask_z = (FILD_z>filter_z_range[0])&(FILD_z<filter_z_range[1])
            FILD_pitch = FILD_pitch[mask_z]
            title = ', z filter'
        elif filter_E_range:
            mask_E = (FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])
            FILD_pitch = FILD_pitch[mask_E]
            title = ', E filter'
        elif filter_time_range:
            if ii == 0:
                plt.figure()
                plt.hist(FILD_time.to_value(), bins=50)
                plt.xlabel('Time [s]')
                plt.ylabel('Counts')
                t_promptlosses = float(input('What is the desired t to filter?'))
            mask_time = (FILD_time<t_promptlosses)
            FILD_pitch = FILD_pitch[mask_time]
            title = ', time filter'

        counts.append(len(FILD_pitch))

        # Plot xi distribution
        if lineplot:
            n, binedg = np.histogram(FILD_pitch, density=normalise, bins=150, range=range) 
            bin_center = binedg[:-1]+np.diff(binedg/2) 
            if normalise:
                axs.plot(bin_center,n/n.max(),alpha=alpha, label=lab, color=c)
            else:
                axs.plot(bin_center,n,alpha=alpha, label=lab, color=c)
        else:
            axs.hist(FILD_pitch, density=normalise, histtype = histtype, bins=300, 
                     range=range, alpha=alpha, label=lab, color=c)
        print('-------------------------')
    
    
    if exptal_data:
        print('Reading and plotting experimental data')
        b = xr.open_dataset(exptal_data)
        # t_list = [0.54,0.55,0.56,0.57,0.58,0.59]
        t_list = [0.55,0.58]
        # fig, ax = plt.subplots()
        for t0 in t_list:
            if normalise:
                factor = (b['integral_over_y'].sel(t=t0,method='nearest')).max()
                ((b['integral_over_y'].sel(t=t0,method='nearest'))/factor).plot(ax=axs,label=str(t0),linestyle='--')
            else:
                ((b['integral_over_y'].sel(t=t0,method='nearest'))/factor).plot(ax=axs,label=str(t0),linestyle='--')

    if xi_deg:
        axs.set_xlabel('Pitch [degree]')
    else:
        axs.set_xlabel('Pitch [adim]')
    if normalise:
        axs.set_ylabel('Normalised counts')
        axs.set_ylim([0,1.02])
    else:
        axs.set_ylabel('Counts')
    axs.legend()
    axs.set_title(title)
    plt.tight_layout()

    return counts


def FILD_countlosses(fn, qid: str = 'active'):
    """
    Routine to count the number of markers that reach FILD and provide some statistics:
    - Total lost markers
    - Percentage of lost markers
    - Number of markers that hit FILD
    - Percentage of markers that hit FILD

    Inputs:
    fn. Path to the ascot file.
    qid. Run's qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q. 
        
    Returns:
    countsfild. Number of markers that get to FILD.
    nlost. Number of markers that hit the wall.
    ntot. Total number of markers.
    """
    print('~~~~~~~~~~~~~~~~~~')
    print('File: ',fn)
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)

    fild_data = read_FILD(h5, sim, variables=['ekin'], state='end')
    # Get number of markers that get to FILD
    countsfild = len(fild_data['ekin'])    
    # Get markers with endcondition = wall
    ec = 'wall'
    state = sim.getstate_markersummary()[0] # [0] to keep only endconds
    nlost = [item[0] for item in state if item[1].lower() == ec.lower()][0]
    # Get total number of markers
    ntot = len(sim.getstate('r',mode='prt',state='ini'))
    
    # print some statistics
    print('-------------------------')
    print(str(nlost)+' particles hit the wall')
    print('That is a '+str(nlost/ntot*100)+' %')
    print(str(countsfild)+' particles hit FILD')
    print('That is a '+str(countsfild/ntot*100)+' %'+' of total markers')
    print('That is a '+str(countsfild/nlost*100)+' %'+' of lost markers')

    return countsfild, nlost, ntot


def FILD_get_power(fn, qid: str = 'active', FILD_flag: str = 'FILD', r_filter: float = None):
    """
    Routine to calculate the power at each triangle in the FILD head and return
    the averaged heat flux at the FILD Plasma Facing Surface (PFS).

    Inputs:
    fn. Path to the ascot file.
    qid. Run's qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    FILD_flag. Name for the flag of the desired FILD. 
    r_filter. R position of the FILD PFS. If none, the averaged
        heat flux in the entire FILD head will be returned.
        In MASTU, if FILD@1.455, r_filter = 1.457 --> R_FILD + 0.002 as rule of thumb
    """
    print('~~~~~~~~~~~~~~~~~~')
    print('File: ',fn)
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)
    print("The run's description is:", sim.get_desc())

    # find elements on the wall that correspond to our selected flag
    wall = sim.wall.read()
    try: # for compatibility with old wall format
        flagnumber = wall['flagIdList'][wall['flagIdStrings'].index(FILD_flag)]
    except KeyError:
        flagnumber = wall['labels'][FILD_flag]
    ids_flag = np.where(wall['flag'] == flagnumber)[0] # ids of the FILD triangles
    ids_wettedtiles, area, ploads, markerloads, _ = sim.getwall_loads(weights=True)
    find = np.in1d(ids_wettedtiles, ids_flag) # Wetted tile IDs that are in FILD
    print('Power getting to ' + str(np.sum(find)) + ' triangles in FILD')
    heatflux = ploads / area
    index_wettedtiles = ids_wettedtiles - 1 # convert IDs to indexes
    # use the indexes for the barycenters and then
    # get the FILD barycenters and filter by R or whatever
    # FILD_heatflux[find] = ploads[find]/area[find]
    barycenters = sim.wall.barycenters(cartesian=False)[index_wettedtiles]
    bary_r_FILD = barycenters[find][:,0] # get the r of the FILD wetted triangles
    if r_filter:
        flag_PFS = bary_r_FILD <= r_filter
        heatflux_PFS = heatflux[find][flag_PFS]
    else:
        heatflux_PFS = heatflux[find]

    return np.mean(heatflux_PFS)
    

def FILD_get_hists(fn, qid: str = 'active',state='end', flag_return: bool = False, xi_deg = True, title=None,
                plot_evsxi: bool = False, plot_xi: bool = False, plot_rz = False, 
                plot_r = False, write_txt: str = None, filter_R_range: list = None,
                filter_z_range: list = None, filter_E_range: list = None, 
                apply_filter_R = False, apply_filter_z = False,
                apply_filter_E = False, apply_filter_both = False,
                exptal_data: str = None, factor = 20):
    """
    Routine to plot a histogram with the losses that reach FILD.

    Inputs:
    path. Path to the ascot file.
    qid. Run's qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    flag_return. To return the FILD parameters calculated.
    xi_deg. If True, the pitch will be plotted in degrees.
    Title. Title for the plots.
    plot_evsxi. Flag to plot the E vs XI histogram.
    plot_xi. Flag to plot the XI histogram.
    plot_rz. Flag to plot the rz histogram.
    plot_r. Flag to plot the r histogram.
    write_txt. Flag to write a txt with the pitch values for comparison with FILDSIM.
        If None, it won't write it. If a fname is given, it will.
        The name should follow: 'pitchASCOT_shot_time_filter.txt'
        For example, for both filters: 'pitchASCOT_47132_300ms_both.txt'
    filter_R_range. If a list is given with a range of R values (e.g. [1.5,1.52]),
    a second histogram on top of the other one with only the particles that end within
    that R range will be plotted. For this, rule of thumb: if FILD @ 1.505, use [1.5,1.51].
    Or if FILD @ 1.49, then [1.485,1.495]
    filter_z_range. If a list is given with a range of z values (e.g. [0.15,0.18]),
    a second histogram on top of the other one with only the particles that end within
    that z range will be plotted.
    filter_E_range. If a list is given with a range of E values in eV (e.g. [50e3,75e3]),
    a second histogram on top of the other one with only the particles that end within
    that E range will be plotted.
    apply_filter_R. Flag to only write to the txt the R filtered values of the pitch.
    apply_filter_E. Flag to only write to the txt the E filtered values of the pitch.
    apply_filter_both. Flag to only write to the txt the R and E filtered values of the pitch.
    exptal_data. String pointing to a netcdf file with experimental data, to compare the pitch distribution.
    factor. Factor by which the experimental data will be divided for an easier comparison.
    """
    print('~~~~~~~~~~~~~~~~~~')
    print('File: ',fn)
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)

    # read fild data
    fild_data = read_FILD(h5, sim, variables=['ekin', 'pitch', 'r', 'z', 'phi', 'ids'], state=state, xi_deg=xi_deg)
    FILD_energy = fild_data['ekin']
    FILD_Pitch = fild_data['pitch'] # CAREFUL, depending on the xi_deg flag this will be in degrees (True) or [0,1] range (False)
    FILD_r = fild_data['r']
    FILD_z = fild_data['z']
    FILD_phi = fild_data['phi']
    FILD_ids = fild_data['ids']


    if filter_R_range is not None and filter_z_range is not None and filter_E_range is not None:
        mask_R_z_E = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])&(FILD_z>filter_z_range[0])&(FILD_z<filter_z_range[1])&(FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])
    elif filter_R_range is not None and filter_E_range is not None:
        mask_R_E = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])&(FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])
    if filter_R_range:
        mask_R = (FILD_r>filter_R_range[0])&(FILD_r<filter_R_range[1])
    if filter_z_range:
        mask_z = (FILD_z>filter_z_range[0])&(FILD_z<filter_z_range[1])
    if filter_E_range:
        mask_E = (FILD_energy>filter_E_range[0])&(FILD_energy<filter_E_range[1])

    if plot_evsxi:
        plt.figure()
        h1 = plt.hist2d(FILD_Pitch, FILD_energy/1000, bins=100, range=[[0,90],[0,95]],cmap=plt.cm.inferno)
        cbar = plt.colorbar(h1[3])
        cbar.set_label('Counts',labelpad=10)
        plt.xlabel('Pitch [deg]')
        plt.ylabel('Energy [keV]')
        plt.tight_layout()
        plt.show()

    if plot_xi:
        plt.figure()
        n, binedg, c = plt.hist(FILD_Pitch, bins=150, range=[0,90], alpha=0.5, color='tomato', label='No filter')
        if filter_R_range:
            n, binedg, c = plt.hist(FILD_Pitch[mask_R], bins=150, range=[0,90], alpha=0.5, color='blue', label='R filter')
        if filter_z_range:
            n, binedg, c = plt.hist(FILD_Pitch[mask_z], bins=150, range=[0,90], alpha=0.5, color='blue', label='z filter')
        if filter_E_range:
            n, binedg, c = plt.hist(FILD_Pitch[mask_E], bins=150, range=[0,90], alpha=0.5, color='green', label='E filter')
        if filter_R_range is not None and filter_z_range is not None and filter_E_range is not None:
            n, binedg, c = plt.hist(FILD_Pitch[mask_R_z_E], bins=150, range=[0,90], alpha=0.5, color='yellow', label='R, z and E filters')
        elif filter_R_range is not None and filter_E_range is not None:
            n, binedg, c = plt.hist(FILD_Pitch[mask_R_E], bins=150, range=[0,90], alpha=0.5, color='yellow', label='R and E filters')
        # bin_center = binedg[:-1]+np.diff(binedg/2)
        # plt.plot(bin_center,n,':b')
        if exptal_data:
            b = xr.open_dataset(exptal_data)
            t_list = [0.3]
            # fig, ax = plt.subplots()
            for t0 in t_list:
                ((b['integral_over_y'].sel(t=t0,method='nearest'))/factor).plot(label=str(t0))

        plt.xlabel('Pitch [degree]')
        plt.ylabel('Counts')
        plt.title(title)
        plt.legend(fontsize=fs-2)
        plt.show()
        plt.tight_layout()

    if plot_r:
        plt.figure()
        n, binedg, c = plt.hist(FILD_r, bins=150, range=[0,2], alpha=0.5, color='tomato', label='No filter')
        if filter_R_range:
            n, binedg, c = plt.hist(FILD_r[mask_R], bins=150, range=[0,2], alpha=0.5, color='blue', label='R filter')
        if filter_E_range:
            n, binedg, c = plt.hist(FILD_r[mask_E], bins=150, range=[0,2], alpha=0.5, color='green', label='E filter')
        if filter_R_range is not None and filter_E_range is not None:
            n, binedg, c = plt.hist(FILD_r[mask_R_E], bins=150, range=[0,2], alpha=0.5, color='yellow', label='R and E filters')
        # bin_center = binedg[:-1]+np.diff(binedg/2)
        # plt.plot(bin_center,n,':b')
        plt.xlabel('R [m]')
        plt.ylabel('Counts')
        plt.legend()
        # plt.title(title)
        plt.show()
        plt.tight_layout()

    
    if write_txt is not None:
        if apply_filter_R:
            np.savetxt(write_txt, FILD_Pitch[mask_R])
        if apply_filter_E:
            np.savetxt(write_txt, FILD_Pitch[mask_E])
        if apply_filter_both:
            np.savetxt(write_txt, FILD_Pitch[mask_R_E])
        else:
            np.savetxt(write_txt, FILD_Pitch)
        

    if plot_rz:
        plt.figure()
        h1 = plt.hist2d(FILD_r, FILD_z, bins=1000, norm=colors.LogNorm(), range=[[0,2],[-2,2]],cmap=plt.cm.inferno)
        cbar = plt.colorbar(h1[3])
        cbar.set_label('Counts',labelpad=10)
        plt.xlabel('R [m]')
        plt.ylabel('z [m]')
        plt.show()
        plt.tight_layout()
        ax = plt.gca()
        ax.set_facecolor('k')


    if flag_return:
        return FILD_time,FILD_energy,FILD_pitch,FILD_Pitch,FILD_r,FILD_z, FILD_ids


def FILD_comparestate_inoutpitchrange(fn, qid: str = 'active', pitchrange = [60,70], state = 'ini', 
                                pitchstate = 'end', xi_deg = True, density=True):
    """
    Routine to compare the initial aspects of populations within a certain pitch range. 
    (Can also be used for the end aspects).
    This will show histograms of the initial/final R,z,phi,E, comparing that for particles
    inside the "pitchrange" and outside of it.

    Inputs:
    fn. Path to the ascot file.
    qid. Run's qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    pitchrange. Range of pitch for histogram comparison. In degrees by default (xi_deg).
    state. Initial or final state of the markers. 'ini' or 'end'.
    pitchstate. Whether the "pitchrange" specified is from the ini or endstate.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)

    FILD_data = read_FILD(h5, sim, variables=['pitch','ids'], state=pitchstate, xi_deg=xi_deg)
    FILD_Pitch = FILD_data['pitch']
    FILD_ids = FILD_data['ids']
    
    # Get the indices of the markers that satisfy the condition
    indpitchin = np.where((FILD_Pitch>=pitchrange[0]) & (FILD_Pitch<=pitchrange[1]))[0]
    indpitchout = np.where((FILD_Pitch<pitchrange[0]) | (FILD_Pitch>pitchrange[1]))[0]
    # careful, IDs ARE NOT RELATED TO INDEXES!!!

    # Read the variables
    variables = sim.getstate('ids','r','z','phimod','ekin',mode='prt',state=state)
    idsvariables = variables[0]
    # keep those inside the range
    idsFILDin = FILD_ids[indpitchin] # ID of markers that satisfy the condition
    mask = np.isin(idsvariables, idsFILDin)
    modified_var_in = [sublist[mask] for sublist in variables]
    idsin, riin, ziin, phiiin, energyiin = modified_var_in
    # those outside the range
    idsFILDout = FILD_ids[indpitchout] # ID of markers that satisfy the condition
    mask = np.isin(idsvariables, idsFILDout)
    modified_var_out = [sublist[mask] for sublist in variables]
    idsout, riout, ziout, phiiout, energyiout = modified_var_out

    fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4)
    ax1.hist(riin.to_value(),bins=50,density=density,histtype='step',color='b',label='Pitch within '+str(pitchrange)+' deg')
    ax1.hist(riout.to_value(),bins=50,density=density,histtype='step',color='orange',label='Pitch outside '+str(pitchrange)+' deg')
    ax1.set_xlabel('R [m]', fontsize=fs)
    ax1.set_ylabel('Counts', fontsize=fs)
    ax1.legend(loc='upper left')
    
    ax2.hist(ziin.to_value(),bins=40,histtype='step',density=density,color='b',label='Pitch within '+str(pitchrange)+' deg')
    ax2.hist(ziout.to_value(),bins=40,histtype='step',density=density,color='orange',label='Pitch outside '+str(pitchrange)+' deg')
    ax2.set_xlabel('z [m]', fontsize=fs)
    ax2.set_ylabel('Counts', fontsize=fs)
    ax2.legend(loc='upper left')
    
    ax3.hist(phiiin.to_value(),bins=40,range=[0,360],density=density,histtype='step',color='b',label='Pitch within '+str(pitchrange)+' deg')
    ax3.hist(phiiout.to_value(),bins=40,range=[0,360],density=density,histtype='step',color='orange',label='Pitch outside '+str(pitchrange)+' deg')
    ax3.set_xlabel('phi [deg]', fontsize=fs)
    ax3.set_ylabel('Counts', fontsize=fs)
    ax3.legend(loc='upper left')

    # ax4.hist(energyiin.to_value()/1e3,bins=50,density=density,histtype='step',color='b',label='Pitch < '+str(philim)+' deg')
    # ax4.hist(energyiout.to_value()/1e3,bins=50,density=density,histtype='step',color='orange',label='Pitch > '+str(philim)+' deg')
    # ax4.set_xlabel('Energy [keV]', fontsize=fs)
    # ax4.set_ylabel('Counts', fontsize=fs)
    # ax4.legend(loc='upper left')

    ax4.hist(energyiin.to_value()/1e3,bins=50,histtype='step',color='b',label='Pitch within '+str(pitchrange)+' deg')
    ax4.hist(energyiout.to_value()/1e3,bins=50,histtype='step',color='orange',label='Pitch outside '+str(pitchrange)+' deg')
    ax4.set_xlabel('Energy [keV]', fontsize=fs)
    ax4.set_ylabel('Counts', fontsize=fs)
    ax4.legend(loc='upper left')


def FILD_comparestate_pitchranges(fn, qid: str = 'active', pitchranges = [[60, 70]], state='ini', 
                                pitchstate='end', xi_deg=True, density=True):
    """
    Routine to compare the initial aspects of populations within specified pitch ranges.
    
    Inputs:
    - fn: Path to the ascot file.
    - qid: Run's qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    - pitchranges: List of pitch ranges for histogram comparison. In degrees by default (xi_deg).
    - state: Initial or final state of the markers. 'ini' or 'end'.
    - pitchstate: Whether the "pitchranges" specified are from the ini or endstate.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    sim = h5.data.__getattribute__(qid)
    
    # Read pitch and IDs for FILD 
    FILD_data = read_FILD(h5, sim, variables=['pitch', 'ids'], state=pitchstate, xi_deg=xi_deg)
    FILD_Pitch = FILD_data['pitch']
    FILD_ids = FILD_data['ids']

    # Read the initial/final state variables
    variables = sim.getstate('ids', 'r', 'z', 'phimod', 'ekin', mode='prt', state=state)
    idsvariables = variables[0]  # IDs

    # Plot histograms for each pitch range
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15, 6))

    for idx, pitchrange in enumerate(pitchranges):
        # Indices of markers within and outside the pitch range
        indpitchin = np.where((FILD_Pitch >= pitchrange[0]) & (FILD_Pitch <= pitchrange[1]))[0]
        idsFILDin = FILD_ids[indpitchin]
        mask_in = np.isin(idsvariables, idsFILDin)

        # Variables within this pitch range
        r_in = variables[1][mask_in]
        z_in = variables[2][mask_in]
        phimod_in = variables[3][mask_in]
        ekin_in = variables[4][mask_in]

        # Plotting for each variable
        color = colorlist[idx % len(colorlist)]  # Cycle through colors

        # R histogram
        ax1.hist(r_in.to_value(), bins=50, range=[1.3, 1.5], density=density, histtype='step', 
                 color=color, label=f'Pitch within {pitchrange} deg')
        ax1.set_xlabel('R [m]', fontsize=12)
        ax1.set_ylabel('Counts', fontsize=12)

        # z histogram
        ax2.hist(z_in.to_value(), bins=40, histtype='step', density=density, color=color, 
                 label=f'Pitch within {pitchrange} deg')
        ax2.set_xlabel('z [m]', fontsize=12)
        ax2.set_ylabel('Counts', fontsize=12)

        # phi histogram
        ax3.hist(phimod_in.to_value(), bins=40, range=[0, 360], density=density, histtype='step', 
                 color=color, label=f'Pitch within {pitchrange} deg')
        ax3.set_xlabel('phi [deg]', fontsize=12)
        ax3.set_ylabel('Counts', fontsize=12)

        # Energy histogram
        ax4.hist(ekin_in.to_value() / 1e3, bins=50, histtype='step', color=color, 
                 label=f'Pitch within {pitchrange} deg')
        ax4.set_xlabel('Energy [keV]', fontsize=12)
        ax4.set_ylabel('Counts', fontsize=12)

    # Adding legends and adjusting layout
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper left')
    ax3.legend(loc='upper left')
    ax4.legend(loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_FILD_FI_dist(file, qid_run: str = 'active', qid_wall: str = 'active',
    FILD_id: str = 'FILD1', ax = None):
    """
    Plot the initial and final distribution of the markers that got to FILD
    @author: Javier Hidalgo-Salaverri <jhsalaverri@us.es>
    NOT CHECKED

    :param file: path to the .h5 file
    :param qid_run: run qid. Can be 'active', 
        the description (in the correct format) or the 10 digit ID number, without the Q.
    :param ax: axis where to plot, if None a new one wil be created
    :param wall_qid: wall qid
    """
    a5 = Ascot(file)
    try:
        int(qid_run) # If it is a number, will not raise an error
        qid_run = 'q'+qid_run
    except ValueError:
        pass
    run_data = a5.__getattribute__(qid_run)

    # Get the affected RFILD index in the wall
    try:
        int(qid_wall) # If it is a number, will not raise an error
        qid_wall = 'q'+qid_wall
    except ValueError:
        pass
    wall = a5.wall.__getattribute__(qid_wall).read()
    id_flag = [flagString == FILD_id for flagString in wall['flagIdStrings']]
    id_flag = int(wall['flagIdList'][id_flag])

    id_marker = []
    for i,m in enumerate(run_data.endstate.read()['walltile']):
        if wall['flag'][m-1] == id_flag and m > 0:
            id_marker.append(i)

    if len(id_marker) > 0:
        # Get the energy and pitch of the desired markers
        energy_beg = run_data.inistate.get('energy', ids = id_marker)/1e3
        pitch_beg = np.arccos(run_data.inistate.get('pitch',
                                                    ids = id_marker))/np.pi*180
        energy_end = run_data.endstate.get('energy', ids = id_marker)/1e3
        pitch_end = np.arccos(run_data.endstate.get('pitch',
                                                    ids = id_marker))/np.pi*180
        print(str(len(id_marker)) +' markers got to FILD')
        
        if ax is None:
            fig, ax = plt.subplots()
        # Get a histogram
        fig, ax = plt.subplots(2,2)
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

        fig.suptitle('FI over FILD')
        plt.tight_layout()
        plt.show(block = False)
    else:
        print('No markers got to FILD')


def read_TPboundary(path="46944_tpbdry.npy"):
    """
    Read pickle with trapped-passing boundary info.

    Inputs:
        path. Path to the file with the dictionary.

    """
    test=np.load(path,allow_pickle=True).item()
#    TP bdry at theta=1000
    plt.figure(); plt.contourf(test['pitch'],test['energy'],test['theta'],levels=500)
    plt.colorbar()
    theta0 = np.array([1000])
    plt.contour(test['pitch'],test['energy'],test['theta'],levels=theta0,colors='k')
    plt.xlabel('Pitch [degree]')
    plt.ylabel('Energy [eV]')


def old_FILD_get_power(fn, qid: str = None, FILD_flag = 1, fname_vtk: str = None):
    """
    Routine to calculate the power at each triangle in the FILD head.

    Inputs:
    fn. Path to the ascot file.
    qid. Run's qid. If not specified, the active run will be used.
    FILD_flag. Number for the flag of the desired FILD. Default is 1.
    fname_vtk. Name for the resulting vtk file. If no name is given, it will not be saved.
    """
    print('~~~~~~~~~~~~~~~~~~')
    print(fn)
    h5 = Ascot(fn).data
    if qid is not None:
        flag_ac = False # Uses spcecified markers.
    else:
        flag_ac = True # Uses active markers.
    print('The active run is being used?:',flag_ac)
    print("The markers' qid is:",qid)

    if flag_ac:
        sim = h5.active
    else:
        sim = h5["run_"+qid]

    # FILD_wall_tile: in which triangle each of the markers that reach FILD hit it
    FILD_energy, FILD_weight, FILD_wall_tile, FILD_time = \
        read_FILD(h5, sim, variables=['energy','weight','walltile', 'time'],state='end')
    FILD_EJ = eV2J(FILD_energy)
    power = FILD_EJ*FILD_weight # power of each marker
    power_pl = FILD_EJ[FILD_time<=0.1e-3]*FILD_weight[FILD_time<=0.1e-3] 
    power_del = FILD_EJ[FILD_time>0.1e-3]*FILD_weight[FILD_time>0.1e-3]
    
    # calculate total elements and FILD elements
    nwall_total = sim.wall.read()['nelements']
    flag = (sim.wall.read()['flag']==FILD_flag).reshape(nwall_total)
    print('Number being used as FILD flag: ', FILD_flag)
    nFILD = np.sum(flag)
    IDstriangs = np.arange(nwall_total) + 1
    triangs_FILD = IDstriangs[flag]

    # We know the power of each marker (power) and the triangle it hits (FILD_wall_tile)
    # We can calculate the power at each triangle with the markers that hit it
    # Looking only at the FILD triangles (triangs_FILD)
    power_triangles = np.zeros(len(triangs_FILD))
    print('Calculating the power at each triangle...')
    for ii, tri in tqdm(enumerate(triangs_FILD)):
        mask = FILD_wall_tile == tri
        power_triangles[ii] = np.sum(power[mask])
    # now we need to normalise by the area of the triangle
    # we need to get the coordinates of the FILD triangles
    xcoord = sim.wall.read()['x1x2x3'][flag]
    ycoord = sim.wall.read()['y1y2y3'][flag]
    zcoord = sim.wall.read()['z1z2z3'][flag]
    # with the x,z,y coordinates of each of the triangle's point, calculate the area
    # be aware that these are 3D coordinates, so we will use the dot product of two vectors
    areas = np.zeros(nFILD)
    print('Calculating the area of each triangle...')
    for ii in tqdm(range(nFILD)):
        vec1 = np.array([xcoord[ii][0]-xcoord[ii][1],ycoord[ii][0]-ycoord[ii][1],zcoord[ii][0]-zcoord[ii][1]])
        vec2 = np.array([xcoord[ii][0]-xcoord[ii][2],ycoord[ii][0]-ycoord[ii][2],zcoord[ii][0]-zcoord[ii][2]])
        areas[ii] = 0.5*np.linalg.norm(np.cross(vec1,vec2))

    # normalise the power by the area
    flux_triangles = power_triangles/areas

    print('Plotting heat flux histogram')
    plt.hist(flux_triangles,bins=50)
    plt.xlabel('flux [MW/m2]')
    plt.ylabel('Counts')
    plt.tight_layout()

    if fname_vtk is not None:
        vtk = sim.wall.toVtk()
        # to create a scalar field with the heat flux at each of the FILD's triangles
        # we need to create a vector with the length of the entire wall, and 
        # fill it with zeros, and then fill the FILD triangles with the heat flux
        flux_all = np.zeros(nwall_total)
        flux_all[triangs_FILD-1] = flux_triangles
        vtk.addScalarField('FILD heat flux',flux_all)
        vtk.writeVtp(fname_vtk)
