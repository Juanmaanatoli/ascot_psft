"""
Script to check the generated plasma input and modify it

@author: L. Velarde - lvelarde@us.es

"""

from a5py import Ascot
import numpy as np
import copy
try:
    import netCDF4 as nc
except ModuleNotFoundError:
    print('Could not import netCDF4')
import matplotlib
# # matplotlib.use('TkAgg')
import matplotlib.pyplot as plt; plt.ion()
from a5py.ascot5io import plasma_1D
from scipy.interpolate import interp1d

fs = 18 #fontsize for all figures and labels
plt.rc('xtick', labelsize=fs)
plt.rc('ytick', labelsize=fs)

colorlist = ['black', 'blue', 'red', 'green', 'darkorange','cyan', 'purple', 'grey',
             'darkgreen', 'pink', 'peru', 'magenta', 'lime', 'gold', 'turquoise', 
             'indianred','slategray','brown']


def plot_profs(fn:str, qid: str = 'active'):
    """
    Plot plasma input.

    Inputs:
        fn. Path to the ascot file.
        qid: ID of the plasma to plot. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(6,7),sharex=True)
    ti = plasma['itemperature']
    te = plasma['etemperature']
    ne = plasma['edensity']
    ni = plasma['idensity']
    rho = plasma['rho']
    ax1.plot(rho,ti,label=desc)
    ax1.set_ylabel('Ti [eV]',fontsize=fs)
    ax2.plot(rho,te,label=desc)
    ax2.set_ylabel('Te [eV]',fontsize=fs)
    ax3.plot(rho,ne,label=desc)
    ax3.set_ylabel('ne [m⁻3]',fontsize=fs)
    ax4.plot(rho,ni,label=desc)
    ax4.set_ylabel('ni [m⁻3]',fontsize=fs)
    plt.xlabel('rho',fontsize=fs)
    plt.tight_layout()
    plt.legend()


def extend_plasma(fn:str, qid:str = 'active', rhoend = 1.45, nrho = 20,
                neend = 1e8, teend = 5, tiend = 5):
    """
    Read a previous plasma input and extend it radially.

    Inputs:
        fn. Path to the ascot file.
        qid: ID of the plasma that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        rhoend. Last value of the exponential fit.
        nrho. Number of tho points that will be included.
        neend, teend, tiend. Values of ne, Te, Ti, that will be extended.

    Returns:
        new. New plasma input.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
    
    te = plasma['etemperature']
    ti = plasma['itemperature']
    ne = plasma['edensity']
    nrho = plasma['nrho']
    nion = plasma['nion']
    anum = plasma['anum'][0]
    znum = plasma['znum'][0]
    mass = plasma['mass'][0]
    charge = plasma['charge'][0]
    rho = plasma['rho']
    desc = desc + '_rhoend'+str(rhoend)

    lastrho = rho[-1]
    print('Rho ends at: ', str(lastrho))
    print('Extending it to: ', str(rhoend))
    rho_new = np.linspace(lastrho,rhoend,nrho)
    ne_new = neend * np.ones(nrho-1)
    te_new = teend * np.ones(nrho-1)
    ti_new = tiend * np.ones(nrho-1)
    # and add the new values
    ne = np.append(ne, ne_new)
    te = np.append(te, te_new)
    ti = np.append(ti, ti_new)
    rho = np.append(rho, rho_new[1:])

    nrho = len(rho)
    idensity = ne.reshape(int(nrho),int(nion))
    new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
        charge, rho, ne, te, idensity, ti, desc=desc)
    h5 = Ascot(fn).data
    h5.plasma[new].activate()
    print('New plasma input created and set as active.')
    print('Plotting it now.')
    plot_profs(fn)

    return new


def mod_plasma_SOL_expdecay(fn:str, qid:str = 'active', rhoini = 1.1, rhoend = 1.6, cut = 1.45,
                            neend = 1e8, teend = 1, tiend = 1):
    """
    Read a previous plasma input and modify the ne, Te and Ti for rho>rhoini by using
    an exponential decay function. This can be used to "realistically" extend the plasma 
    input when ascot complains.

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        version. Version of the modified input. Int value!
        rhoini. Initial value at which the exponential will begin.
        rhoend. Last value of the exponential fit.
        cut. Rho value at which the profile will be cut.
        neend / teend / tiend. Value of ne/te/ti at rhoend.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
    
    te = plasma['etemperature']
    ti = plasma['itemperature']
    ne = plasma['edensity']
    nrho = plasma['nrho']
    nion = plasma['nion']
    # nion = plasma['nion']
    anum = plasma['anum'][0]
    znum = plasma['znum'][0]
    mass = plasma['mass'][0]
    charge = plasma['charge'][0]
    rho = plasma['rho']
    # desc = 'ModSOLfrom_q'+qid+'_expdecay_v'+str(version)
    desc = 'ModSOLfrom_q'+qid+'_expdecay_rhoini'+str(rhoini)+'rhoend'+str(rhoend)

    print('\nPlotting the modified profile')
    # Plot before modifying    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(6,7),sharex=True)
    ni = plasma['idensity']
    ax1.plot(rho,ti,'k',label=desc_orig)
    ax2.plot(rho,te,'k',label=desc_orig)
    ax3.plot(rho,ne,'k',label=desc_orig)
    ax4.plot(rho,ni,'k',label=desc_orig)

    def exp_decay(x, yBeg, yEnd,):
        b = np.log(yEnd/yBeg)/(x[0]-x[-1])
        a = yEnd/np.exp(-b*x[-1])
        y_out = a*np.exp(-b*x)
        return y_out

    indrho = np.argmin(np.abs(rho-rhoini))
    rho_new = np.linspace(rho[indrho],rhoend,50)
    ne_new = exp_decay(rho_new, ne[indrho], neend)
    te_new = exp_decay(rho_new, te[indrho], teend)
    ti_new = exp_decay(rho_new, ti[indrho], tiend)
    
    # cut the original rho, ne, te and ti arrays until rho = 1
    ne = ne[rho<rho[indrho]]
    te = te[rho<rho[indrho]]
    ti = ti[rho<rho[indrho]]
    rho = rho[rho<rho[indrho]]
    # and add the new values
    ne = np.append(ne, ne_new)
    te = np.append(te, te_new)
    ti = np.append(ti, ti_new)
    rho = np.append(rho, rho_new)
    # finally, cut them again at rho = cut
    if rhoend > cut:
        ne = ne[rho<cut]
        te = te[rho<cut]
        ti = ti[rho<cut]
        rho = rho[rho<cut]
    
    nrho = len(rho)
    idensity = ne.reshape(int(nrho),int(nion))
    # Plot after modifying
    ax1.plot(rho,ti,'b--',label=desc)
    ax1.set_ylabel('Ti [eV]',fontsize=fs)
    ax2.plot(rho,te,'b--',label=desc)
    ax2.set_ylabel('Te [eV]',fontsize=fs)
    ax3.plot(rho,ne,'b--',label=desc)
    ax3.set_ylabel('ne [m⁻3]',fontsize=fs)
    ax4.plot(rho,idensity,'b--',label=desc)
    ax4.set_ylabel('ni [m⁻3]',fontsize=fs)
    plt.xlabel('rho',fontsize=fs)
    plt.tight_layout()
    plt.legend()

    print('The description for this profile would be ', desc)
    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, ne, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        qid_new = new[-10:]
        print('New plasma input created and set as active.')
        print('Plotting it now to compare against the original one.')
        compare_profs(fn, [qid,qid_new])
    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def scan_ne_edge(fn:str, qid:str = 'active', rhoini = 0.9, factor = 2,
                neend = 1e8, rhoend = 1.1):
    """
    Read a previous plasma input and modify the density at the edge.

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        rhoini. First value of the ne that will be modified.
        factor. factor by which the ne will be 
        neend. Value at which to fix ne when it goes below 0.
        rhoend. rho from which the ne will be brought to 'zero' (neend).
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
        
    print('\nThe current plasma profile input is ', desc)
    print('Plotting current plasma profiles')
    plot_profs(fn, qid)

    ne = plasma['edensity']
    rho = plasma['rho']
    if factor > 1:
        desc = desc + '_needgex'+str(factor)
    else:
        desc = desc + '_needgex'+str(factor)+'_rhoend'+str(rhoend)

    indrho = np.argmin(np.abs(rho-rhoini))
    # indrhoend = np.argmin(np.abs(rho-rhoend))
    # mask = (rho>rho[indrho]) & (rho<rho[indrhoend])
    mask = (rho>rho[indrho])
    # keep only the edge ne
    ne_edge = ne[mask][:-1]
    m = np.diff(ne_edge)
    new_m = m * factor
    new_array = np.zeros_like(ne_edge)
    new_array[0] = ne_edge[0]  # Keep the first value same as original

    for i in range(1, len(ne_edge)):
        new_array[i] = new_array[i-1] + new_m[i-1]

    
    # Copy the original ne and replace the new values
    ne_new = np.copy(ne)
    idx = np.where(mask)[0]
    idx_rpl = idx[:-1] # we will leave out the last value
    ne_new[idx_rpl] = new_array.reshape(len(new_array),1)

    # check for negative values and substitute
    ne_new[ne_new < 0] = neend
    # bring to "zero" for large rhos
    ne_new[rho>rhoend] = neend

    print('\nPlotting the modified profile')
    plt.figure()
    plt.plot(rho,ne,'k', label='Original profile')
    plt.plot(rho,ne_new, 'b--', label='New profile')
    plt.ylabel('ne [m⁻3]',fontsize=fs)
    plt.xlabel('rho',fontsize=fs)
    plt.tight_layout()
    plt.legend()

    print('The description for this profile would be ', desc)
    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        # get the rest of the plasma data
        te = plasma['etemperature']
        ti = plasma['itemperature']
        nion = plasma['nion']
        anum = plasma['anum'][0]
        znum = plasma['znum'][0]
        mass = plasma['mass'][0]
        charge = plasma['charge'][0]
        nrho = len(rho)
        idensity = ne_new.reshape(int(nrho),int(nion))
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, ne_new, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        print('New plasma input created and set as active.')
        print('Plotting it now.')
        plot_profs(fn)

    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def move_ne_radially(fn:str, qid:str = 'active', rhotras = 0.02, Rtras = None):
    """
    Read a previous plasma input and move it radially.

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        rhotras. The amount of traslation in rho.
        Rtras. The amount of traslation in R. CAREFUL. Only one of them can be defined.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
        
    print('\nThe current plasma profile input is ', desc)
    print('Plotting current plasma profiles')
    plot_profs(fn, qid)

    ne = plasma['edensity']
    rho = plasma['rho']
    if Rtras is None:
        desc = desc + '_nemoverho' + str(rhotras)
    elif rhotras is None:
        desc = desc + '_nemoveR' + str(Rtras)
        print('Not implemented yet. Will exit.')
        return
    else:
        print('Only rho or R can be used. Not both.')
        return

    # Shift rho values by the desired
    shifted_rho = rho - rhotras

    # Extend the initial ne value to fill the range starting from 0
    if rhotras < 0:
        # Filter out ne values where shifted_rho < 0
        valid_indices = shifted_rho <= rho[-1]
        shifted_rho = shifted_rho[valid_indices]
        shifted_ne = ne[valid_indices]
        extended_rho = np.insert(shifted_rho, 0, 0)
        extended_ne = np.insert(shifted_ne, 0, ne[0])
    else:
        # Filter out ne values where shifted_rho < 0
        valid_indices = shifted_rho >= 0
        shifted_rho = shifted_rho[valid_indices]
        shifted_ne = ne[valid_indices]

        # Extend the last ne value to fill the range to the last rho
        rhoextension = np.linspace(shifted_rho[-1], rho[-1], 4).reshape(4)
        neextension = np.linspace(shifted_ne[-1], ne[-1], 4).reshape(4)
        extended_rho = np.concatenate((shifted_rho,rhoextension))
        extended_ne = np.concatenate((shifted_ne,neextension))
    
    # Interpolate to old rho grid
    # interp_func = interp1d(extended_rho, extended_ne)
    interp_func = interp1d(extended_rho, extended_ne, bounds_error=False, fill_value="extrapolate")
    interpolated_ne = interp_func(rho)

    # Plot the results
    print('\nPlotting the modified profile')
    plt.figure()
    plt.plot(rho, ne,'bo', label='Original', markersize=0.5)
    plt.plot(shifted_rho, shifted_ne, 'k', label='Shifted')
    plt.plot(extended_rho, extended_ne,'g--', label='Extended')
    plt.plot(rho, interpolated_ne,'r-.', label='Interpolated')
    plt.xlabel('rho',fontsize=fs)
    plt.ylabel('ne',fontsize=fs)
    plt.title(('rhotras = ' + str(rhotras)))
    plt.legend()
    plt.tight_layout()

    print('The description for this profile would be ', desc)
    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        # get the rest of the plasma data
        te = plasma['etemperature']
        ti = plasma['itemperature']
        nion = plasma['nion']
        anum = plasma['anum'][0]
        znum = plasma['znum'][0]
        mass = plasma['mass'][0]
        charge = plasma['charge'][0]
        nrho = len(rho)
        idensity = interpolated_ne.reshape(int(nrho),int(nion))
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, interpolated_ne, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        print('New plasma input created and set as active.')
        print('Plotting it now.')
        plot_profs(fn)

    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def multiply_ne(fn:str, qid:str = 'active', factor = None):
    """
    Read a previous plasma input and multiply it.

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        factor. Factor by which the ne will be multiplied.
            If None, will ask as input, once the ne(rho=0) has been printed.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
    
    print('\nThe current plasma profile input is ', desc)
    print('Plotting current plasma profiles')
    plot_profs(fn, qid)

    ne = plasma['edensity']
    rho = plasma['rho']
    if factor is None:
        print('ne at the core is: ', ne[0][0])
        indrho = np.argmin(np.abs(rho-1))
        print('ne at the separatrix is: ', ne[indrho][0])
        dum = input('Do you know the factor you want? [y/n] ')
        if dum == 'y':
            factor = input('Please specify the factor: ')
        elif dum == 'n':
            dum2 = input('Do you know the target value at the core? [y/n] ')
            if dum2 == 'y':
                ne_core = float(input('Please specify the target ne_core (without the 1e19): '))
                factor = ne_core * 1e19 / ne[0][0]
            else:
                print('Then think about it before running this again :)')
                return
    else:
        print('The density profile will be multiplied by: ', str(factor))

    desc = desc + '_netotx' + str(factor)
    ne_new = ne * factor
    
    # Plot the results
    print('\nPlotting the modified profile')
    plt.figure()
    plt.plot(rho, ne,'bo', label='Original', markersize=0.5)
    plt.plot(rho, ne_new, 'k', label='Multiplied by ' + str(factor))
    plt.xlabel('rho',fontsize=fs)
    plt.ylabel('ne',fontsize=fs)
    plt.legend()
    plt.tight_layout()

    print('The description for this profile would be ', desc)
    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        # get the rest of the plasma data
        te = plasma['etemperature']
        ti = plasma['itemperature']
        nion = plasma['nion']
        anum = plasma['anum'][0]
        znum = plasma['znum'][0]
        mass = plasma['mass'][0]
        charge = plasma['charge'][0]
        nrho = len(rho)
        idensity = ne_new.reshape(int(nrho),int(nion))
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, ne_new, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        print('New plasma input created and set as active.')
        print('Plotting it now.')
        plot_profs(fn)

    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def flatten_ne_core(fn:str, qid:str = 'active'):
    """
    Read a previous plasma input and multiply it.

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
        
    print('\nThe current plasma profile input is ', desc)
    print('Plotting current plasma profiles')
    plot_profs(fn, qid)

    ne = plasma['edensity']
    ne_new = ne.copy()
    rho = plasma['rho']

    indrho = np.argmin(np.abs(rho-0.95))
    print('ne at rho=0.95 is: ', ne_new[indrho][0])
    dum = input('Do you know the value you want? [y/n] ')
    if dum == 'y':
        value = float(input('Please specify the factor: '))
    else:
        print('Then think about it before running this again :)')
        return

    desc = desc + '_neflat' + str(value)
    neind = np.argmin(np.abs(ne_new-value))    
    ne_new[rho<rho[neind]] = ne_new[neind]

    # Plot the results
    print('\nPlotting the modified profile')
    plt.figure()
    plt.plot(rho, ne,'b', label='Original', markersize=0.5)
    plt.plot(rho, ne_new, 'k', label='Flattened to ' + str(value))
    plt.xlabel('rho',fontsize=fs)
    plt.ylabel('ne',fontsize=fs)
    plt.legend()
    plt.tight_layout()

    print('The description for this profile would be ', desc)
    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        # get the rest of the plasma data
        te = plasma['etemperature']
        ti = plasma['itemperature']
        nion = plasma['nion']
        anum = plasma['anum'][0]
        znum = plasma['znum'][0]
        mass = plasma['mass'][0]
        charge = plasma['charge'][0]
        nrho = len(rho)
        idensity = ne_new.reshape(int(nrho),int(nion))
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, ne_new, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        print('New plasma input created and set as active.')
        print('Plotting it now.')
        plot_profs(fn)

    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def increase_kineticvar(fn:str, qid:str = 'active', var_names = ['edensity'], addvalues: list = None, 
                        vals_out = [1e8], fnomf = None):
    """
    Read a previous plasma input and increase/decrease the desired variables in it,
    the amount indicated by addvalues.
    Can be compared agains the profiles from the OMFIT SLICE file (only TS data).

    Inputs:
        fn. Path to the ascot file.
        qid. ID of the plasma input that will be modified. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        vars. List of variables that will be modified:
            'edensity' for electron density
            'etemperature' for electron temperature
            'itemperature' for ion temperature
        addvalue. List of amounts that will be added to each kinetic profile (float).
            If None, will ask as input, once the var(rho=0) has been printed.
        vars_out. Values for the data that goes negative.
        fnomf. Path to the OMFIT TS file.
    """
    h5 = Ascot(fn)
    try:
        int(qid) # If it is a number, will not raise an error
        qid = 'q'+qid
    except ValueError:
        pass
    plasma = h5.data.plasma.__getattribute__(qid).read()
    desc = h5.data.plasma.__getattribute__(qid).get_desc()
        
    print('\nThe current plasma profile input is ', desc)
    print('Plotting current plasma profiles')
    plot_profs(fn, qid)

    rho = plasma['rho']
    variables = {}
    # Read the variables and plot before modifying   
    for ii, key in enumerate(var_names):
        variables[key] = plasma[key]

    # Copy for later plotting
    desc_orig = desc
    variables_orig = copy.deepcopy(variables)
    # Calculate the value that will be added if not specified
    if addvalues is None:
        addvalues = [None] * len(var_names)
        # Go trhough the variables and calculate their addvalue
        for ii, key in enumerate(var_names):
            print(key+' at the core is: '+ str(variables[key][0][0]))
            indrho = np.argmin(np.abs(rho-1))
            print(key+' at the separatrix is: '+ str(variables[key][indrho][0]))
            dum = input('Do you know the value you want to add? [y/n] ')
            if dum == 'y':
                if key == 'edensity':
                    addvalues[ii] = float(input('Please specify the value (without the 1e19): '))
                    addvalues[ii] *= 1e19
                else:
                    addvalues[ii] = float(input('Please specify the value: '))
            elif dum == 'n':
                dum2 = input('Do you know the target value at the core? [y/n] ')
                if dum2 == 'y':
                    if key == 'edensity':
                        core_val = float(input('Please specify the target ne_core (without the 1e19): '))
                        addvalues[ii] = core_val * 1e19 - variables[key][0][0]
                    else:
                        core_val = float(input('Please specify the target core value: '))
                        addvalues[ii] = core_val - variables[key][0][0]
                else:
                    print('Then think about it before running this again :)')
                    # return
    
    fig, axes = plt.subplots(len(var_names),1,figsize=(6,7),sharex=True)
    if len(var_names) == 1:
        axes = [axes]
    # plot the errorbars from the omfit file if included
    if fnomf:
        ds = nc.Dataset(fnomf) 
        rhoom = np.sqrt(ds['psi_n'][:].data)
        Teom = ds['T_e'][:].data.reshape(len(rhoom))
        neom = ds['n_e'][:].data.reshape(len(rhoom))
        neom_err = ds['n_e__uncertainty'][:].data.reshape(len(rhoom))
        Teom_err = ds['T_e__uncertainty'][:].data.reshape(len(rhoom))

        for ii, key in enumerate(var_names):
            if key == 'edensity':
                var_om = neom
                var_om_err = neom_err
            axes[ii].errorbar(rhoom, neom, yerr=var_om_err, color='r')
    # Add the corresponding value and correct negative values, then plot
    for ii, key in enumerate(var_names):
        print('The '+ key + ' profile will be added: ' + str(addvalues[ii]))
        variables[key] += addvalues[ii]
        variables[key][variables[key] < 0] = vals_out[ii]
        axes[ii].plot(rho,variables_orig[key],'ko',label=desc_orig)
        axes[ii].plot(rho,variables[key],'b',label='Added '+str(addvalues[ii]))
        axes[ii].set_ylabel(key,fontsize=fs)
        desc = desc + key + str(addvalues[ii])
    plt.xlabel('rho',fontsize=fs)
    plt.legend()
    plt.tight_layout()
    plt.show()


    ans = input('\nDo you wish to write it as a new input? [y/n] ')

    if ans == 'y':
        print('The description for this profile would be ', desc)
        dum = input('Do you wish to specify a different description? [y/n] ')
        if dum == 'y':
            desc = input('Please write it here: ')
        # get the rest of the plasma data
        nion = plasma['nion']
        anum = plasma['anum'][0]
        znum = plasma['znum'][0]
        mass = plasma['mass'][0]
        charge = plasma['charge'][0]
        nrho = len(rho)
        if 'edensity' not in var_names:
            ne = plasma['edensity']
        else:
            ne = variables['edensity']
        if 'etemperature' not in var_names:
            te = plasma['etemperature']
        else:
            te = variables['etemperature']
        if 'itemperature' not in var_names:
            ti = plasma['itemperature']
        else:
            ti = variables['itemperature']
        idensity = ne.reshape(int(nrho),int(nion))
        new = plasma_1D.write_hdf5(fn, nrho, nion, anum, znum, mass,
            charge, rho, ne, te, idensity, ti, desc=desc)
        h5 = Ascot(fn).data
        h5.plasma[new].activate()
        print('New plasma input created and set as active.')
        print('Plotting it now.')
        plot_profs(fn)

    elif ans == 'n':
        print('\nOkay. Maybe try again with new settings :)')


def compare_profs(fn:str, qidlist, labellist = None, fnomf = None, factor=1):
    """
    Plot together different plasma inputs.
    Can also plot the raw data from the TS with the error bars.

    Inputs:
        fn. Path to the ascot file.
        qidlist. IDs of the plasma inputs that will be plotted. List of str with the 10 digit ID number.
        labelist. List of labels for the legend.
        fnomf. Path to omfit file. If given, it will also plot omfit raw data.
        factor. Factor by which the experimental (omfit) ne will be multiplied for the plot.
    """
    h5 = Ascot(fn).data
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(6,7),sharex=True)
    if fnomf:
        ds = nc.Dataset(fnomf) 
        rhoom = np.sqrt(ds['psi_n'][:].data)
        Teom = ds['T_e'][:].data.reshape(len(rhoom))
        neom = ds['n_e'][:].data.reshape(len(rhoom))*factor
        neom_err = ds['n_e__uncertainty'][:].data.reshape(len(rhoom))*factor
        Teom_err = ds['T_e__uncertainty'][:].data.reshape(len(rhoom))

        ax2.errorbar(rhoom, Teom, yerr=Teom_err, color='gray', ls='none', 
                    fmt='o', capsize=3, markersize=2, alpha=0.5)
        ax3.errorbar(rhoom, neom, yerr=neom_err, color='gray', ls='none', 
                    fmt='o', capsize=3, markersize=2, alpha=0.5)

    for ii,qid in enumerate(qidlist):
        desc = h5.plasma['q'+qid].get_desc()
        if labellist:
            label = labellist[ii]
        else:
            label = desc
        plasma = h5.plasma['q'+qid].read()
        ti = plasma['itemperature']
        te = plasma['etemperature']
        ne = plasma['edensity']
        ni = plasma['idensity']
        rho = plasma['rho']
        ax1.plot(rho,ti,label=label,color=colorlist[ii])
        ax1.set_ylabel('Ti [eV]',fontsize=fs)
        ax2.plot(rho,te,label=label,color=colorlist[ii])
        ax2.set_ylabel('Te [eV]',fontsize=fs)
        ax3.plot(rho,ne,label=label,color=colorlist[ii])
        ax3.set_ylabel('ne [m⁻3]',fontsize=fs)
        ax4.plot(rho,ni,label=label,color=colorlist[ii])
        ax4.set_ylabel('ni [m⁻3]',fontsize=fs)
        plt.xlabel('rho',fontsize=fs)
        plt.tight_layout()
        plt.legend()    # plot the errorbars from the omfit file if included


def compare_neprof(fn:str, qidlist, labellist = None, fnomf = None, factor=1):
    """
    Plot together different density profiles.
    Can also plot the raw data from the TS with the error bars.

    Inputs:
        fn. Path to the ascot file.
        qidlist. IDs of the plasma inputs that will be plotted. List of str with the 10 digit ID number.
        labelist. List of labels for the legend.
        fnomf. Path to omfit file. If given, it will also plot omfit raw data.
        factor. Factor by which the experimental (omfit) ne will be multiplied for the plot.
    """
    h5 = Ascot(fn).data
    
    fig, ax = plt.subplots(1,1,figsize=(6,3))
    if fnomf:# plot the errorbars from the omfit file if included
        ds = nc.Dataset(fnomf) 
        rhoom = np.sqrt(ds['psi_n'][:].data)
        neom = ds['n_e'][:].data.reshape(len(rhoom))*factor
        neom_err = ds['n_e__uncertainty'][:].data.reshape(len(rhoom))*factor

        ax.errorbar(rhoom, neom, yerr=neom_err, color='gray', ls='none', 
                    fmt='o', capsize=3, markersize=2, alpha=0.5)
    for ii,qid in enumerate(qidlist):
        desc = h5.plasma['q'+qid].get_desc()
        if labellist:
            label = labellist[ii]
        else:
            label = desc
        plasma = h5.plasma['q'+qid].read()
        ne = plasma['edensity']
        rho = plasma['rho']
        ax.plot(rho,ne,label=label,color=colorlist[ii])
        ax.set_ylabel('ne [m⁻3]',fontsize=fs)
        plt.xlabel('rho',fontsize=fs)
        plt.tight_layout()
        plt.xlim([0.75,1.1])
        plt.legend(loc='upper right') 


def check_change_rho(fn: str, R1: float, bfqid: str = None, z = 0, phi = 0):
    """
    Check the change in rho that corresponds to the TS radial error (1cm at the LFS in MAST-U,
    this has been harcoded)
    Will assume t=0

    Inputs:
    fn. Path to ascot file.
    R1. Radial coordinate of the separatrix
    bfqid. If the bfield is not the active one, the qid is required.
    z. Vertical coordinate at which rho will be evaluated.
    phi. Toroidal coordinate at which rho will be evaluated.
    """
    h5 = Ascot(fn)
    if bfqid:
        h5.input_init(bfield=bfqid)
    else:
        h5.input_init(bfield=True)
    R2 = R1 + 0.01 # data provided by Rory Scanell
    R3 = R1 - 0.01 # data provided by Rory Scanell
    rho1 = h5.input_eval(R1,phi,z,0,'rho')
    rho2 = h5.input_eval(R2,phi,z,0,'rho')
    rho3 = h5.input_eval(R3,phi,z,0,'rho')
    print('rho1 is: ', rho1)
    print('rho2 is: ', rho2)
    print('rho3 is: ', rho3)
    print('The error (+) would be: ', rho2-rho1)
    print('The error (-) would be: ', rho1-rho3)
    if bfqid:
        h5.input_free(bfield=bfqid)
    else:
        h5.input_free(bfield=True)
    

def plot_omfit_profs(fnomf: str, factorne=1, fnascot = None, 
                     qid = None, play = None):
    """
    Plot the profiles directly from the OMFIT files. 
    Include errorbars if possible 

    Inputs:
        fnomf. Path to the OMFIT file.
        type. Type of data:
            'fits' if file where fits are stored, 
            'TS' if SLICE file from TS data, 
            'CXRS' if SLICE file from CXRS data.
        factorne. If the experimental data of the ne needs to be multiplied by some value.
        fnascot. Path to ASCOT file. If included, will also compare to ASCOT input profile.
        qid. ID of the ASCOT plasma input to plot and play with. Can be 'active', 
            the description (in the correct format) or the 10 digit ID number, without the Q.
        play. To play with the fit.
            If play == 'radial', the fit will be shifted radially.
            If play == 'incrase', the fit will be shifted upwards/downwards.
    """
    ds = nc.Dataset(fnomf) 
    rho = np.sqrt(ds['psi_n'][:].data)
    Te = np.zeros(len(rho)); Ti = np.zeros(len(rho)); ne = np.zeros(len(rho))
    Te_err = np.zeros(len(rho)); Ti_err = np.zeros(len(rho)); ne_err = np.zeros(len(rho))

    try:
        Te = ds['T_e'][:].data.reshape(len(rho))
        ne = ds['n_e'][:].data.reshape(len(rho))
        Ti = ds['T_12C6'][:].data.reshape(len(rho))
    except IndexError:
        pass
    try:
        ne_err = ds['n_e__uncertainty'][:].data.reshape(len(rho))
        Te_err = ds['T_e__uncertainty'][:].data.reshape(len(rho))
        Ti_err = ds['T_12C6__uncertainty'][:].data.reshape(len(rho))
    except IndexError:
        pass

    fig, axes = plt.subplots(3,1,figsize=(8,6))

    axes[2].errorbar(rho, ne*factorne, yerr=ne_err*factorne, color='k')
    axes[1].errorbar(rho, Te, yerr=Te_err, color='k')
    axes[0].errorbar(rho, Ti, yerr=Ti_err, color='k')
    axes[2].set_ylabel('ne [m-3]', fontsize=fs)
    axes[1].set_ylabel('Te [eV]', fontsize=fs)
    axes[0].set_ylabel('Ti [eV]', fontsize=fs)
    axes[-1].set_xlabel('rho', fontsize=fs)

    if fnascot:
        h5 = Ascot(fnascot)
        if qid is None:
            raise ValueError('Please specify the plasma input ID')
        try:
            int(qid) # If it is a number, will not raise an error
            qid = 'q'+qid
        except ValueError:
            pass
        plasma = h5.data.plasma.__getattribute__(qid).read()
        arho = plasma['rho']
        ane = plasma['edensity']
        ate = plasma['etemperature']
        ati = plasma['itemperature']
        axes[2].plot(arho, ane, 'r--')
        axes[1].plot(arho, ate, 'r--')
        axes[0].plot(arho, ati, 'r--')
        plt.tight_layout()

        if play == 'increase':
            ans = input('Do you know how much you want to modify the ne? [y/n] ')
            if ans == 'y':
                print('Please input the value you want to add to the ne (without the 1e19)')
                neadd = float(input('Write a 0 if you do not want to modify the ne: '))
            elif ans == 'n':
                ans2 = input('Do you know a max/min value and the corresponding fit value? [y/n] ')
                if ans2 == 'y':
                    nemaxmin = float(input('Please indicate the max/min (without the 1e19): '))
                    nefit = float(input('Please indicate the corresponding fit value (without the 1e19): '))
                    neadd = (nemaxmin - nefit) 
                else:
                    print('Sorry, cannot help you :(')
                    return
            ans = input('Do you know how much you want to modify the Te? [y/n] ')
            if ans == 'y':
                print('Please input the value you want to add to the Te')
                Teadd = float(input('Write a 0 if you do not want to modify the Te: '))
            elif ans == 'n':
                ans2 = input('Do you know a max/min value and the corresponding fit value? [y/n] ')
                if ans2 == 'y':
                    Temaxmin = float(input('Please indicate the max/min (without the 1e19): '))
                    Tefit = float(input('Please indicate the corresponding fit value (without the 1e19): '))
                    Teadd = Temaxmin - Tefit
                else:
                    print('Sorry, cannot help you :(')
                    return

            print(str(neadd) + '*1e9 will be added to the ne prof')
            print(str(Teadd) + ' will be added to the Te prof')
            ane_new = ane + neadd * 1e19
            ate_new = ate + Teadd
            # Plot it
            axes[2].plot(arho, ane_new, 'darkorange')
            axes[1].plot(arho, ate_new, 'darkorange')            


        elif play == 'radial':
            rhotras = float(input('How much do you want to shift rho? (Can be negative): '))
            shifted_rho = arho - rhotras
            # Plot it
            axes[2].plot(shifted_rho, ane, 'darkorange')
            axes[1].plot(shifted_rho, ate, 'darkorange')


        elif play == 'gradient':
            print('Be aware that this only works for the ne')
            rhoini = float(input('From what rho do you want to change the gradient? '))
            indrho = np.argmin(np.abs(arho-rhoini))
            mask = (arho>arho[indrho])
            print('How much do you want to change the gradient?')
            factor = float(input('(Positive will multiply, negative divide): '))
            # keep only the edge ne
            ne_edge = ane[mask][:-1]
            m = np.diff(ne_edge)
            if factor < 0:
                factor = -1 / factor
            new_m = m * factor
            new_array = np.zeros_like(ne_edge)
            new_array[0] = ne_edge[0]  # Keep the first value same as original
            for i in range(1, len(ne_edge)):
                new_array[i] = new_array[i-1] + new_m[i-1]
            # Copy the original ne and replace the new values
            ne_new = np.copy(ane)
            idx = np.where(mask)[0]
            idx_rpl = idx[:-1] # we will leave out the last value
            ne_new[idx_rpl] = new_array.reshape(len(new_array),1)
            negne = float(input('What should the negative values be equal too?: '))
            ne_new[ne_new < 0] = negne
            # Plot it
            axes[2].plot(arho, ne_new, 'darkorange')