"""
To calculate some figures of merits from the ASCOT file.

@author: Lina Velarde <lvelarde@us.es>
"""

import numpy as np
from a5py import Ascot
from ascot_psft.postprocessing import wall_hits as wh



def get_shinethrough(path, qid: str = None):
    """
    Routine to calculate the shine-through from the bbnbi run.

    Inputs:
    path. Path to the ascot file.
    qid. BBNBI run's qid. Required.
    """
    h5 = Ascot(path)
    print("The run's qid is:",qid)
    sim = h5.data["bbnbi_"+qid]
    print("The run's description is:",sim.get_desc())
    
    bbstate = sim.getstate_markersummary()[0] # [0] to keep only endconds
    wall = [item[0] for item in bbstate if item[1] == 'WALL'][0]
    ion = [item[0] for item in bbstate if item[1] == 'IONIZED'][0]

    st = (1 - ion / (ion + wall)) * 100
    print('Shine-through is: ' + str(st))

    return st
    

def get_ionised(path, qid: str = None):
    """
    Routine to get the number of ionised markers from the bbnbi run.

    Inputs:
    path. Path to the ascot file.
    qid. BBNBI run's qid. Required.
    """
    h5 = Ascot(path)
    print("The run's qid is:",qid)
    sim = h5.data["bbnbi_"+qid]
    print("The run's description is:",sim.get_desc())
    
    bbstate = sim.getstate_markersummary()[0] # [0] to keep only endconds
    ion = [item[0] for item in bbstate if item[1] == 'IONIZED'][0]
    print('Number of ionised markers: ' + str(ion))

    return ion


def get_number_markers_ec(path, qidbbnbi: str, qid: str = None, 
                          ec = 'wall', flag_fild = None):
    """
    Routine to get a summary of the markers with a certain endcondition.

    Inputs:
    path. Path to the ascot file.
    qidbbnbi. BBNBI run's qid. To look for total markers. 
    qid. Run's qid. If not specified, the active run will be used.
    ec. End condition to return. Default is wall, but it can be "none" for aborted
        particles, or tmax for confined ones. Also polmax or tormax.
    flag_fild. If specified, will calculate the percentage of lost markers that
               get to FILD. Typical MU flag would be just "FILD", for AUG "FILD1".

    Returns:
        count. Number of markers with our endcond
        perc. Percentage (count/ionised markers)
        percFILD. Percentage of markers that reach FILD (lossesFILD/counts)
    """
	
    h5 = Ascot(path)
    if qid is not None:
        flag_ac = False # Uses spcecified run.
    else:
        flag_ac = True # Uses active run.
    print('The active run is being used?:',flag_ac)
    print("The run's qid is:",qid)
    print("The BBNBI run's qid is:",qidbbnbi)

    bbnbirun = h5.data["bbnbi_"+qidbbnbi]
    if flag_ac:
        sim = h5.data.active
    else:
        sim = h5.data["run_"+qid]

    print("The run's description is:",sim.get_desc())
    print("The endcondtion is:", ec)
    # Get total ionized markers
    bbstate = bbnbirun.getstate_markersummary()[0] # [0] to keep only endconds
    ion = [item[0] for item in bbstate if item[1] == 'IONIZED'][0]
    print('Number of ionised markers: ' + str(ion))
    # Get markers with our endcondition
    state = sim.getstate_markersummary()[0] # [0] to keep only endconds
    count = [item[0] for item in state if item[1].lower() == ec.lower()][0]
    print('Number of ' + ec + ' markers: ' + str(count))
    perc = count / ion * 100
    print('That is the ' + str(perc)+ ' %')

    if flag_fild:
        if ec.lower() != 'wall':
            raise ValueError('WALL must be used as endcondition')

        idsFILD, lossesFILD = wh.get_FILD_load(path,flag=flag_fild,run_id=qid)
        # print('Number of markers that reach FILD: ' + str(lossesFILD))
        percFILD = lossesFILD / count * 100
        print('That is the ' + str(percFILD)+ ' %')


    return count, perc, percFILD


def calc_perc_change(pathlist, qidbbnbilist, qidlist, flag_fild = 'FILD'):
    """
    Calculate the change in losses from one scenario to another.
    It will be pathlist[1] - pathlist[0]

    Inputs:
        pathlist. List of the 2 paths to be compared.
        qidbbnbilist. List of the 2 bbnbi's qids to be compared.
        qidlist. List of the 2 run's qids to be compared.
    """
    old_losses, old_perc, old_percFILD = \
    get_number_markers_ec(pathlist[0], qidbbnbilist[0], qidlist[0], flag_fild=flag_fild)
    new_losses, new_perc, new_percFILD = \
    get_number_markers_ec(pathlist[1], qidbbnbilist[1], qidlist[1], flag_fild=flag_fild)
    
    change_losses = ((new_perc - old_perc) / old_perc) * 100
    change_lossesFILD = ((new_percFILD - old_percFILD) / old_percFILD) * 100
    print('The change in global losses is: ' + str(change_losses) + ' %')
    print('The change in FILD losses is: ' + str(change_lossesFILD) + ' %')
