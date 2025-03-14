"""
Perform simulation in batches using the same .h5 file in MARCONI

DISCLAIMER: this script does not work if the total simulations time is over
10 minutes as the Marconi interactive maximum time is reached

by Javier Hidalgo Salaverri (jhsalaverri@us.es)

"""

import numpy as np
from a5py.ascot5io.ascot5 import Ascot
import subprocess

def bbnbi_batch(file_id, bbnbi5, plasma_id, bfield_id, nbi_id, descriptions,
    nmarkers: int = 1e6):
    """
    Run a series of bbnbi5 simulations

    @params file_id: path to the .h5 file
    @params bbnbi5: path to the bbnbi5.cmd file
    @params plasma_id: array with the ids
    @params bfield_id: array with the ids. Same dimensions as plasma_id
    @params nbi_id: array with the ids. Same dimensions as plasma_id
    @params descriptions: array with the descriptions of the produced markers
    @params nmarkers: number of markers for the simulation
    """
    a5 = Ascot(file_id)
    nsim = len(plasma_id)
    cmd = 'sbatch ' + bbnbi5 + ' ' + file_id + ' ' + str(int(nmarkers))
    for i in range(nsim):
        # Set the inputs
        a5.reload()
        a5.bfield.__getattribute__(bfield_id[i]).set_as_active()
        a5.plasma.__getattribute__(plasma_id[i]).set_as_active()
        a5.nbi.__getattribute__(nbi_id[i]).set_as_active()
        a5.reload()

        # Run simulation
        subprocess.call(cmd, shell = True)
        # sim_time = code_utilities.time_in_seconds()
        code_utilities.wait4sim(return_time = False)

        # Get the run name
        a5.reload()
        qid = code_utilities.get_last_run_id(file_id, 'marker')
        # time_since_last_run = sim_time-\
        #     code_utilities.date2sec(a5.marker.__getattribute__(qid).get_date())
        # if time_since_last_run > 10.0:
            # raise NameError('Last created marker is too old (> 10 seconds)')

        # Change the descriptions
        a5.marker.__getattribute__(qid).set_desc(descriptions[i])

def ascot_batch(file_id, ascot5, plasma_id, bfield_id, nbi_id, marker_id,
    descriptions):
    """
    Run a series of ascot5 simulations

    @params file_id: path to the .h5 file
    @params ascot5: path to the ascot5.cmd file
    @params plasma_id: array with the ids
    @params bfield_id: array with the ids. Same dimensions as plasma_id
    @params nbi_id: array with the ids. Same dimensions as plasma_id
    @params marker_id: array with the ids. Same dimensions as plasma_id
    @params descriptions: array with the descriptions of the produced markers

    @# TODO: Check if nbi is needed
    """
    a5 = Ascot(file_id)
    nsim = len(plasma_id)
    cmd = 'sbatch ' + ascot5 + ' ' + file_id
    for i in range(nsim):
        # Set the inputs
        a5.reload()
        a5.bfield.__getattribute__(bfield_id[i]).set_as_active()
        a5.plasma.__getattribute__(plasma_id[i]).set_as_active()
        a5.nbi.__getattribute__(nbi_id[i]).set_as_active()
        a5.marker.__getattribute__(marker_id[i]).set_as_active()
        a5.reload()

        # Run simulation
        subprocess.call(cmd, shell = True)
        # sim_time = code_utilities.time_in_seconds()
        code_utilities.wait4sim(return_time = False)

        # Get the run name
        a5.reload()
        qid = code_utilities.get_last_run_id(file_id)
        # time_since_last_run = sim_time-\
        #     code_utilities.date2sec(a5.__getattribute__(qid).get_date())
        # if time_since_last_run > 10.0:
        #     raise NameError('Last created run is too old (> 10 seconds)')

        # Change the descriptions
        a5.__getattribute__(qid).set_desc(descriptions[i])
