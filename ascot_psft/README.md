# ascot_psft
Libraries of utilities to work with the ASCOT code.

To use the repository, you will likely need to add the folder where you cloned it to the python path, so you can import the "ascot_psft" module. This can be achieved by doing:
sys.path.append(`path`) in ipython before starting to work with the repository. In Marenostrum for example, if you cloned it in the home, `path` will be something like '/home/us/username/'.

Alternatively, you can use the export PYTHONPATH method outside of ipython:

export PYTHONPATH='`path`'

If you then do $PYTHONPATH, you should see the path to your home directory. 

This option is in general better as you can add that "export" line to you bash file so it will be exported automatically every time you log into MareNostrum.

## cmd files
Some example of cmd files have been added to the cmdfiles folder. They all assume that you created a soft link called:

"bbnbi5" if running BBNBI.

"ascot5_mpi" if running ASCOT from the folder where it was compiled with MPI=1. 

"ascot5_noMPI" if running ASCOT from the folder where it was compiled without MPI=1.

(The soft link should point to the corresponding ascot5 folder in each case)


## MyRoutines folder
The MyRoutines folder is a place where you can keep your routines that you do not wish to upload to the repository. The changes in this folder are not tracked by git.


## Use in MareNostrum
If you encounter any errors of the type "ModuleNotFoundError", please raise an issue here as we need to ask the IT service to install the modules for us.
Currently, the "tqdm", "xarray" and the "netCDF4" modules have been installed. 
The MAST-U vessel module needs to be copied from Marconi so it cannot be used yet.

For interactive plotting, use the "-Y" flag when doing the ssh, and if you encounter problems
when trying to close the figures using the close button, try closing them from the terminal instead.
