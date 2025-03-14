# 0.1.1 Adapted to other machines
- Fixed paths
- Commented the tkAgg as it does not work on tok
- Added routine for vessel plotting
- Added cmd files for various machines

# 0.1.0 MareNostrum ready
- Modified some imports for interactive plotting in Marenostrum
- Finally, will keep the check_bfield script with all bf plots
- Renamed plots_ascot to plots_outputs 

# 0.0.10 New MyRoutines + opened to other users
- Included MyRoutines folder for personal scripts
- Included a general "manage_wall" file with move_wall and wall_to_stl routines
- Updated all "devices" variables to either "MU" or "AUG"
- Included "plotSDD" script to work with distributions & slowing down distr. Needs fixing.

## 0.0.9 New FILD routines (previous commits)
- Scripts to analyse losses in FILD
- Scripts to look at wall loads

## 0.0.8 Transport and FOM
- Scripts to analyse DeltaP_phi
- Script to calculate some figures of merit

## 0.0.7 Included MAST-U scripts
- Routines for input generation MAST-U specific
- Routines for plotting markers from the simulation
- markers_grid routine for resonances studies
- Routines to check the bf and the plasma inputs

## 0.0.6 QOL improvements
- Added a routine to check if markers are inside the vessel
- Improved wall_hits

## 0.0.5 Added 2D AUG wall
- Added 2D AUG wall to vessel_mesh.h5
- Fixed a bug on write_inputs in the wall section

## 0.0.4 Update to be coherent with new ASCOT functions
- Updated files so they can work with the new ASCOT version (5.5.3)
- Remove the batch functionalities as they were not working properly
- Added functionalities to add the NBI and insert FILD (no longer support natively by ASCOT)
- insert_FILD now allows to add different flags for FILD and the vessel for easy postprocessing
- Added a master vessel_mesh as repository. Added the AUG wall in 3D

## 0.0.3 Naming coherence and ax support
- Kept name coherence along the code
- Added support for external axes in the plotting routines

## 0.0.2 Added __init__ file and extra utilities
- Libraries can be added by import ascot_psft
- Plot utilities have been added

## 0.0.1 Creation of the git and first scripts
- Added scripts to plot inputs, write inputs, and perform simulation in a batch
