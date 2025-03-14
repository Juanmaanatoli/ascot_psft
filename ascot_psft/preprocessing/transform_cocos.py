"""
Transform cocos of a file.
This should be copied to MyRoutines and modified there.
"""
from a5py.physlib import cocos as cocosmod
from freeqdsk import geqdsk

cocos_in = 5
cocos_out = 15
fn_in = '/marconi_work/FUA38_ASCOT_3D/lvelarde/mastu/fild/MU03/48302/MAST_48302_550ms.geqdsk'
phiclockwise = False # phi clockwise from above
weberperrad = True # psi normalised to 2pi (divided)
fn_out = '/marconi_work/FUA38_ASCOT_3D/lvelarde/mastu/fild/MU03/48302/MAST_48302_550ms_COCOS15.geqdsk'

# read file
with open(fn_in, "r") as f:
    eqd = geqdsk.read(f)
# Check the COCOS_in is correct
cocos_check = cocosmod.assign(
    eqd["qpsi"][0], eqd["cpasma"], eqd["bcentr"], eqd["simagx"],
    eqd["sibdry"], phiclockwise, weberperrad)
if cocos_check != cocos_in:
    raise ValueError('The input COCOS is not correct (or the other flags)')
eqdout = cocosmod.fromCocosNtoCocosM(eqd, cocos_n = cocos_in, cocos_m = cocos_out)
# write it to new file
with open(fn_out, "w") as f:
    geqdsk.write(eqdout, f)
