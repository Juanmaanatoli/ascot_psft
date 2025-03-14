"""
Script to process distributions and plot slowing down distr with custom options.
These routines are not great, and not adapted to the new ASCOT, and probably work mostly for MAST-U. 
But they are a starting point.
In any case, it is probably worth checking the new built-in functions that ASCOT has for working
with distributions before trying to fix this.

@author: L. Velarde - lvelarde@us.es

NOT CHECKED
"""

import numpy as np
import a5py
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt; plt.ion()
from scipy import interpolate
from a5py import Ascot
from marker.alias import get as alias
import ascot_psft.marker.interpret as interpret
import warnings
warnings.filterwarnings('ignore')
import a5py.dist as distmod
from matplotlib.colors import LinearSegmentedColormap
cmapidl = LinearSegmentedColormap.from_list(
    'mycmap', ['black', 'blue', 'red', 'yellow', 'white'], N=256)

fs = 18 #fontsize for all figures and labels
plt.rc('xtick',labelsize=fs)
plt.rc('ytick',labelsize=fs)
# plt.savefig('distrbirth_rz.png',dpi=700, bbox_inches = "tight")   



def get_dist(self, dist=None, **kwargs):
        """
        Return distribution dictionary.

        The ordinate is density which is integrated over the dimensions which
        are given in kwargs.

        Args:
            dist : dict_like, optional <br>
               Give input distribution explicitly instead of reading one from
               HDF5 file.
            kwargs : Name(s) (R, phi, z, ppa, ppe, time, charge) of those
               dimensions along which the distribution is either sliced or
               integrated over. Names are given as keyword value pairs where
               a tuple value (a, b) are indices for slicing and array [a, b]
               are indices for integration. A scalar zero means whole
               dimension is integrated.

        Returns:
            Distribution dictionary.
        """
        if not dist:
            dist = distmod.histogram2distribution(self.read())
        distmod.squeeze(dist, **kwargs)

        return dist


def plot_dist_2D(dist, *args, logscale=False, equal=False, axes=None, cmap=True, log=False):
    """
    Plot distribution as a 2D plot (pcolormesh).

    This function assumes the given distribution is squeezed so that only two
    dimensions remain.

    Args:
        dist : dict_like <br>
            Distribution dictionary.
        args : str, str <br>
            Name of the x and y coordinates e.g. "R", "z".
        axes : Axes, optional <br>
            Axes to which the distribution is plotted. If None, a new figure is
            created and displayed.
        cmap : optional <br>
            # If True, the cmap gnuplot is used
            If True, the cmap from IDL is used
    """

    newfig = axes is None
    if newfig:
        plt.figure()
        axes = plt.gca()

    # dist = h5.active.dist5d.read()
    ordinate = None
    if "distribution" in dist:
        ordinate = dist["distribution"]
    elif "histogram" in dist:
        ordinate = dist["histogram"]

    if len(args) == 0:
        x = dist["abscissae"][0]
        y = dist["abscissae"][1]
    else:
        x = args[0]
        y = args[1]

    if x == dist["abscissae"][0]:
        ordinate = np.transpose(ordinate)

    if logscale:
        ordinate = np.log10(ordinate)

    ordinate = np.ma.masked_invalid(ordinate)

    print(x, y)

    distx = dist[x]
    if x == 'energy':
        distx = distx/1e3 #convert to keV
        label = 'Energy [keV]'
    else:
        label = x

    if log: # we can plot the logarithm of the energy. is this useful though?
        print(dist[x][7:])
        print(np.log(dist[x][7:]))
    
    if cmap:
        if log:
            mesh = axes.pcolormesh(np.log(distx[7:]), dist[y], ordinate[:,7:], cmap=plt.cm.gnuplot  , vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))
        else:
            mesh = axes.pcolormesh(distx, dist[y], ordinate, cmap=cmapidl, vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))
            # mesh = axes.contourf(distx, dist[y], ordinate, cmap=cmapidl, vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))
            # mesh = axes.contour(distx, dist[y], ordinate, levels = 40, cmap=cmapidl, vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))
    else:
        mesh = axes.pcolormesh(distx, dist[y], ordinate, vmin=np.nanmin(ordinate), vmax=np.nanmax(ordinate))

    # https://stackoverflow.com/a/16125413/190597 (Joe Kington)
    # and https://stackoverflow.com/a/35905483
    axes.patch.set(hatch='x', edgecolor=[0.9, 0.9, 0.9])
    cbar = plt.colorbar(mesh, ax=axes)
    cbar.set_label('Counts', fontsize=18)

    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', direction='out')
    axes.set_xlabel(label, fontsize=fs)
    axes.set_ylabel(y, fontsize=fs)

    if equal:
        axes.axis("image")
    else:
        axes.axis("tight")

    if newfig:
        plt.show(block=False)



def plot_E_xi_dist(fn, qid=None, E_edges=np.linspace(-15e3,100e3,100), xi_edges=70,
                       logscale=False, equal=False, axes=None, dist=None):
        """
        Convert (ppa, ppe) to (E, xi) and plot the distribution.

        Either makes a 1D or a 2D plot depending on number of input arguments.

        Args:
            fn: path to file
            args : str, str <br>
                Name of the x-coordinate, and optionally y-coordinate if a 2D
                plot is desired.
            axes : Axes, optional <br>
                Axes where plotting happens. If None, a new figure will be
                created.
            equal : bool, optional <br>
                Make axes equal.
            dist : dict_like, optional <br>
               Give input distribution explicitly instead of reading one from
               HDF5 file. Dimensions that are not x or y are integrated over.
        """
        args = ['energy', 'pitch'] # we define here the arguments, as this is to 
                                    # convert and plot to xi,E
        abscissae = {"r" : 0, "phi" : 0, "z" : 0, "energy" : 0,
                     "pitch" : 0, "time" : 0, "charge" : 0}

        x = alias(args[0])
        del abscissae[x]
        y = None
        if len(args) > 1:
            y = alias(args[1])
            del abscissae[y]

        h5 = Ascot(fn)
        if qid is not None:
            flag_ac = False # Uses spcecified markers.
        else:
            flag_ac = True # Uses active markers.

        print('The active marker is being used?:',flag_ac)
        print("The markers' qid is:",qid)

        if flag_ac:
            dist5d = h5.active.dist5d
        else:
            dist5d = h5['q'+qid].dist5d


        dist = dist5d.get_E_xi_dist(E_edges=E_edges, xi_edges=xi_edges)

        for k in list(abscissae.keys()):
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        plot_dist_2D(dist, x, y, logscale=logscale, equal=equal, axes=axes)


