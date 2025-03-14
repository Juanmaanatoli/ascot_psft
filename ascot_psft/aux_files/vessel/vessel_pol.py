"""
Plots the AUG or MAST-U poloidal view of the vessel on top of input "ax" or in
active figure if none

Author: J.F. Rivero-Rodriguez (University of Seville)

File: vessel.py
"""

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

import scipy.io as spio

paths = [
        '/home/us/us776713/ascot_psft/aux_files/vessel',
        '/lustre/home/lvelarde/ascot_psft/aux_files/vessel',
        '/tokp/work/pcano/lina/ascot_psft/aux_files/vessel'
        ]

def AUG(ax=None,color='k',linewidth=0.7):

    for path in paths:
        if os.path.exists(path):
            print('Vessel read from: ', path)
            filepath = [path+'/vesselAUG.mat']
            
    # filepath = [s+"/vesselAUG.mat" for s in sys.path if os.path.isfile(s+"/vesselAUG.mat")]
    # if len(filepath) == 0:
    #     filepath = ["vesselAUG.mat"]
    ves=spio.loadmat(filepath[0])
    
    if ax == None:
        for i in range(ord('a'),ord('z')+1):
            plt.plot(ves["ves"][chr(i)][0][0][:,0],ves["ves"][chr(i)][0][0][:,1],
                     color,linewidth=linewidth)
        for i in range(ord('a'),ord('g')+1):
            plt.plot(ves["ves"][chr(i)+chr(i)][0][0][:,0],ves["ves"][chr(i)+chr(i)][0][0][:,1],
                     color,linewidth=linewidth)
    else:
        for i in range(ord('a'),ord('z')+1):
            ax.plot(ves["ves"][chr(i)][0][0][:,0],ves["ves"][chr(i)][0][0][:,1],
                    color,linewidth=linewidth)
        for i in range(ord('a'),ord('g')+1):
            ax.plot(ves["ves"][chr(i)+chr(i)][0][0][:,0],ves["ves"][chr(i)+chr(i)][0][0][:,1],
                    color,linewidth=linewidth)

def MASTu(ax=None,FILD_R=1.6,color='k',linewidth=0.7):
    
    for path in paths:
        if os.path.exists(path):
            print('Vessel read from: ', path)
            filepath = [path+'/vesselMASTu.mat']

    # filepath = [s+"/vesselMASTu.mat" for s in sys.path if os.path.isfile(s+"/vesselMASTu.mat")]
    # if len(filepath) == 0:
    #     filepath = ["vesselMASTu.mat"]
    
    ves=spio.loadmat(filepath[0],struct_as_record=False)#squeeze_me=True,
    # return ves
    for key in ves.keys():
        if key != "__header__" and key != "__globals__" and key != "__version__":
            plot_MASTu_square(ves[key][0][0], ax, color,linewidth=linewidth)

    plot_MASTu_FILD(FILD_R, ax, color)

def plot_MASTu_square(ves, ax=None, color='k',linewidth=0.7):
    
    for i in range(ves.centreR.size):
        if ves.shapeAngle1[0][i] == 0 and ves.shapeAngle2[0][i] == 0:
            x = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
                          ves.centreR[0][i]+ves.dR[0][i]/2,
                          ves.centreR[0][i]-ves.dR[0][i]/2,
                          ves.centreR[0][i]-ves.dR[0][i]/2,
                          ves.centreR[0][i]+ves.dR[0][i]/2])
            y = np.array([ves.centreZ[0][i]-ves.dZ[0][i]/2,
                          ves.centreZ[0][i]+ves.dZ[0][i]/2,
                          ves.centreZ[0][i]+ves.dZ[0][i]/2,
                          ves.centreZ[0][i]-ves.dZ[0][i]/2,
                          ves.centreZ[0][i]-ves.dZ[0][i]/2])
            # A = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2])
            # B = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2])
            # C = np.array([ves.centreR[0][i]-ves.dR[0][i]/2,
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2])
            # D = np.array([ves.centreR[0][i]-ves.dR[0][i]/2,
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2])

        elif ves.shapeAngle1[0][i] == 0:
            x = np.array([ves.centreR[0][i]+ves.dR[0][i]/2
                          -ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
                          ves.centreR[0][i]+ves.dR[0][i]/2
                          +ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
                          ves.centreR[0][i]-ves.dR[0][i]/2
                          +ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
                          ves.centreR[0][i]-ves.dR[0][i]/2
                          -ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
                          ves.centreR[0][i]+ves.dR[0][i]/2
                          -ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i]))])
            y = np.array([ves.centreZ[0][i]-ves.dZ[0][i]/2,
                          ves.centreZ[0][i]+ves.dZ[0][i]/2,
                          ves.centreZ[0][i]+ves.dZ[0][i]/2,
                          ves.centreZ[0][i]-ves.dZ[0][i]/2,
                          ves.centreZ[0][i]-ves.dZ[0][i]/2])
            # A = np.array([ves.centreR[0][i]+ves.dR[0][i]/2
            #               -ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2])
            # B = np.array([ves.centreR[0][i]+ves.dR[0][i]/2
            #               +ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2])
            # C = np.array([ves.centreR[0][i]-ves.dR[0][i]/2
            #               +ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2])
            # D = np.array([ves.centreR[0][i]-ves.dR[0][i]/2
            #               -ves.dZ[0][i]/2/np.tan(np.deg2rad(ves.shapeAngle2[0][i])),
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2])
        else:
            x = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
                          ves.centreR[0][i]+ves.dR[0][i]/2,
                          ves.centreR[0][i]-ves.dR[0][i]/2,
                          ves.centreR[0][i]-ves.dR[0][i]/2,
                          ves.centreR[0][i]+ves.dR[0][i]/2])
            y = np.array([ves.centreZ[0][i]-ves.dZ[0][i]/2
                          +ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i])),
                          ves.centreZ[0][i]+ves.dZ[0][i]/2
                          +ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i])),
                          ves.centreZ[0][i]+ves.dZ[0][i]/2
                          -ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i])),
                          ves.centreZ[0][i]-ves.dZ[0][i]/2
                          -ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i])),
                          ves.centreZ[0][i]-ves.dZ[0][i]/2
                          +ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i]))])
            # A = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2
            #               +ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i]))])
            # B = np.array([ves.centreR[0][i]+ves.dR[0][i]/2,
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2
            #               +ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i]))])
            # C = np.array([ves.centreR[0][i]-ves.dR[0][i]/2,
            #               ves.centreZ[0][i]+ves.dZ[0][i]/2
            #               -ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i]))])
            # D = np.array([ves.centreR[0][i]-ves.dR[0][i]/2,
            #               ves.centreZ[0][i]-ves.dZ[0][i]/2
            #               -ves.dR[0][i]/2*np.tan(np.deg2rad(ves.shapeAngle1[0][i]))])
        if ax == None:
            plt.plot(x,y,color,linewidth=linewidth)
        else:
            ax.plot(x,y,color,linewidth=linewidth)

def plot_MASTu_FILD(R,ax=None, color='k'):

    vstruct = spio.matlab.mio5_params.mat_struct()
 
    vstruct.dR = np.array([[1.8-R,0.2, 0.04, 0.04,0.015]])
    vstruct.dZ = np.array([[0.15,0.15,0.115,0.09,0.002]])
    vstruct.centreR = np.array([[0.1+(1.8-R)/2,0,-0.1+0.04/2,-0.1+0.04/2,-0.1+0.04-0.015/2]])
    vstruct.centreZ = np.array([[0,0,0.075-0.115/2,0.075-0.115+0.090/2,0.075-0.115+0.090+0.002/2]])
    
    vstruct.centreR = vstruct.centreR+R+0.1
    vstruct.centreZ = vstruct.centreZ+0.109
    
    vstruct.shapeAngle1 =np.array([[0,0,0,0,0]])
    vstruct.shapeAngle2 =np.array([[0,0,0,0,0]])

    plot_MASTu_square(vstruct, ax, color)
