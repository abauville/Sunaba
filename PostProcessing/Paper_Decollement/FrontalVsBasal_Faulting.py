#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 10:52:45 2018

@author: abauville
"""

# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan, arctan
from PaperDecollement_Utils import getColormap, get_XYandPattern
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
from matplotlib.colors import Normalize, to_rgba
#from scipy import interpolate
#from scipy import ndimage

from Units import *
import time




## Create the folder tree
# ================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/NoTopo/Beta0/Weak10/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
#if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]
#rootFolder += subFolder + "/"
DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")


rootFolders = [''] * len(superDirList) # initialize to the right size
nStepsList = np.zeros(len(superDirList),dtype='int16');
for i in range(len(superDirList)):
    rootFolders[i] = superRootFolder + superDirList[i] + "/Output/"
#    subFoldersList = os.listdir(rootFolders[i])
#    try:
#        subFoldersList.remove('.DS_Store')
#    except ValueError:
#        Dummy=0
#    rootFolders[i] += subFoldersList[0] + "/"
#    
    stepFolderList = os.listdir(rootFolders[i])
    try:
        stepFolderList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
#Setup = Output.readInput(rootFolder +  'Input/input.json')
#s = Setup.Description
#Char = Setup.Char

pushVel = 10.0*cm/yr

Compute = False
    
nSim = 2#len(superDirList)
#    nSim = 11
iSim0 = 0#nSim-1
if Compute:    
    it=-1
#    goldenRatio = (1.0+np.sqrt(5)/2.0)
#    renderMatplotlib = 0
#    renderMayavi = 1
#    renderer = 0 # 0 Matplotlib; 1 Mayavi
#    if renderer == renderMatplotlib:
#        # set figure        
#        cm2inch = 0.393701
#        figW = 2.0*18.0 * cm2inch
#        figH = figW/goldenRatio
#        plt.close()
#        plt.figure(1,figsize=(figW,figH))
#        mngr = plt.get_current_fig_manager()
#    #     to put it into the upper left corner for example:
#        FigCoord = mngr.window.geometry().getRect()
#        mngr.window.setGeometry(FigCoord[0]-2500,FigCoord[1]-200,FigCoord[2],FigCoord[3])
#        Ax1 = plt.axes([0,0,1,1])
#
#                
#    # define the limits of the axis
#    padAx = 0.1
#    ypMin = -padAx
#    ypMax = 6.0
#    xpMin = -(ypMax-ypMin)*goldenRatio+padAx
#    xpMax = padAx            
    
    
    # set a font dict
    font = {'family': 'Montserrat',
    'weight': 'bold',
    'size': 16
    }
    
    

    Hsed = 2.0 * km
    
    
#    iStep = i0
    
        
    Hgrid = 64 # Thickness H in terms of grid points
    
    
    ix0_surf = 0
    ix1_surf = Setup.Grid.nxC
    iy0_surf = Hgrid-5
    iy1_surf = Hgrid+5
    
    
    ix0_base = ix0_surf
    ix1_base = ix1_surf
    iy0_base = 0
    iy1_base = 5


    iSub = 0
    Colors = ([1.0,0.2,0.4],[0.2,0.4,1.0])
    plt.clf()
    xFront = []
    xBase = []
    timeList = []
    
    strainFront = []
    strainBase = []
    
    
    
    
    for iSim in range(iSim0,nSim):
        nSteps = nStepsList[iSim]
        i0 = 0#jump-2
        jump = 5
        
        nProcessedSteps = len(np.arange(i0,nSteps,jump))
        xFront.append(np.zeros(nProcessedSteps))
        xBase.append(np.zeros(nProcessedSteps))
        timeList.append(np.zeros(nProcessedSteps))
        
        strainFront.append(np.zeros((nProcessedSteps,ix1_surf-ix0_surf)))
        strainBase.append(np.zeros((nProcessedSteps,ix1_surf-ix0_surf)))
        
        frame = int(i0/jump)
    
        it = 0
        for iStep in range(i0,nSteps,jump):
        
            tic = time.time()
            
            
            
            ## Set Folder
            # ================================
            
#            outFolder = os.listdir(rootFolder)[-1]
            outFolder = 'Out_%05d' % (iStep)
            dataFolder = rootFolders[iSim] + outFolder + "/"
            #print(rootFolders[iSim])
            if (it%10==0):
                print("iStep = %i/%i" % (iStep, nSteps-1))
               
            ## Get Data and Pattern
            # ================================
            Setup = Output.readInput(rootFolders[iSim] +  'Input/input.json')
            Char = Setup.Char
            
            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            #phase = Output.getData(dataFolder + 'phase.bin').data
            #mask = phase == 0
            
            #strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
            strain  = Output.getData(dataFolder + 'strain.bin',True).data
            
            
            
            A = np.where(strain[ix0_surf:ix1_surf,iy0_surf:iy1_surf]>0.25)
            B = np.where(strain[ix0_base:ix1_base,iy0_base:iy1_base]>0.25)
            
            strainFront[iSub][it] = strain[ix0_surf:ix1_surf,iy0_surf+1]
            strainBase [iSub][it] = strain[ix0_base:ix1_base,iy0_base+1]
            
            try:
                xFront[iSim][it] = ix1_surf-A[0][0]
                xBase[iSim][it]  = ix1_surf-B[0][0]
            except IndexError:
                daijoubu = 1
                    
            timeList[iSim][it] = timeSim/(365.25*24*3600)
            it+=1
                    
        # end iStep
        iSub+=1
        
    #end iSim
#    for iSim in range(iSim0,nSim):
#        plt.plot(timeList[iSim],xFront[iSim],"-",color=Colors[iSim])
#        plt.plot(timeList[iSim],xBase[iSim],"--",color=Colors[iSim])
#    plt.legend(("Front %s" % superDirList[0], "Base %s" % superDirList[0],"Front %s" % superDirList[1], "Base %s" % superDirList[1]))
    
    
    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/FrontalVsBasal.npz",
             strainFront = strainFront,
             strainBase = strainBase,
             xFront = xFront,
             xBase = xBase
             )
    
else: #if Compute   
    
    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/FrontalVsBasal.npz");
    strainFront = loadedData["strainFront"][()]
    strainBase = loadedData["strainBase"][()]
    xFront = loadedData["xFront"][()]
    xBase = loadedData["xBase"][()]
    
    
#end Compute
    
    
    
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
               (0.5 , 0.0, 0.0),
               (0.75, 1.0, 1.0),
               (1.0 , 1.0, 1.0)),

     'green': ((0.0 , 1.0, 1.0),
               (0.25, 0.0, 0.0),
               (0.75, 0.0, 0.0),
               (1.0 , 1.0, 1.0)),

     'blue':  ((0.0 , 1.0, 1.0),
               (0.25 , 1.0, 1.0),
               (0.5 , 0.0, 0.0),
               (1.0 , 0.0, 0.0))
    }
CMAP = LinearSegmentedColormap('Test', cdict1)
plt.register_cmap(cmap=CMAP)
plt.set_cmap("Test")
plt.clf()
for iSim in range(iSim0,nSim):

    plt.subplot(2,2,iSim*2+1)
    plt.pcolor(strainFront[iSim])
    plt.colorbar()
    
    plt.subplot(2,2,iSim*2+1+1)
    plt.pcolor(strainBase[iSim])
    plt.colorbar()
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
    