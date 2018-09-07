#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 15:45:40 2018

@author: abauville
"""

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
from matplotlib.colors import Normalize
import matplotlib.colors
#from scipy import interpolate
#from scipy import ndimage

from Units import *
import time




## Create the folder tree
# ================================
Weak = 10
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/NoTopo/Beta0/Weak%i/" % Weak
#superRootFolder = "/work/G10501/abauville/Paper_Decollement/NoTopo/Beta0/Weak%i/" % Weak
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"
#superDirList = os.listdir(superRootFolder)
superDirList = ["Weak60/Lambda60"]
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

Compute = True
    
nSim = len(superDirList)
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
    
    Setup = Output.readInput(rootFolders[0] +  'Input/input.json')
    ix0_surf = 1
    ix1_surf = Setup.Grid.nxC+1
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
    
    slope = []
    
    
    for iSim in range(iSim0,nSim):
        nSteps = nStepsList[iSim]
        i0 = 0#0#jump-2
        jump = 1
        
        nProcessedSteps = len(np.arange(i0,nSteps,jump))
        xFront.append(np.zeros(nProcessedSteps))
        xBase.append(np.zeros(nProcessedSteps))
        slope.append(np.zeros(nProcessedSteps))
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
            rawData  = Output.getData(dataFolder + 'phase.bin',True)
            phase = rawData.data
            H = 1.0
            ny = rawData.ny
            ymin = rawData.ymin/H
            ymax = rawData.ymax/H
            nx = rawData.nx
            xmin = rawData.xmin/H
            xmax = rawData.xmax/H
            dx = (xmax-xmin)/(nx-1)
            ymin = rawData.ymin/H
            ymax = rawData.ymax/H
            dy = (ymax-ymin)/(ny-1)
            #mask = phase == 0
            
            #strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
            strain  = Output.getData(dataFolder + 'strain.bin',True).data
            
            

            A = np.where(strain[ix0_surf:ix1_surf,iy0_surf:iy1_surf]>0.25)
            B = np.where(strain[ix0_base:ix1_base,iy0_base:iy1_base]>0.25)
            
            
            
            strainFront[iSub][it] = np.max(strain[ix0_surf:ix1_surf,iy0_surf:iy1_surf],1)
            strainBase [iSub][it] = np.max(strain[ix0_base:ix1_base,iy0_base:iy1_base],1)
            
            iFront = 0
            try:
                iFront           = A[0][0]
                xFront[iSub][it] = A[0][0]
                xBase[iSub][it]  = B[0][0]
            except IndexError:
                daijoubu = 1
                    
            timeList[iSim][it] = timeSim
            
            
            Y = np.ones(phase.shape) * dy
            Y = np.cumsum(Y,1) + ymin
            Y[phase==0] = 0.0
            Topo = np.max(Y,1)
            iMaxTopo = np.argmax(Topo)
            if iMaxTopo == 0:
                iMaxTopo = 1
            if iMaxTopo<=iFront:
                iMaxTopo = nx
            TopoPrism = Topo[iFront:iMaxTopo]
            x = np.arange(iMaxTopo-iFront) * dx
            
            p = np.polyfit(x,TopoPrism,1)
            slope[iSub][it] = np.arctan(p[0])
            
            
#            phase  = Output.getData(dataFolder + 'phase.bin',True).data
            
            it+=1
        # end iStep
        iSub+=1
        
    #end iSim
#    for iSim in range(iSim0,nSim):
#        plt.plot(timeList[iSim],xFront[iSim],"-",color=Colors[iSim])
#        plt.plot(timeList[iSim],xBase[iSim],"--",color=Colors[iSim])
#    plt.legend(("Front %s" % superDirList[0], "Base %s" % superDirList[0],"Front %s" % superDirList[1], "Base %s" % superDirList[1]))
    
    
#    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/SurfaceAngle.npz",
    np.savez("./Data/SurfaceAngle_Weak%i_%s.npz" % (Weak,superDirList[0]),
             strainFront = strainFront,
             strainBase = strainBase,
             xFront = xFront,
             xBase = xBase, 
             timeList = timeList,
             slope = slope
             )
    
else: #if Compute   
    
#    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/SurfaceAngle.npz");
    loadedData = np.load("./Data/SurfaceAngle_Weak%i_%s.npz" % (Weak,superDirList[0]));
    strainFront = loadedData["strainFront"][()]
    strainBase = loadedData["strainBase"][()]
    xFront = loadedData["xFront"][()]
    xBase = loadedData["xBase"][()]
    timeList = loadedData["timeList"][()]
    slope = loadedData["slope"][()]
    
#end Compute
    
    
    
cdictRed = {'red':  ((0.0 , 1.0, 0.0),
                    (1.0 , 1.0, 0.0)),

          'green': ((0.0 , 1.0, 1.0),
                    (1.0 , 0.0, 1.0)),
 
          'blue':  ((0.0 , 1.0, 1.0),
                    (1.0 , 0.0, 0.0))
    }
cdictGray = {'red':  ((0.0 , 1.0, 0.0),
                     (1.0 , 0.0, 0.0)),

            'green': ((0.0 , 1.0, 1.0),
                      (1.0 , 0.0, 1.0)),

            'blue':  ((0.0 , 1.0, 1.0),
                      (1.0 , 0.0, 0.0))
    }


CMAP = np.zeros((256,4));
CMAP[:,0] = np.linspace(1.0,1.0,256);
CMAP[:,1] = np.linspace(1.0,0.0,256);
CMAP[:,2] = np.linspace(1.0,0.0,256);
CMAP[:,3] = np.linspace(0.25,1.0,256);
CMAP[0,3] = 0.0;


myCMAP =  matplotlib.colors.ListedColormap(CMAP,"RedTransparent",256);
#CMAP = LinearSegmentedColormap('RedTransparent', cdictRed)
plt.register_cmap(cmap=myCMAP)

CMAP = np.zeros((256,4));
CMAP[:,0] = np.linspace(1.0,0.0,256);
CMAP[:,1] = np.linspace(1.0,0.0,256);
CMAP[:,2] = np.linspace(1.0,0.0,256);
CMAP[:,3] = np.linspace(0.25,1.0,256);
CMAP[0,3] = 0.0;


myCMAP =  matplotlib.colors.ListedColormap(CMAP,"GrayTransparent",256);
#CMAP = LinearSegmentedColormap('RedTransparent', cdictRed)
plt.register_cmap(cmap=myCMAP)


plt.clf()
plt.set_cmap("GrayTransparent")

Hsed = 1.0e3
dx = (Setup.Grid.xmax-Setup.Grid.xmin)/Setup.Grid.nxC / Hsed

for iSim in range(iSim0,nSim):
    plt.plot(timeList[iSim],slope[iSim]*180.0/np.pi,'.')
    smoothSlope = np.zeros(slope[iSim].shape)
    smoothWindowSize = 3
    for i in range(smoothSlope.size):
        i0 = i-smoothWindowSize
        i1 = i+smoothWindowSize+1
        
        if i0<0:
            i0 = 0
        if i1>smoothSlope.size:
            i1 = smoothSlope.size
        
        smoothSlope[i] = np.mean(slope[iSim][i0:i1])
        
    plt.plot(timeList[iSim],slope[iSim]*180.0/np.pi,'.')
    plt.plot(timeList[iSim],smoothSlope*180.0/np.pi,'-')
    
    
    
    
    
    
    
#    short = pushVel*timeList[iSim]+1e-12
#    normDispBase = np.zeros(strainBase[iSim].shape)
#    normDispFront = np.zeros(strainBase[iSim].shape)
#    for i in range(len(short)):
#        normDispBase[i,:] = strainBase[iSim][i,:]*dx#/short[i]
#        normDispFront[i,:] = strainFront[iSim][i,:]*dx#/short[i]
##    plt.subplot(2,2,iSim*2+1)
#    plt.subplot(2,1,iSim+1)
##    plt.imshow(strainBase[iSim],vmin=0.5,vmax=20.0)
#    plt.imshow(normDispBase,vmin=0.01,vmax=0.5)
##    plt.pcolormesh(strainBase[iSim],vmin=0.5,vmax=30.0)
#    plt.set_cmap("GrayTransparent")
#    plt.colorbar()
#    
##  
##    plt.subplot(2,2,iSim*2+1+1)
##    plt.subplot(2,1,1)
##    plt.imshow(strainFront[iSim],vmin=0.5,vmax=10.0)
#    plt.imshow(normDispFront,vmin=0.025,vmax=0.25)
##    plt.pcolormesh(strainFront[iSim],vmin=0.5,vmax=10.0,edgecolors=(0.0,0.0,0.0,0.0))
#    plt.set_cmap("RedTransparent")
#    plt.colorbar()
#    plt.axis("auto")
##     
            
            
            
##    plt.set_cmap("GrayTransparent")
#    CMAP = plt.get_cmap()
#    CMAP._lut[:-1,3] = np.linspace(0.0,1.0,len(CMAP._lut[:-1,3]))
            
            
            
            
            
            
            
            
            
            
            
            
    