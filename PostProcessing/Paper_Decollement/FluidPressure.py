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


# Colormap Strain rate
# c-b-k-r-y
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
CMAP = LinearSegmentedColormap('Perso', cdict1)
plt.register_cmap(cmap=CMAP)

## Create the folder tree
# ================================
Weak = 20.0
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/NoTopo/Beta0/Weak%i/" % (Weak)
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
superDirList = ['Hc1.000_Lambda90']
    
    
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




nSim = 1#len(superDirList)
#    nSim = 11
iSim0 = nSim-1
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
        i0 = nSteps-1#jump-2
        jump = 1
        
        
        frame = int(i0/jump)
    
        it = 0
        for iStep in range(i0,nSteps,jump):
        
            tic = time.time()
            
            
            
            ## Set Folder
            # ================================
            
#            outFolder = os.listdir(rootFolder)[-1]
#            outFolder = 'Out_%05d' % (iStep)
            outFolder = os.listdir(rootFolders[iSim])[-1]
            dataFolder = rootFolders[iSim] + outFolder + "/"
            #print(rootFolders[iSim])
            if (it%10==0):
                print("iStep = %i/%i" % (iStep, nSteps-1))
               
            ## Get Data and Pattern
            # ================================
            Setup = Output.readInput(rootFolders[iSim] +  'Input/input.json')
            Char = Setup.Char
            CharExtra = Input.CharExtra(Char)
            
            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            rawData  = Output.getData(dataFolder + 'phase.bin',False)
            phase = rawData.data
            H = 1.0
            ny = rawData.ny
            ymin = rawData.ymin/H * Char.length
            ymax = rawData.ymax/H * Char.length
            nx = rawData.nx
            xmin = rawData.xmin/H * Char.length
            xmax = rawData.xmax/H * Char.length
            dx = (xmax-xmin)/(nx-1)
            dy = (ymax-ymin)/(ny-1)
            
            strain  = Output.getData(dataFolder + 'strain.bin').data
            Pressure  = Output.getData(dataFolder + 'P.bin',True).data
#            plt.pcolor(Pressure.T/1e6,vmin=0.0)
#            plt.axis("equal")
#            plt.colorbar()
#            plt.set_cmap("Perso")
            
            
            # Extract a profile
            ix = 700
            P_prof = Pressure[ix,:]
            phase_prof = phase[ix,:]
            strain_prof = strain[ix,:]
            P_prof[P_prof<0.0] = 0.0
#            P_prof = np.flipud(P_prof)
#            phase_prof = np.flipud(phase_prof)
#            strain_prof = np.flipud(strain_prof)
            phase_prof = np.flipud(phase_prof)
            H_prof = np.cumsum(phase_prof)*dy
            phase_prof = np.flipud(phase_prof)
            H_prof = np.flipud(H_prof)
            
            rho     = 2500.0
            rho_w   = 1000.0
            g = 9.81
            if superDirList[iSim][-2] == "a":
                Lambda = float(superDirList[iSim][-1:]) / 100.0
            else:
                Lambda = float(superDirList[iSim][-2:])  / 100.0
            
            
            WeakFac = Weak/100.0
            Sediment = Setup.MatProps["1"]


            Plitho = rho*g*H_prof
            Plitho_C = rho*g*H_prof + 0.5 * Sediment.cohesion * np.cos(Sediment.frictionAngle)
            Phydro = rho_w*g*H_prof
           
            
            Lambda_Weak = (1.0-WeakFac)*(  Lambda ) + WeakFac
            preFac = 0.05
            C1 = Sediment.strainWeakEnd
            C0 = Sediment.strainWeakStart
            Fac = 1.0-preFac-(strain_prof-C0)/(C1-C0)
            Fac[Fac<0] = 0.0
            Fac[Fac>1.0] = 1.0
            Lambda_prof = Lambda*Fac + (1.0-Fac)*Lambda_Weak
            
            Pf = P_prof*Lambda_prof
            Pf2 = Plitho*Lambda
            
            
            plt.clf()
            plt.plot(Phydro/1e6,H_prof,'-b')
            plt.plot(Plitho/1e6,H_prof,'-m')
            plt.plot(Plitho_C/1e6,H_prof,'--m')
            plt.plot(P_prof/1e6,H_prof,'-r')
            plt.plot(Pf/1e6,H_prof,'-k')
            plt.plot(Pf2/1e6,H_prof,'--k')
            
            plt.axis([0.0,np.max(Plitho_C/1e6),np.max(H_prof),0.0])
            
#            plt.clf()
#            plt.plot(strain_prof,H_prof)
#            plt.plot(Lambda_prof,H_prof)
#            plt.axis([0.0,1.0,0.0,5000])
#            plt.plot(strain_prof,H_prof)
#            plt.plot(phase_prof,H_prof)
            
            
            
        # end iStep
        iSub+=1
        
    # end iSim
# end if Compute            
            

            
            
            
            
            
            
            
    