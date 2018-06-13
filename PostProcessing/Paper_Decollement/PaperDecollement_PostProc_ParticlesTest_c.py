# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan, arctan
from PaperDecollement_SegmentColorMap import getColormap
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
from matplotlib.colors import Normalize
#from scipy import interpolate
from scipy import ndimage
import time

from Units import *





#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Test/Output/"
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Beta0/"

superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    

nSim = len(superDirList)
#superDirList.remove('Input')

#ResFacList = np.zeros(nSim)
#for i in range(nSim):
#    ResFacList[i] = float(superDirList[i][-7:])
#
#I = ResFacList.argsort()
#I.astype(np.int)
#I = I[::-1]
#ResFacList = ResFacList[I]
#superDirList = [ superDirList[i] for i in I]

rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[1]
rootFolder += subFolder + "/"








# Read parameters of this simulation
# =====================
#inputFile = Output.readJson(rootFolder + simFolder + inFolder + '/input.json')
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char


DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")



it=-1
iSub = 0

Compute = True
if Compute:    
#    plt.close()
    plt.figure(1)
    plt.clf()
    nSim = 6
    iSim0 = nSim-1
    Hsed = 2.0 * km
    nSteps = 746
    i0 = nSteps-1#jump-2
    iStep = i0
    jump = 1
    #for iStep in range(i0,nSteps,jump):
    
    for iSim in range(iSim0,nSim):
        iSub+=1
        it+=1
        rootFolder = superRootFolder + superDirList[iSim] + "/Output/"
        
        subFolder = os.listdir(rootFolder)
        try:
            subFolder.remove('.DS_Store')
        except ValueError:
            Dummy=0

        rootFolder += subFolder[0] + "/"
        print(rootFolder)
        outFolder = os.listdir(rootFolder)[-1]
        print("iStep = %i/%i" % (iStep, nSteps-1))
           
        dataFolder = rootFolder + outFolder + "/"
#        State       = Output.readState(dataFolder + "modelState.json")
#        timeSim   = (State.time+ State.dt) * Setup.Char.time
        
        
        # Load and subsample
        sr = 50 # sample rate        
        PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sr]/Hsed
        PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sr]/Hsed
        PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sr]/Hsed
        PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sr]/Hsed

        PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sr]

        
        # Geometrical pattern
        xmax = 0.0
        xmin = -22#np.min(PartXIni)

        # Create horizontal layer pattern with 2 values 0 or 1        
        nLayersY = 7
        YPattern = np.cos(PartYIni*1.0*np.pi*nLayersY+np.pi/2.0)
        maxVal= np.max(YPattern)
        YPattern = np.round((YPattern+maxVal)/(2.0*maxVal))
        
        # Create vertical layer pattern with 2 values 0 or 1
        nLayersX = 11
        nColors = nLayersX
        XPattern = PartXIni
#        XPattern -= xmin
        XPattern /= xmax-xmin
        XPattern = -XPattern
#        XPattern[XPattern<0.0] = 0.0
        PartPattern= 4.0*np.round(XPattern * (nLayersX) )
        PartPattern+= 2.0*YPattern
        
        
        # Create the colormap many random colors
        if iSub==1:
            CMAP = getColormap(nColors,"myColorMap")
            plt.register_cmap(cmap=CMAP)
            plt.set_cmap("myColorMap")
        
            
        
        
        # Map strain of the colormap
        maxStrain = 5.0
        minStrain = 1.0
        PartStrain -= minStrain
        PartStrain[PartStrain<0.0] = 0.0
        PartStrain /= maxStrain-minStrain
        
        PartStrainLogical = PartStrain>1.0
        PartStrain[PartStrainLogical] = 1.0        
        PartPattern += 1.0*PartStrain

        
#        plt.subplot(np.ceil(nSim/2.0),2,iSim+1)
        plt.scatter(PartX,PartY,c=PartPattern,s=15.0,vmin=0.0,vmax=4*nColors-1)

#        plt.set_cmap("xrandCol_DarkBands")
        
        plt.colorbar()
             
        plt.axis("equal")
        plt.axis([-10,0,0,2.5])   
       
        plt.title(superDirList[iSim])
        plt.pause(0.0001)
        
        
#plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/EndMembers",dpi=800)
    