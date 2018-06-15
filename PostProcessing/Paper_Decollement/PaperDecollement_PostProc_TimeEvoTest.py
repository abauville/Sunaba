or#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 16:05:51 2018

@author: abauville
"""

# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan, arctan
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

#import MaterialsDef as Material


##             Units
## =====================================
m       = 1.0
s       = 1.0
K       = 1.0
kg      = 1.0

cm      = 0.01      * m
km      = 1000.0    * m

mn      = 60        * s
hour    = 60        * mn
day     = 24        * hour
yr      = 365       * day
Kyr     = 1e3       * yr
Myr     = 1e6       * yr

Pa = 1.0
MPa = 1e6*Pa

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
CMAP = LinearSegmentedColormap('Pressure', cdict1)
plt.register_cmap(cmap=CMAP)

# k-r-y
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5 , 1.0, 1.0),
                   (1.0 , 1.0, 1.0)),

         'green': ((0.0 , 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (1.0 , 1.0, 1.0)),

         'blue':  ((0.0 , 0.0, 0.0),
                   (1.0 , 0.0, 0.0))
        }
CMAP = LinearSegmentedColormap('kry', cdict1)
plt.register_cmap(cmap=CMAP)

# k-r-y
cdict1 = {'red':  ((0.0 , 1.0, 1.0),
                   (0.33, 1.0, 1.0),
                   (0.66, 1.0, 1.0),
                   (1.0 , 0.0, 0.0)),

         'green': ((0.0 , 1.0, 1.0),
                   (0.33, 1.0, 1.0),
                   (0.66, 0.0, 0.0),
                   (1.0 , 0.0, 0.0)),

         'blue':  ((0.0 , 1.0, 1.0),
                   (0.33, 0.0, 0.0),
                   (0.66, 0.0, 0.0),
                   (1.0 , 0.0, 0.0))
        }
CMAP = LinearSegmentedColormap('wyrk', cdict1)
plt.register_cmap(cmap=CMAP)


cdict1 = {'red':  ((0.0 , 1.0, 1.0),
                   (0.25, 0.25, 0.25),
                   (0.5 , 1.0, 1.0),
                   (0.75, 1.0, 1.0),
                   (1.0 , 0.25, 0.0)),

         'green': ((0.0 , 1.0, 1.0),
                   (0.25, 1.0, 1.0),
                   (0.5 , 1.0, 1.0),
                   (0.75, 0.0, 0.0),
                   (1.0 , 0.0, 0.0)),

         'blue':  ((0.0 , 1.0, 1.0),
                   (0.25, 1.0, 1.0),
                   (0.5 , 0.25,0.25),
                   (0.75, 0.0, 0.0),
                   (1.0 , 0.0, 0.0))
        }
CMAP = LinearSegmentedColormap('wcyrk', cdict1)
plt.register_cmap(cmap=CMAP)






# Create the colormap many random colors
nRandColors  = 12
condensedCMAPtemp = np.random.rand(nRandColors,4) * 0.5+.1
condensedCMAPtemp[:,0] += .4
condensedCMAPtemp[:,1] += .25
#condensedCMAPtemp[:,3] += .1
condensedCMAP = np.zeros((2*nRandColors,4))
condensedCMAP[0::2,:] = condensedCMAPtemp
condensedCMAP[1::2,:] = condensedCMAPtemp*0.25
#condensedCMAP = np.concatenate(condensedCMAP,0.5*condensedCMAP)
condensedCMAP[:,-1] = 1.0

nColors = condensedCMAP.shape[0]



# Create the equivalent matplotlib colormap

#cdict1 = {'red':  np.zeros((nColors,3)),
#              'green': np.zeros((nColors,3)),
#              'blue': np.zeros((nColors,3))
#              }
#for i in range(nColors):
#    cdict1["red"]    [i,:] = np.array([i/(nColors-1),    condensedCMAP[i,0],   condensedCMAP[i,0]  ])
#    cdict1["green"][i,:] = np.array([i/(nColors-1),    condensedCMAP[i,1],   condensedCMAP[i,1]  ])
#    cdict1["blue"]  [i,:] = np.array([i/(nColors-1),    condensedCMAP[i,2],   condensedCMAP[i,2]  ])
#
#CMAP = LinearSegmentedColormap('xrandCol_DarkBands"', cdict1)
#plt.register_cmap(cmap=CMAP)

#condensedCMAP = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

CMAP = LinearSegmentedColormap.from_list('xrandCol_DarkBands',condensedCMAP)
plt.register_cmap(cmap=CMAP)







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
CharExtra = Input.CharExtra(Char)


# Read strain rate data
# =====================
#dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')


#dataType = 'P'
#iSim = 0
#rootFolder = superRootFolder + superDirList[iSim] + "/"

DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")


dataSet     = Output.getData(rootFolder + DirList[0] + '/P.bin')

C = Setup.MatProps['1'].cohesion
phi = Setup.MatProps['1'].frictionAngle

xmin    = dataSet.xmin * Setup.Char.length
xmax    = dataSet.xmax * Setup.Char.length
ymin    = dataSet.ymin * Setup.Char.length
ymax    = dataSet.ymax * Setup.Char.length
nx      = dataSet.nx
ny      = dataSet.ny
Wbox       = xmax-xmin 
Hbox       = ymax-ymin
dx      = Wbox/(nx-1)
dy      = Hbox/(ny-1)

xmin = xmin-dx/2.0
xmax = xmax-dx/2.0
ymin = ymin-dy/2.0
ymax = ymax-dy/2.0



nSteps = 747# len(DirList)
jump = 1
nSteps = int(nSteps/jump)


# Define grid
# =====================
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)

Hmatrix = 9000.0;


# Analytical yield stress
# =====================  
S3 = Setup.Physics.Pback
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
Sy_back = (S1-S3)/2.0
P_lim = (S1+S3)/2.0

S1 = 1.0/(1.0-sin(phi)) * (  2*0*cos(phi) + S3*(1+sin(phi))  )
Sy_back_C0 = (S1-S3)/2.0

# Units and characteristic values
# =====================
EII = np.abs(Setup.BC.Stokes.backStrainRate)
eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
G = Setup.MatProps['1'].G

t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
charTime = t

timeUnit = charTime
stressUnit = 1.0*MPa#Setup.Physics.Pback
strainRateUnit = 1.0# np.abs(Setup.BC.Stokes.backStrainRate)
strainUnit = 1.0


#i0 = 746
jump = 1


#
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    os.system("mkdir /home/abauvill/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
#elif sys.platform == "darwin":
#    # OS X
#    os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
    

iStep = 1
nSim = len(superDirList)



#nSteps = 5000
nSubPlots = 6
jump= 1#int(np.floor((nSteps-1)/(nSubPlots)))
i0 = nSteps-1#jump-2
nt = len(range(i0,nSteps,jump))

xminP = xmin*1.0
xmaxP = xmax
yminP = ymin
ymaxP = ymax*0.6

# ColorScales
# Strain Rate    
cAx_srMin = -13
cAx_srMax = -11
# Strain
cAx_sMin = 0.0
cAx_sMax = 10.0

axDict = dict()
nSim = len(superDirList)
iSub = 0
Hsed = 2e3


#plt.colorbar()
it=-1
iSub = 0

Compute = True
if Compute:    
    i0 = 100
    jump = 100
    nSteps = 701
    plt.figure(2)
    plt.clf()
    nSim = 2
    
for iStep in range(i0,nSteps,jump):
    for iSim in range(0,nSim):
        iSub+=1
        it+=1
        if iSim == 0:    
            ISim = 1
        else:
            ISim = 3
        rootFolder = superRootFolder + superDirList[ISim] + "/Output/"
        subFolder = os.listdir(rootFolder)
        try:
            subFolder.remove('.DS_Store')
        except ValueError:
            Dummy=0

        rootFolder += subFolder[0] + "/"
#            subFolder = os.listdir(rootFolder)[0]
        print(rootFolder)
        outFolder = "Out_%05d" % iStep #DirList[iStep]
#        outFolder = os.listdir(rootFolder)[-1]
        print("iStep = %i/%i" % (iStep, nSteps-1))
        # index of the first node that register the minimum of khi (where khi is active)
        # Set file
        # =====================
        
        
        Setup = Output.readInput(rootFolder +  'Input/input.json')
        #Char = Setup.Char
        #CharExtra = Input.CharExtra(Char)
    
        dataFolder = rootFolder + outFolder + "/"
        State       = Output.readState(dataFolder + "modelState.json")
        timeSim   = (State.time+ State.dt) * Setup.Char.time
        
#            phase = Output.getData(dataFolder + 'phase.bin').data
#            mask = phase == 0
        
#            strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
        
        PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data/Hsed
        PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data/Hsed
        
        PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data/Hsed
        PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data/Hsed
        
#            PartPhase  = Output.getParticleData(dataFolder + 'particles_phase.bin',True).data
        PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data
#            PartTimeLastPlastic  = Output.getParticleData(dataFolder + 'particles_timeLastPlastic.bin',True).data / yr


        # subsampling
        sr = 100 # sample rate
        PartX    = PartX[0::sr]
        PartY    = PartY[0::sr]
        PartXIni = PartXIni[0::sr]
        PartYIni = PartYIni[0::sr]
        PartStrain    = PartStrain[0::sr]


        # Geometrical pattern
                
        xmax = 0.0
        xmin = -20#np.min(PartXIni)
        PartXPat = PartXIni
        PartXPat -= xmin
        PartXPat /= xmax-xmin
        
#        nLayers = 8
#        YPattern = np.cos(PartYIni*1.0*np.pi*nLayers)
        nLayers = 1
        #XPattern = np.sin(PartXIni*1.0*np.pi*nLayers)
#        XPattern = np.sin(PartXPat*0.5*np.pi*nLayers)
        XPattern = PartXPat
        XPattern[XPattern<0.0] = 0.0
        PartPattern= 2.0*np.round(XPattern * (nRandColors-1))
        #PartPattern+= YPattern
        
        #plt.scatter(PartX,PartY,c=PartPattern)


        # Map strain of the colormap
        maxStrain = 5.0
        minStrain = 1.0
        PartStrain -= minStrain
        PartStrain[PartStrain<0.0] = 0.0
        PartStrain /= maxStrain-minStrain
        
        PartStrainLogical = PartStrain>1.0
        PartStrain[PartStrainLogical] = 1.0
        #PartPattern += 1.0*PartStrainLogical
        
        PartPattern += 1.0*PartStrain

        
        plt.subplot(np.ceil((nSteps-i0)/jump),2,iSub)
#            plt.pcolor(xv,yv,strainRate)
        
#            plt.scatter(PartX,PartY,s=PartStrainLogical*100,c=PartTimeLastPlastic)
        plt.scatter(PartX,PartY,c=PartPattern,s=15.0,vmin=0.0)
#        plt.set_cmap("wyrk")
        plt.set_cmap("xrandCol_DarkBands")
#        plt.set_cmap("Pastel1")
#            plt.scatter(PartX,PartY,s=PartPhase,c=PartYIni,vmin=0.0,vmax=1.0)
#        plt.colorbar()
             
        plt.axis("equal")
        plt.axis([-10,0,0,2.5])   
        
        if (iStep==i0):
            plt.title(superDirList[ISim])
        plt.pause(0.0001)
        
        
#plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/EndMembers",dpi=800)
    