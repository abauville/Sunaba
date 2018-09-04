
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
#from scipy import interpolate
#from scipy import ndimage

from Units import *



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
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Static2/Beta0/"
Weak = 50
wWater = ""#"_w_Water"
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/NoTopo/Beta0/Weak%i%s/" % (Weak, wWater)
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    
ProductionMode = False
#weakList = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
#weakList = [10, 20, 30, 40, 50, 60, 70, 80, 90]
weakList = [1, 5, 10, 20, 40, 60, 80]
#LambdaList = [40,50,60,80,90]
LambdaList = [00,40,60,80]


if ProductionMode:
    sampleRate = 1
    pointSize = 0.0002
else:
    sampleRate = 100
    pointSize = 0.4

superDirList = []
for weak in weakList:
    for Lambda in LambdaList:
        superDirList.append("Weak%02d/Lambda%02d" % (weak, Lambda))




rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]


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
    
    stepFolderList = os.listdir(rootFolders[i])
    try:
        stepFolderList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char

pushVel = 10.0*cm/yr


iSim = 8

outFolder = os.listdir(rootFolders[iSim])[-1]

dataFolder = rootFolders[iSim] + outFolder + "/"
print(rootFolders[iSim])
#print("iStep = %i/%i" % (iStep, nSteps-1))
print(outFolder)  
## Get Data and Pattern
# ================================
Char = Output.readInput(rootFolders[iSim] +  'Input/input.json').Char
timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
plt.figure(2)
PartX = []
PartY = []
PartPattern = []
plt.clf()
#strain  = Output.getData(dataFolder + 'strain.bin',True).data
phase = Output.getData(dataFolder + 'phase.bin',True).data
mask = phase==0
strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask=mask).data
#sII  = Output.getData(dataFolder + 'sigma_II.bin',True,mask=mask).data
strain  = Output.getData(dataFolder + 'strain.bin',True,mask=mask).data

#plt.pcolormesh(np.log10(strainRate.T),cmap="Perso",shading='Gouraud')
plt.contourf(np.log10(strainRate.T),1000,cmap="Perso")
#plt.contourf(np.log10(strainRate.T),[-14,-13,-12,-11,-10],cmap="Perso")
#plt.contourf(sII.T,1000,cmap="Perso")
#plt.contourf(strain.T,1000,cmap="Perso")
ax = plt.gca()
#ax.invert_yaxis()
plt.colorbar()

plt.axis([600,1200,0,300])
plt.axis("equal")


            
         