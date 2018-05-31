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


superRootFolder = "/Users/abauville/Output/Paper_Decollement/Test/Output/"

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

rootFolder = superRootFolder + superDirList[0] + "/"






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
iSim = 0
rootFolder = superRootFolder + superDirList[iSim] + "/"

DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")


dataSet     = Output.getData(rootFolder + DirList[0] + '/khi.bin')

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



nSteps = len(DirList)
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


i0 = 0
jump = 1


#
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    os.system("mkdir /home/abauvill/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
#elif sys.platform == "darwin":
#    # OS X
#    os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
    

iStep = 10
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
nSim = 1
iSub = 0

plt.clf()

#plt.colorbar()
it=-1
iSub = 0

Compute = True
if Compute:    
    for iStep in range(i0,nSteps,jump):
        for iSim in range(0,nSim):
            iSub+=1
            it+=1
#            rootFolder = superRootFolder + superDirList[iSim] + "/"
            rootFolder = superRootFolder + superDirList[0] + "/"
            if iSim == 1:
                rootFolder = rootFolder.replace("_phib30","_phib20")
            print(rootFolder)
            outFolder = "Out_%05d" % iStep #DirList[iStep]
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
            
            PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data
            PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data
            
            Hsed = 2.0*km
            PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data/Hsed
            PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data/Hsed
            
            PartPhase  = Output.getParticleData(dataFolder + 'particles_phase.bin',True).data
            PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data
            PartTimeLastPlastic  = Output.getParticleData(dataFolder + 'particles_timeLastPlastic.bin',True).data / yr


            PartStrainLogical = PartStrain>1.0

            plt.clf()
#            plt.pcolor(xv,yv,strainRate)
            
            plt.scatter(PartX,PartY,s=PartStrainLogical*100,c=PartTimeLastPlastic)
#            plt.scatter(PartX,PartY,c=PartTimeLastPlastic)
            plt.set_cmap("wyrk")
#            plt.set_cmap("Pastel1")
#            plt.scatter(PartX,PartY,s=PartPhase,c=PartYIni,vmin=0.0,vmax=1.0)
            plt.colorbar()
            plt.axis("equal")

            plt.pause(0.0001)
            
            
#plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/EndMembers",dpi=800)
    