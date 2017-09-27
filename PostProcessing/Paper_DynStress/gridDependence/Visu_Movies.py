# Output reading test

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
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

#rootFolder = "/Users/abauville/StokesFD_Output/ViscoElasticBuildUp/"
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
#superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/gridDependence/Test/"
superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_Stronger_Seed_Save/"
superDirList = os.listdir(superRootFolder)
superDirList.remove('.DS_Store')
nSim = len(superDirList)
#superDirList.remove('Input')

ResFacList = np.zeros(nSim)
for i in range(nSim):
    ResFacList[i] = float(superDirList[i][-7:])

I = ResFacList.argsort()
I.astype(np.int)
I = I[::-1]
ResFacList = ResFacList[I]
superDirList = [ superDirList[i] for i in I]

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
iSim = nSim-1
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
dx      = Wbox/nx
dy      = Hbox/ny



nSteps = len(DirList)
jump = 1
nSteps = int(nSteps/jump)
time_t = np.zeros(nSteps)
TauII_t = np.zeros(nSteps)
P_t = np.zeros(nSteps)
EII_t = np.zeros(nSteps)


#nyTotal = 100 # to exclude the sticky air
#jumpy = 51
#ny = int(np.ceil(nyTotal/jumpy))
#Ix_y = np.zeros(ny);

#TauII_t_y = np.zeros([nSteps,ny])
#P_t_y = np.zeros([nSteps,ny])
#EII_t_y = np.zeros([nSteps,ny])


#
#P_dict      = dict()
#EII_dict    = dict()
#TauII_dict  = dict()
#time_dict   = dict()

# Extract Data from file
# =======================
# Get Ix for each Sim

# Test
#nSim = 3
#nSim = 9


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
stressUnit = Setup.Physics.Pback
strainRateUnit = np.abs(Setup.BC.Stokes.backStrainRate)


plt.figure(6)
# plot
ResFac = ResFacList[iSim]
iyCell = int(100/2*ResFac)+1
plt.clf()
plt.ion()
i0 = 90000
jump = 10
time_t = np.zeros(len(range(i0,nSteps,jump)))
Pfault_t    = np.zeros(len(range(i0,nSteps,jump)))
Pfar_t      = np.zeros(len(range(i0,nSteps,jump)))
it=-1

os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
for iStep in range(i0,nSteps,jump):
    it+=1
    outFolder = "Out_%05d" % iStep #DirList[iStep]
    print("iStep = %i/%i" % (iStep, nSteps-1))
    # index of the first node that register the minimum of khi (where khi is active)
    # Set file
    # =====================
    Setup = Output.readInput(rootFolder +  'Input/input.json')
    Char = Setup.Char
    CharExtra = Input.CharExtra(Char)

    dataFolder = rootFolder + outFolder + "/"
    thisData = Output.getData(dataFolder + 'P.bin')
    P = thisData.data * CharExtra.stress
    ixCell = np.argmin(thisData.data[:,iyCell]) # find the minimum khi
    
    
    thisData = Output.getData(dataFolder + 'strainRate.bin')
    strainRate = thisData.data * 1.0/Char.time
    
    State       = Output.readState(dataFolder + "modelState.json")
    timeSim   = (State.time+ State.dt) * Setup.Char.time
    
#    if (np.mod(iStep,100)==0):
#        print("iStep = %i/%i" % (iStep, nSteps))
#        print("ixCell = %i" %ixCell)
    
    plt.clf()
   
    
    cAx_PMin = 1.0
    cAx_PMax = 3.0
    
    cAx_srMin = -2.0
    cAx_srMax = +2.0

    plt.subplot(2,2,1)
    plt.pcolor(xv - dx/2.0,yv - dy/2.0,P/stressUnit,vmin=cAx_PMin,vmax=cAx_PMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.colorbar()
    plt.axis('equal')
    plt.set_cmap("Pressure")
    
    plt.subplot(2,2,2)
    plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.colorbar()
    plt.axis('equal')
    plt.set_cmap("Pressure")
    
    time_t[it] = timeSim
    Pfault_t [it]   = P[ixCell, iyCell]
    Pfar_t [it]     = P[-2, iyCell]
    
    plt.subplot(2,1,2)
    
    
    plt.plot(time_t[:it+1]/timeUnit, Pfault_t[:it+1]/stressUnit,'.k',markerSize=1)
    plt.plot(time_t[:it+1]/timeUnit, Pfar_t[:it+1]/stressUnit,'.b',markerSize=1)
    plt.plot([0,1.5],np.array([P_lim,P_lim])/stressUnit,'--r')
    plt.axis([0,1.5,0,3])
#    plt.show()
#    plt.pause(0.0001)
    plt.savefig("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim] + "/Frame_%06d.png" % it, dpi=300)
    
    
#    
    
##end iStep
        
        
        
        
        
        
        
        
        
        

