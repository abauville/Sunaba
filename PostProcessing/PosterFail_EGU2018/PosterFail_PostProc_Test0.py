# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
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

#rootFolder = "/Users/abauville/StokesFD_Output/ViscoElasticBuildUp/"
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
#superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/gridDependence/Test/"
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    superRootFolder = "/home/abauvill/StokesFD_Output/Paper_DynStress/Output/dtDependence/Test_Stronger_Seed_NoSubGridDiff/"
#elif sys.platform == "darwin":
#    # OS X
##    superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/gridDependence/Test/"
#    superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_NoAdv_NoInterp_adaptative/"


superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity/Output/"
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
dx      = Wbox/nx
dy      = Hbox/ny



nSteps = len(DirList)
jump = 1
nSteps = int(nSteps/jump)
#nSteps = 1
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
stressUnit = 1.0*MPa#Setup.Physics.Pback
strainRateUnit = 1.0# np.abs(Setup.BC.Stokes.backStrainRate)
strainUnit = 1.0

plt.figure(1)
# plot
#ResFac = ResFacList[iSim]
#iyCell = int(100/2*ResFac)+1
plt.clf()
plt.ion()
i0 = 0
jump = 1
time_t = np.zeros(len(range(i0,nSteps,jump)))
Pfault_t    = np.zeros(len(range(i0,nSteps,jump)))
Pfar_t      = np.zeros(len(range(i0,nSteps,jump)))
maxStrain_Sim = np.zeros(nSim)
it=-1

#
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    os.system("mkdir /home/abauvill/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
#elif sys.platform == "darwin":
#    # OS X
#    os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
    

iStep = 100
nSim = len(superDirList)
plt.figure(1)
plt.clf()
plt.figure(2)
plt.clf()
#for iStep in range(i0,nSteps,jump):
for iSim in range(0,nSim):
    it+=1
    rootFolder = superRootFolder + superDirList[iSim] + "/"
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
    P = np.ma.masked_array(P, P < 0.75*MPa) # Hide the air
    
    mask = P < .75*MPa
    
    thisData = Output.getData(dataFolder + 'strainRate.bin')
    strainRate = thisData.data * 1.0/Char.time
    strainRate = np.ma.masked_array(strainRate, mask) # Hide the air, note: should be done with phase but I forgot to save it
    
    thisData    = Output.getData(dataFolder + 'strain.bin')
    strain      = thisData.data
    strain      = np.ma.masked_array(strain, mask) # Hide the air, note: should be done with phase but I forgot to save it
#    
#    if iSim==0:
#        refStrain = strain
#        
#    diffStrain = strain-refStrain
#    
    State       = Output.readState(dataFolder + "modelState.json")
    timeSim   = (State.time+ State.dt) * Setup.Char.time
    
#    if (np.mod(iStep,100)==0):
#        print("iStep = %i/%i" % (iStep, nSteps))
#        print("ixCell = %i" %ixCell)
    maxStrain_Sim[iSim] = np.max(strain)
       
    print("interp")
    n = 1000
    x = np.linspace(-2250,0,n)
    y = np.linspace(0,2500,n)
    xScaled = (x-xmin)/Wbox * nx
    yScaled = (y-ymin)/Hbox * ny
#    f = interpolate.interp2d(xv - dx/2.0, yv - dy/2.0, strain, kind='linear')
#    strainInterp = interpolate.griddata((x,y), strain, (xv - dx/2.0, yv - dy/2.0), method='nearest')
    strainInterp = ndimage.map_coordinates(strain, [xScaled,yScaled], order=1)
    print("finished interp")
   
    # ColorScales
    # Pressure
    cAx_PMin = 0.0
    cAx_PMax = 50
    # Strain Rate    
    cAx_srMin = -14
    cAx_srMax = -10
    # Strain
    cAx_sMin = 0.0
    cAx_sMax = 3.0
    # diffStrain
    cAx_dsMin = -5.0
    cAx_dsMax =  5.0


#    plt.figure(1)
#    plt.subplot(nSim,1,iSim+1)
#    plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
##    plt.pcolor(xv - dx/2.0,yv - dy/2.0,strain/strainUnit,vmin=cAx_sMin,vmax=cAx_sMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
##    plt.pcolor(xv - dx/2.0,yv - dy/2.0,diffStrain/strainUnit,vmin=cAx_dsMin,vmax=cAx_dsMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#    plt.colorbar()
#    plt.axis('equal')
#    plt.set_cmap("Pressure")
#    
#    plt.plot(x,y,'-m')
    
    plt.pause(0.01)

#    Color = ['b','y','g']
    plt.figure(2)
    plt.subplot(3,1,1)
    plt.plot(strainInterp)

    plt.subplot(3,1,2)
    plt.plot(np.cumsum(np.flipud(strainInterp))*dx/2.0)
    
    # extract topo
    topo = np.zeros(nx)
    for ix in range(nx):
        topo[ix] = ny-np.count_nonzero(mask[ix,:])
    
    
    plt.subplot(3,1,3)
    plt.plot(topo*dy)
    
#plt.figure(2)
#plt.clf()
#plt.plot(np.arange(0,nSim),maxStrain_Sim)
#        
        
#        
#plt.figure(2)
#plt.subplot(2,1,1)
#plt.legend(np.arange(nSim))
#
#plt.subplot(2,1,2)
#plt.legend(np.arange(nSim))
#        
#        
#        
        
        

