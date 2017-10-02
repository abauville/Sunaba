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


#rootFolder = "/Users/abauville/StokesFD_Output/ViscoElasticBuildUp/"
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_NoAdv_NoInterp_adaptative/"
superDirList = os.listdir(superRootFolder)
superDirList.remove('.DS_Store')
nSim = len(superDirList)
#superDirList.remove('Input')

dt_stressFacList = np.zeros(nSim)
for i in range(nSim):
    dt_stressFacList[i] = float(superDirList[i][-7:])

I = dt_stressFacList.argsort()
I.astype(np.int)
I = I[::-1]
dt_stressFacList = dt_stressFacList[I]
superDirList = [ superDirList[i] for i in I]

rootFolder = superRootFolder + superDirList[0] + "/"

DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")




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


nyTotal = 100 # to exclude the sticky air
jumpy = 51
ny = int(np.ceil(nyTotal/jumpy))
Ix_y = np.zeros(ny);

TauII_t_y = np.zeros([nSteps,ny])
P_t_y = np.zeros([nSteps,ny])
EII_t_y = np.zeros([nSteps,ny])



P_dict      = dict()
EII_dict    = dict()
TauII_dict  = dict()
time_dict   = dict()

# Extract Data from file
# =======================
# Get Ix for each Sim

# Test
#nSim = 3
#nSim = 9


ExtractData=True
if ExtractData:
    iyCell = 51
    ixCell = 85
    for iSim in range(nSim):
        
#        dt_stressFac = dt_stressFacList[iSim]
        rootFolder = superRootFolder + superDirList[iSim] + "/"
        iStep = 0
        
        # Get the list of directories for this sim
        DirList = os.listdir(rootFolder)
        try:
            DirList.remove('.DS_Store')
        except ValueError:
            print("dummy print: no .DS_Store")
            
        try:
            DirList.remove('Input')
        except ValueError:
            print("dummy print: no Input")
        
        
        nSteps = len(DirList)
        P_dict[superDirList[iSim]]      = np.zeros(nSteps)
        TauII_dict[superDirList[iSim]]  = np.zeros(nSteps)
        EII_dict[superDirList[iSim]]    = np.zeros(nSteps)
        time_dict[superDirList[iSim]]   = np.zeros(nSteps)
        
        
        Setup = Output.readInput(rootFolder +  'Input/input.json')
        Char = Setup.Char
        CharExtra = Input.CharExtra(Char)
        
        
        
#        # First pass: Find which cell to monitor
#        # ===================================
#        iStep = 0
#        Ix_t = np.zeros(nSteps)
#        jump = int(nSteps/100)
#        print("iSim = %i, nSteps = %i, jump = %i" % (iSim,nSteps, jump))
#        for i in range(0,nSteps,jump):
#            outFolder = DirList[i]
#            # Set file
#            # =====================
#            #outFolder = "Out_00009"
#            dataFolder = rootFolder + outFolder + "/"
#            
#            # Read data
#            # =====================
#            thisData = Output.getData(dataFolder + 'khi.bin')
#            #iyCell = 85#thisData.ny-2
#            Ix = int(np.argmin(thisData.data[:,iyCell]) )
#            Ix_t[iStep] = Ix
#            if (Ix>0):
#                break
#            iStep += 1;
#        
#        ixCell = Ix#int(Ix_t[int(np.argmax(Ix_t>0))])
#        print("ixCell = %i" %ixCell)
        
        
        for iStep in range(0,nSteps):
            outFolder = "Out_%05d" % iStep #DirList[iStep]
            
            
            # index of the first node that register the minimum of khi (where khi is active)
            # Set file
            # =====================
            dataFolder = rootFolder + outFolder + "/"
            
            thisData = Output.getData(dataFolder + 'P.bin')
            ixCell = np.argmin(thisData.data[:,iyCell]) # find the minimum khi
            if (np.mod(iStep,100)==0):
                print("iStep = %i/%i" % (iStep, nSteps-1))
                print("ixCell = %i" %ixCell)
            
            State       = Output.readState(dataFolder + "modelState.json")
            
            
            P_dict[superDirList[iSim]][iStep]     =  Output.getData_OneCell(dataFolder + 'P.bin', ixCell,iyCell) * CharExtra.stress
            TauII_dict[superDirList[iSim]][iStep] =  Output.getData_OneCell(dataFolder + 'sigma_II.bin', ixCell,iyCell) * CharExtra.stress
            EII_dict[superDirList[iSim]][iStep]   =  Output.getData_OneCell(dataFolder + 'strainRate.bin', ixCell,iyCell) * 1.0/Char.time
            time_dict[superDirList[iSim]][iStep]  = (State.time+ State.dt) * Setup.Char.time #
        #end iStep
    #end iSim
    np.savez("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/dtDependenceAdaptative_minP",
             P_dict = P_dict,
             TauII_dict = TauII_dict,
             EII_dict = EII_dict,
             time_dict = time_dict
             )
else:
    loadedData = np.load("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/dtDependenceAdaptative_minP.npz");
    P_dict     = loadedData["P_dict"][()]
    TauII_dict = loadedData["TauII_dict"][()]
    EII_dict   = loadedData["EII_dict"][()]
    time_dict  = loadedData["time_dict"][()]
        
# Analytical yield stress
# =====================  
S3 = Setup.Physics.Pback
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
Sy_back = (S1-S3)/2.0
P_lim = (S1+S3)/2.0

# Units and characteristic values
# =====================
EII = np.abs(Setup.BC.Stokes.backStrainRate)
eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
G = Setup.MatProps['1'].G

t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
charTime = t

timeUnit = charTime
stressUnit = Setup.Physics.Pback

#plt.close("all")
plt.figure(5)
plt.clf()
iSim0 = 0
#for iSim in range(iSim0,nSim):
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, P_dict[superDirList[iSim]]/stressUnit,'.')
for iSim in range(iSim0,nSim):
#for iSim in range(nSim-1,iSim0-1,-1):
    plt.plot(time_dict[superDirList[iSim]]/timeUnit+dt_stressFacList[iSim], TauII_dict[superDirList[iSim]]/stressUnit,'.',markersize=1)

#for iSim in range(nSim-1,nSim-0):
#    P = P_dict[superDirList[iSim]]
#    TauII = TauII_dict[superDirList[iSim]]
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, P/stressUnit,'.')
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, TauII/stressUnit,'.')
#    Tau_y = C*cos(phi) + P*sin(phi)
#    S1 = P+TauII
#    S3 = P-TauII
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, Tau_y/stressUnit,'-')
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, S1/stressUnit,'-')
#    plt.plot(time_dict[superDirList[iSim]]/timeUnit, S3/stressUnit,'-')
##    

plt.plot([0,2],np.array([Sy_back, Sy_back])/stressUnit)



#        
plt.plot([0,2],np.array([P_lim, P_lim])/stressUnit)
plt.axis([0,2.5,0,4.0])
plt.legend(superDirList[iSim0:])
#plt.plot(time_dict[superDirList[iSim]]/timeUnit, TauII_dict[superDirList[iSim]]/stressUnit,'-')
#plt.plot(time_dict[superDirList[9]][0::100]/timeUnit, TauII_dict[superDirList[9]][0::100]/stressUnit,'-')       

#plt.plot(np.diff(time_dict[superDirList[9]]))

        
        
        
        
        
        
        
        
        
        
        

