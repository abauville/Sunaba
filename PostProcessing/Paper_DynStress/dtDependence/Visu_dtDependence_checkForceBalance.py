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
superRootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_WeakInclusion_NoAdv_NoInterp_adaptative_NEW/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except:
    print("dummy print: no Input")
        
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
jumpy = 31
ny = int(np.ceil(nyTotal/jumpy))
Ix_y = np.zeros(ny);

TauII_t_y = np.zeros([nSteps,ny])
P_t_y = np.zeros([nSteps,ny])
EII_t_y = np.zeros([nSteps,ny])



P_dict      = dict()
EII_dict    = dict()
TauII_dict  = dict()
time_dict   = dict()
NormRes_dict = dict()
dt_dict = dict()
Sxx_dict = dict()
Sxy_dict = dict()

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
    for iSim in range(0,nSim-1):
        
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
        NormRes_dict[superDirList[iSim]]   = np.zeros(nSteps)
        dt_dict[superDirList[iSim]]   = np.zeros(nSteps)
        Sxx_dict[superDirList[iSim]]   = np.zeros(nSteps)
        Sxy_dict[superDirList[iSim]]   = np.zeros(nSteps)
        
        
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
        
        step0 = nSteps-1
        for iStep in range(0,nSteps):
#        for iStep in range(step0,step0+1):
            outFolder = "Out_%05d" % iStep #DirList[iStep]
            print(outFolder)
            
            # index of the first node that register the minimum of khi (where khi is active)
            # Set file
            # =====================
            dataFolder = rootFolder + outFolder + "/"
            
            thisData = Output.getData(dataFolder + 'P.bin')
            
            State       = Output.readState(dataFolder + "modelState.json")
            Char.time = State.Char_time
            Char.length = State.Char_length
            Char.mass = State.Char_mass
            CharExtra = Input.CharExtra(Char)
            P = thisData.data[:,1:-1] * CharExtra.stress
            ixCell = 8 + np.argmin(thisData.data[8:,iyCell]) # find the minimum khi # beyond the inclusion limits, only relevant for the lowest iyCell
            ixCell = 192
            if (np.mod(iStep,100)==0):
                print("iStep = %i/%i" % (iStep, nSteps-1))
                print("ixCell = %i" %ixCell)
            
            
            thisData = Output.getData(dataFolder + 'sigma_xx.bin')
            Sxx = thisData.data[:,1:-1] * CharExtra.stress# I don't need the first and last row (y coordinate)
            thisData = Output.getData(dataFolder + 'sigma_xy_node.bin')
            Sxy = thisData.data * CharExtra.stress
            
            thisData = Output.getData(dataFolder + 'sigma_II.bin')
            SII = thisData.data [:,1:-1] * CharExtra.stress
            
            dSxx_dx = np.diff(Sxx,1,0)/dx
            dSxy_dy = np.diff(Sxy,1,1)/dy
            dSxy_dx = np.diff(Sxy,1,0)/dx
            dP_dx   = np.diff(P  ,1,0)/dx
            
            #Fbx = dSxx_dx + dSxy_dy - dP_dx
            TotalSxx = (P-Sxx - Setup.Physics.Pback)/2.0
            S3 = P-SII
            intTotalSxx_y =  np.sum(TotalSxx,1)*dy#TotalSxx.shape[1]
            intSxx_y = - np.sum(Sxx,1)*dy#TotalSxx.shape[1]
            intS3_y =  np.sum(S3,1)*dy#TotalSxx.shape[1]
            intP_y =  np.sum(P,1)*dy
            intPo_y =  np.sum((P-Setup.Physics.Pback),1)*dy
#            Fbx = dSxy_dy
            #intFbx = np.sum(Fbx,0)*dx
            #intFby = np.sum(Fbx,1)*dy
            
            time_dict[superDirList[iSim]][iStep]  = (State.time + State.dt) * Char.time #
            dt_dict[superDirList[iSim]][iStep]  = (State.dt) * Char.time #
            NormRes_dict[superDirList[iSim]][iStep] = State.residual
            NormRes_dict[superDirList[iSim]][iStep] = P[2,-2]

            # Units and characteristic values
            # =====================
            EII = np.abs(Setup.BC.Stokes.backStrainRate)
            eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
            G = Setup.MatProps['1'].G
            
            t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
            charTime = t
            
            timeUnit = charTime
            stressUnit = Setup.Physics.Pback
            intstress_yUnit = Setup.Physics.Pback*Setup.Grid.nyC*dy
            
            plt.clf()
            plt.plot(intTotalSxx_y/intstress_yUnit)
            plt.plot(intSxx_y/intstress_yUnit)
            #plt.plot(intSxx_y/intPo_y[-1],'.')
            
            plt.plot(intPo_y/intstress_yUnit)
            #plt.plot(intPo_y/intPo_y[-1],'.')
            
            plt.plot(intS3_y/intstress_yUnit)
            #plt.plot(intS3_y/intPo_y[-1],'.')
            
            #np.sum(intPo_y/intPo_y[-1])*dx
            plt.axis([0,nx, 0.5, 1.5 ])
            
            plt.legend(["TotalSxx_y","Sxx_y", "Po_y", "S3_y"])

            plt.pause(0.001)
        
        
        
        
        
        
plt.figure()
plt.plot(NormRes_dict[superDirList[iSim]])
