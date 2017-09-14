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
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
DirList = os.listdir(rootFolder)
DirList.remove('.DS_Store')
DirList.remove('Input')
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
time_t = np.zeros(nSteps)
TauII_t = np.zeros(nSteps)
P_t = np.zeros(nSteps)
EII_t = np.zeros(nSteps)


ny = 100 # to exclude the sticky air
Ix_y = np.zeros(ny);

S1_min_y = np.zeros(ny)
S1_max_y = np.zeros(ny)
S3_min_y = np.zeros(ny)
S3_max_y = np.zeros(ny)
TauII_min_y = np.zeros(ny)
TauII_max_y = np.zeros(ny)
P_min_y = np.zeros(ny)
P_max_y = np.zeros(ny)
EII_min_y = np.zeros(ny)
EII_max_y = np.zeros(ny)

time_S1_min_y = np.zeros(ny)
time_S1_max_y = np.zeros(ny)
time_S3_min_y = np.zeros(ny)
time_S3_max_y = np.zeros(ny)
time_TauII_min_y = np.zeros(ny)
time_TauII_max_y = np.zeros(ny)
time_P_min_y = np.zeros(ny)
time_P_max_y = np.zeros(ny)
time_EII_min_y = np.zeros(ny)
time_EII_max_y = np.zeros(ny)

Compute=True

if Compute:
    for iyCell in range(0,ny):
        print("iyCell = %i" % iyCell)
        Ix_t = np.zeros(nSteps)
        # First pass: Find which cell to monitor
        # ===================================
        iStep = 0
        for outFolder in DirList:
            # Set file
            # =====================
            #outFolder = "Out_00009"
            dataFolder = rootFolder + outFolder + "/"
            
            # Read data
            # =====================
            thisData = Output.getData(dataFolder + 'khi.bin')
            #iyCell = 85#thisData.ny-2
            Ix = int(np.argmin(thisData.data[:,iyCell]) )
            Ix_t[iStep] = Ix
            iStep += 1;
        
        Ix_y[iyCell] = Ix_t[int(np.argmax(Ix_t>0))]
        
        # Second pass: Get the time sequence of P, TauII, EII
        # ===================================
        iStep = 0
        for outFolder in DirList:
             # index of the first node that register the minimum of khi (where khi is active)
            # Set file
            # =====================
            #outFolder = "Out_00009"
            dataFolder = rootFolder + outFolder + "/"
            
            # Read data
            # =====================
            state = Output.readState(dataFolder + 'modelState.json')
            # Write data
            # =====================
            time_t[iStep] = state.time
            thisData = Output.getData(dataFolder + 'P.bin')
            P_t[iStep] = thisData.data[Ix,iyCell] * CharExtra.stress
        
            thisData = Output.getData(dataFolder + 'Sigma_II.bin')
            TauII_t[iStep] = thisData.data[Ix,iyCell] * CharExtra.stress
        
            thisData = Output.getData(dataFolder + 'strainRate.bin')
            EII_t[iStep] = thisData.data[Ix,iyCell] * 1.0/Char.time
        
            iStep += 1;
    
    
        # Transform the data to get S1, S3, Sy as well
        # ===================================
        S1_t = P_t + TauII_t
        S3_t = P_t - TauII_t
        Sy_t = C*cos(phi) + P_t*sin(phi)
    
    
        # Extract the maximum, minimum and respective times
        # ===================================
        
        Imax = int(np.argmax(S1_t))
#        Imin = Imax + int(np.argmin(S1_t[Imax:]))
        Imin = Imax + int(np.argmin(S3_t[Imax:]))
        S1_max_y[iyCell]    = S1_t[Imax]
        S1_min_y[iyCell]    = S1_t[Imin]
        time_S1_max_y[iyCell]    = time_t[Imax]
        time_S1_min_y[iyCell]    = time_t[Imin]
        
#        Imax = int(np.argmax(S3_t))
#        Imin = Imax + int(np.argmin(S3_t[Imax:]))
        S3_max_y[iyCell]    = S3_t[Imax]
        S3_min_y[iyCell]    = S3_t[Imin]
        time_S3_max_y[iyCell]    = time_t[Imax]
        time_S3_min_y[iyCell]    = time_t[Imin]
        
#        Imax = int(np.argmax(TauII_t))
#        Imin = Imax + int(np.argmin(TauII_t[Imax:]))
        TauII_max_y[iyCell]    = TauII_t[Imax]
        TauII_min_y[iyCell]    = TauII_t[Imin]
        time_TauII_max_y[iyCell]    = time_t[Imax]
        time_TauII_min_y[iyCell]    = time_t[Imin]
        
#        Imax = int(np.argmax(P_t))
#        Imin = Imax + int(np.argmin(P_t[Imax:]))
        P_max_y[iyCell]    = P_t[Imax]
        P_min_y[iyCell]    = P_t[Imin]
        time_P_max_y[iyCell]    = time_t[Imax]
        time_P_min_y[iyCell]    = time_t[Imin]
        
#        Imax = int(np.argmax(EII_t))
#        Imin = Imax + int(np.argmin(EII_t[Imax:]))
        EII_max_y[iyCell]    = EII_t[Imax]
        EII_min_y[iyCell]    = EII_t[Imin]
        time_EII_max_y[iyCell]    = time_t[Imax]
        time_EII_min_y[iyCell]    = time_t[Imin]

        

        np.savez("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/Stress_vs_y",
                 Ix_y = Ix_y,
                 S1_min_y = S1_min_y,
                 S1_max_y = S1_max_y,
                 S3_min_y = S3_min_y,
                 S3_max_y = S3_max_y,
                 TauII_min_y = TauII_min_y,
                 TauII_max_y = TauII_max_y,
                 P_min_y = P_min_y,
                 P_max_y = P_max_y,
                 EII_min_y = EII_min_y,
                 EII_max_y = EII_max_y,
                 
                 time_S1_min_y = time_S1_min_y,
                 time_S1_max_y = time_S1_max_y,
                 time_S3_min_y = time_S3_min_y,
                 time_S3_max_y = time_S3_max_y,
                 time_TauII_min_y = time_TauII_min_y,
                 time_TauII_max_y = time_TauII_max_y,
                 time_P_min_y = time_P_min_y,
                 time_P_max_y = time_P_max_y,
                 time_EII_min_y = time_EII_min_y,
                 time_EII_max_y = time_EII_max_y
                 )
             
else:
    loadedData = np.load("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/Stress_vs_y.npz");
    Ix_y = loadedData["Ix_y"]
    S1_min_y = loadedData["S1_min_y"]
    S1_max_y = loadedData["S1_max_y"]
    S3_min_y = loadedData["S3_min_y"]
    S3_max_y = loadedData["S3_max_y"]
    TauII_min_y = loadedData["TauII_min_y"]
    TauII_max_y = loadedData["TauII_max_y"]
    P_min_y = loadedData["P_min_y"]
    P_max_y = loadedData["P_max_y"]
    EII_min_y = loadedData["EII_min_y"]
    EII_max_y = loadedData["EII_max_y"]
    
    time_S1_min_y = loadedData["time_S1_min_y"]
    time_S1_max_y = loadedData["time_S1_max_y"]
    time_S3_min_y = loadedData["time_S3_min_y"]
    time_S3_max_y = loadedData["time_S3_max_y"]
    time_TauII_min_y = loadedData["time_TauII_min_y"]
    time_TauII_max_y = loadedData["time_TauII_max_y"]
    time_P_min_y = loadedData["time_P_min_y"]
    time_P_max_y = loadedData["time_P_max_y"]
    time_EII_min_y = loadedData["time_EII_min_y"]
    time_EII_max_y = loadedData["time_EII_max_y"]
    
#endif
    



# Analytical yield stress
# =====================  
S3 = Setup.Physics.Pback
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
Sy_back = (S1-S3)/2.0
P_Lim = (S1+S3)/2.0

# Units and characteristic values
# =====================
EII = np.abs(Setup.BC.Stokes.backStrainRate)
eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
G = Setup.MatProps['1'].G

t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
charTime = t

timeUnit = charTime
stressUnit = Setup.Physics.Pback



# Plot
# =====================
Y = -dy/2 + dy*np.arange(0,ny)

# min, max Stress vs y
plt.figure(1)
plt.clf()
plt.plot(Y,S1_max_y/stressUnit,'.g')
plt.plot(Y,S1_min_y/stressUnit,'xg')
plt.plot(Y[[0,-1]],np.array([S1,S1])/stressUnit,'--g')

plt.plot(Y,S3_max_y/stressUnit,'.y')
plt.plot(Y,S3_min_y/stressUnit,'xy')
plt.plot(Y[[0,-1]],np.array([S3,S3])/stressUnit,'--y')

plt.plot(Y,P_max_y/stressUnit,'.b')
plt.plot(Y,P_min_y/stressUnit,'xb')
plt.plot(Y[[0,-1]],np.array([P_Lim,P_Lim])/stressUnit,'--b')

plt.plot(Y,TauII_max_y/stressUnit,'.r')
plt.plot(Y,TauII_min_y/stressUnit,'xr')
plt.plot(Y[[0,-1]],np.array([Sy_back,Sy_back])/stressUnit,'--r')

# time of min, max stress vs y
plt.figure(2)
plt.clf()
plt.plot(Y,time_TauII_max_y/timeUnit,'.r')
plt.plot(Y,time_TauII_min_y/timeUnit,'xr')

plt.plot(Y,time_S1_max_y/timeUnit,'.g')
plt.plot(Y,time_S1_min_y/timeUnit,'xg')

plt.plot(Y,time_S3_max_y/timeUnit,'.y')
plt.plot(Y,time_S3_min_y/timeUnit,'xy')

plt.plot(Y,time_P_max_y/timeUnit,'.b')
plt.plot(Y,time_P_min_y/timeUnit,'xb')





## Mohr Diagram
# =====================
fig = plt.figure(3)
plt.clf()
ax = fig.add_subplot(1, 1, 1)

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')

# Eliminate upper and right axes
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')




iy = 50
P = P_max_y[iy] / stressUnit
TauII = TauII_max_y[iy] / stressUnit
nPoints = 50
psi = np.linspace(0,1*np.pi,nPoints)

# Coulomb enveloppe
xEnveloppe = np.linspace(0.0,1.2*(P+TauII),2)
Sy = C/stressUnit + tan(phi)*xEnveloppe
plt.plot(xEnveloppe,Sy,'-k')

# Mohr Circle
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='r')

P = P_min_y[iy] / stressUnit
TauII = TauII_min_y[iy] / stressUnit
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='y')


##plt.plot(Sigma1,0,"xr")
##plt.plot(Sigma3,0,"xb")
##plt.plot(P+Tauxx, Tauxy, "xg")
##plt.plot(P+Tauyy, Tauyx, "xm")
#
#plt.plot([P,P+TauII*np.cos(2*psi)],[0,TauII*np.sin(2*psi)])
#
plt.axis("equal")
#
#plt.show()
#
#
#

