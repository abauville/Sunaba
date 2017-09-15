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
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_Stronger_Seed_10timesWeaker/dt_stressFac_1.0e-03/"

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


TauII_t_y = np.zeros([nSteps,ny])
P_t_y = np.zeros([nSteps,ny])
EII_t_y = np.zeros([nSteps,ny])


ExtractData=False



# Extract Data from file
# =======================
# Data is the time series of a point on the fault at each position along y
if ExtractData:
    iy = -1
    for iyCell in range(0,nyTotal,jumpy):
        iy+=1
        print("iyCell = %i" % iyCell)
        print("0")
        Ix_t = np.zeros(nSteps)
        print("A")
        t = time.time()
        # First pass: Find which cell to monitor
        # ===================================
        iStep = 0
#        for outFolder in DirList:
        for i in range(0,nSteps,jump):
            outFolder = DirList[i]
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
        
#        Ix_y[iyCell] = Ix_t[int(np.argmax(Ix_t>0))]
        Ix = int(Ix_t[int(np.argmax(Ix_t>0))])
        elapsed = time.time() - t
        print("B,  %.1f s" % elapsed)
        # Second pass: Get the time sequence of P, TauII, EII
        # ===================================
        iStep = 0
        
        t = time.time()


#        for outFolder in DirList:
        for i in range(0,nSteps,jump):
            outFolder = DirList[i]
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
            
            P_t_y[iStep,iy] = Output.getData_OneCell(dataFolder + 'P.bin', Ix,iyCell) * CharExtra.stress
            TauII_t_y[iStep,iy] = Output.getData_OneCell(dataFolder + 'Sigma_II.bin',  Ix,iyCell) * CharExtra.stress
            EII_t_y[iStep,iy] = Output.getData_OneCell(dataFolder + 'strainRate.bin',  Ix,iyCell) * 1.0/Char.time

        
            iStep += 1;
            
        ###############
        
        elapsed = time.time() - t
        print("C, %.1f s" % elapsed)
        t = time.time()
        
    np.savez("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/Stress_vs_yc",
             Ix_y = Ix_y,
             P_t_y = P_t_y,
             TauII_t_y = TauII_t_y,
             EII_t_y = EII_t_y
             )
         
else:
    loadedData = np.load("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/Stress_vs_yc.npz");
    Ix_y        = loadedData["Ix_y"]
    P_t_y       = loadedData["P_t_y"]
    TauII_t_y   = loadedData["TauII_t_y"]
    EII_t_y     = loadedData["EII_t_y"]
    
    
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



# Compute min, max etc...
# =====================
Compute = True
if Compute:
    S1_t_y = P_t_y + TauII_t_y
    S3_t_y = P_t_y - TauII_t_y
    Sy_t_y = C*cos(phi) + P_t_y*sin(phi)
    
    Imax_t = np.argmax(S1_t_y,0)
#    Imax_t = np.argmin(np.abs(Sy_t_y-TauII_t_y),0)
    
    Ihalf_t = Imax_t/2
    Iquar_t = Imax_t/4
    Ihalf_t = Ihalf_t.astype(Imax_t.dtype)
    Iquar_t = Iquar_t.astype(Imax_t.dtype)
    Imin_t = np.zeros(ny, dtype=np.int)
    for i in range(S1_t_y.shape[1]):    
        Imin_t[i] = Imax_t[i] + np.argmin(S3_t_y[Imax_t[i]:,i],0)
#        Imin_t[i] = nSteps-1
        
        
    Imax = np.ravel_multi_index([Imax_t,np.arange(ny)],S1_t_y.shape)
    Iquar = np.ravel_multi_index([Iquar_t,np.arange(ny)],S1_t_y.shape)
    Ihalf = np.ravel_multi_index([Ihalf_t,np.arange(ny)],S1_t_y.shape)
    Imin = np.ravel_multi_index([Imin_t,np.arange(ny)],S1_t_y.shape)
    
    
    
    S1max_y   = S1_t_y.flat[Imax]         / stressUnit
    S1half_y  = S1_t_y.flat[Ihalf]        / stressUnit
    S1min_y   = S1_t_y.flat[Imin]         / stressUnit
    
    S3max_y   = S3_t_y.flat[Imax]         / stressUnit
    S3half_y  = S3_t_y.flat[Ihalf]        / stressUnit
    S3min_y   = S3_t_y.flat[Imin]         / stressUnit
    
    Pmax_y   = P_t_y.flat[Imax]           / stressUnit
    Pquar_y  = P_t_y.flat[Iquar]          / stressUnit
    Phalf_y  = P_t_y.flat[Ihalf]          / stressUnit
    Pmin_y   = P_t_y.flat[Imin]           / stressUnit
    
   
    TauIImax_y   = TauII_t_y.flat[Imax]   / stressUnit
    TauIIquar_y  = TauII_t_y.flat[Iquar]  / stressUnit
    TauIIhalf_y  = TauII_t_y.flat[Ihalf]  / stressUnit
    TauIImin_y   = TauII_t_y.flat[Imin]   / stressUnit
    
    

    
    
#    Imax = np.ravel_multi_index([Imax_t,np.arange(ny)],S1_t_y.shape)
    
    
#    S1 = np.amax(S1_t_y,0)
#    S1max = S1_t_y.flat[Imax]
#    S1b = S1_t_y.flat[Imax]
#    Ibuild_t = Imax_t/2
#        Imin = Imax + int(np.argmin(S1_t[Imax:]))
#    Imin_t = Imax_t + int(np.argmin(S3_t_y[Imax_t:,:]),0)



# Plot
# =====================
Y = -dy/2 + dy*np.arange(0,ny)

#plt.close("all")
plt.figure(1)
plt.clf()
plt.plot(Y,S1max_y,'.g')
plt.plot(Y,S1half_y,'xg')
plt.plot(Y,S1min_y,'sg')

plt.plot(Y,S3max_y,'.y')
plt.plot(Y,S3half_y,'xy')
plt.plot(Y,S3min_y,'sy')

plt.plot(Y,Pmax_y,'.b')
plt.plot(Y,Phalf_y,'xb')
plt.plot(Y,Pmin_y,'sb')

plt.plot(Y,TauIImax_y,'.r')
plt.plot(Y,TauIIhalf_y,'xr')
plt.plot(Y,TauIImin_y,'sr')


#plt.plot(Y,S1b,'xr')










## Mohr Diagram
# =====================
fig = plt.figure(2)
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




#iy = int(ny/2)
iy = 1
C = C/stressUnit

nPoints = 50
psi = np.linspace(0,1*np.pi,nPoints)

# Coulomb enveloppe
xEnveloppe = np.linspace(0.0,4.0,2)
Sy = C + tan(phi)*xEnveloppe
plt.plot(xEnveloppe,Sy,'-k')

# Mohr Circle





P = Pquar_y[iy]
TauII = TauIIquar_y[iy]
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='y')
plt.plot(P,0.0,'.y')

P = Phalf_y[iy]
TauII = TauIIhalf_y[iy]
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='b')
plt.plot(P,0.0,'.b')

P = Pmax_y[iy] 
TauII = TauIImax_y[iy]
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='g')
plt.plot(P,0.0,'.g')

A = TauII*sin(phi)
B = TauII*sin(phi)*cos(phi)

plt.plot([P-A,P-A],[0.0,3.0],'--k')
plt.plot([P-A-B,P-A-B],[0.0,3.0],'--k')

P = Pmin_y[iy]
TauII = TauIImin_y[iy]
plt.plot(P+TauII*np.cos(psi), TauII*np.sin(psi), linewidth=1.5, color='r')
plt.plot(P,0.0,'.r')

plt.axis("equal")
fig = plt.figure(3)
plt.clf()
plt.subplot(2,1,1)
plt.plot(P_t_y[:,iy]/stressUnit,'-b')
plt.plot(S3_t_y[:,iy]/stressUnit,'-y')
plt.subplot(2,1,2)
plt.plot(np.diff(P_t_y[:,iy]/stressUnit),'-b')
plt.plot(np.diff(S3_t_y[:,iy]/stressUnit),'-y')

#plt.plot((P_t_y[:,iy]-S3_t_y[:,iy])/stressUnit,'-k')



#plt.axis("equal")


