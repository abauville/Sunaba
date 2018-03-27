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




cdict1 = {'red':  ((0.0 , 1.0, 1.0),
                   (0.25, 0.25, 0.25),
                   (0.5 , 1.0, 1.0),
                   (0.75, 1.0, 1.0),
                   (1.0 , 0.0, 0.0)),

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


superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNewInvType1/FixedDt_Method1/Output/"
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



nSteps = len(DirList)
jump = 20
nSteps = int(nSteps/jump)
nSteps = 240
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

#plt.figure(1)
# plot
#ResFac = ResFacList[iSim]
#iyCell = int(100/2*ResFac)+1
#plt.clf()
#plt.ion()
i0 = 0
#jump = 1
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
    

#iStep = 10
nSim = len(superDirList)
plt.figure(3)
plt.clf()
#plt.figure(2)
#plt.clf()

iSim = 0
#plt.figure(2)
#plt.clf()
for iStep in range(i0,nSteps,jump):
#for iSim in range(0,nSim):
    it+=1
    rootFolder = superRootFolder + superDirList[iSim] + "/"
    outFolder = "Out_%05d" % iStep #DirList[iStep]
    print("iStep = %i/%i" % (iStep, nSteps-1))
    # index of the first node that register the minimum of khi (where khi is active)
    # Set file
    # =====================
    
    
    #Setup = Output.readInput(rootFolder +  'Input/input.json')
    #Char = Setup.Char
    #CharExtra = Input.CharExtra(Char)

    dataFolder = rootFolder + outFolder + "/"
    State       = Output.readState(dataFolder + "modelState.json")
    timeSim   = (State.time+ State.dt) * Setup.Char.time
    
    
    phase = Output.getData(dataFolder + 'phase.bin').data
    mask = phase == 0
    
    strainRate = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
#    strain = Output.getData(dataFolder + 'strain.bin',True,mask).data
    Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True,mask).data
    Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True,mask).data
    SII = Output.getData(dataFolder + 'sigma_II.bin',True,mask).data
    Pressure   = Output.getData(dataFolder + 'P.bin',True,mask).data
#    strainRate = thisData.data * 1.0/Char.time
#    strainRate = np.ma.masked_array(strainRate, mask) # Hide the air, note: should be done with phase but I forgot to save it



#    
    
#    maxStrain_Sim[iSim] = np.max(strain)
#       
#    print("interp")
#    n = 1000
#    x = np.linspace(-2250,0,n)
#    y = np.linspace(0,2500,n)
#    xScaled = (x-xmin)/Wbox * nx
#    yScaled = (y-ymin)/Hbox * ny
##    f = interpolate.interp2d(xv - dx/2.0, yv - dy/2.0, strain, kind='linear')
##    strainInterp = interpolate.griddata((x,y), strain, (xv - dx/2.0, yv - dy/2.0), method='nearest')
#    strainInterp = ndimage.map_coordinates(strain, [xScaled,yScaled], order=1)
#    print("finished interp")
#   
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


    
    
#    
#    # Extract Topo
    yvMasked = yv.copy()#np.ma.masked_array(yv, mask)
    yvMasked[mask] = 0.0
    yvMasked = yvMasked 
#    Topo = np.max(yvMasked,1)
    Topo_iy = np.argmax(yvMasked,1)
    Topo_iy -= 2# to be sure that I don't sample the air, because phase and Vx don't have the same number of points
#    
#
#    
    ixS = np.argmin(np.abs(xv[:,0]-(-2000.0)))
#    iyE = 120


    plt.subplot((round(nSteps/jump)+1)/4,4,it+1)
    plt.pcolor(xv[ixS:,:] - dx/2.0,yv[ixS:,:] - dy/2.0,np.log10(strainRate[ixS:,:]/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
##    plt.pcolor(xv - dx/2.0,yv - dy/2.0,strain/strainUnit,vmin=cAx_sMin,vmax=cAx_sMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
##    plt.pcolor(xv - dx/2.0,yv - dy/2.0,diffStrain/strainUnit,vmin=cAx_dsMin,vmax=cAx_dsMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    
    
    plt.axis('equal')
    plt.axis([-2000,0,0,1100])
#    plt.set_cmap("wcyrk")
#    plt.set_cmap("Pressure")
    plt.set_cmap("gist_gray")
    plt.axis("off")
#    plt.title("time=%.2e yr" % (timeSim/yr))
#    plt.colorbar()
#    plt.tight_layout()
    
    refP = np.mean(Pressure)*10.0
    
    Segment = np.array([[-1,1],[0,0]]) * Hbox / 10.0
    step = 10
    for iy in range(0,ny,step):
        for ix in range(ixS,nx,step):
            iCell =ix + nx*iy
            thisSxy = Sxy[ix,iy]
            thisSxx = Sxx[ix,iy]
            if thisSxy!=0:
                Tau = thisSxx/thisSxy
            else:
                Tau = 10000000000.0;
           
            if thisSxy<0.0:
                psi = arctan(-Tau+np.sqrt(Tau*Tau+1))
            else:
                psi = arctan(-Tau-np.sqrt(Tau*Tau+1));

            x = xmin + dx*ix
            y = ymin + dy*iy
            
            # Plot the first principal stress direction
            rot = np.array([[cos(psi), -sin(psi)],[sin(psi), cos(psi)]])
#            scale = (Pressure[ix,iy]+SII[ix,iy])/refP
            scale = +SII[ix,iy]/refP
            rotSegment = np.matmul(rot,Segment)
            
            SglyphOptions = dict(linestyle='-',color=[1,1,0],linewidth=1.5)
            plt.plot(x+rotSegment[0,:]*scale,y+rotSegment[1,:]*scale,**SglyphOptions)
#    
#            # Plot the thrid principal stress direction
#            psi += np.pi/2.0
#            rot = np.array([[cos(psi), -sin(psi)],[sin(psi), cos(psi)]])
##            scale = (Pressure[ix,iy]-SII[ix,iy])/refP
#            scale = (-SII[ix,iy])/refP
#            rotSegment = np.matmul(rot,Segment)
#            plt.plot(x+rotSegment[0,:]*scale,y+rotSegment[1,:]*scale,**SglyphOptions)
    
#    
#    
            
#    plt.axis('equal')
            
    plt.pause(0.01)



plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/FaultInitiationSnapshots.png",r=500)


