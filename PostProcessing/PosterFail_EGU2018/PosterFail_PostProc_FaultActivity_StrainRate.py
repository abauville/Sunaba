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


superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNew/FixedDt_Method1/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNew/FixedDt_Method0/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/Test00/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/Corotational/FixedDt_Method0/Output/"
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
iSim = 1
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


#
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    os.system("mkdir /home/abauvill/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
#elif sys.platform == "darwin":
#    # OS X
#    os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
    

iStep = 10
nSim = len(superDirList)
plt.figure(1)
plt.clf()
#plt.figure(2)
#plt.clf()

i0 = 00
#nSteps = 500
jump= 1

nt = len(range(i0,nSteps,jump))
maxStrainRate_sub0 = np.zeros(nt)
maxStrainRate_sub1 = np.zeros(nt)
maxStrainRate_sub2 = np.zeros(nt)
avStrainRate_sub0 = np.zeros(nt)
avStrainRate_sub1 = np.zeros(nt)
avStrainRate_sub2 = np.zeros(nt)
posMaxStrainRate_sub0 = np.zeros(nt)


Compute = False

if Compute:
    
    # Extract a subset of the data
    # coord of the subset0
    x0 = xmin+Wbox*0.3
    x1 = xmin+Wbox*0.8
    y0 = ymin+Hbox*0.3
    y1 = 900#ymin+Hbox*0.325
    #ix0_sub0 = np.argmin(np.abs(x-x0))
    #ix1_sub0 = np.argmin(np.abs(x-x1))
    ix0_sub0  = 0
    ix1_sub0  = nx
    iy0_sub0 = np.argmin(np.abs(y-y0))
    iy1_sub0 = iy0_sub0+1
    #    iy1_sub0 = np.argmin(np.abs(y-y1))
    
    
    
    #    
    # coord of the subset1    
    x0 = xmin+Wbox*0.6
    x1 = xmin+Wbox*0.9
    y0 = ymin+Hbox*-0.05
    y1 = ymin+Hbox*0.025
    #    ix0_sub1 = np.argmin(np.abs(x-x0))
    #    ix1_sub1 = np.argmin(np.abs(x-x1))
    ix0_sub1 = ix0_sub0
    ix1_sub1 = ix1_sub0
    iy0_sub1 = np.argmin(np.abs(y-y0))
    #    iy1_sub1 = np.argmin(np.abs(y-y1))
    iy1_sub1 = iy0_sub1+1
    
    
    ## coord of the subset12 
    #x0 = xmin+Wbox*0.6
    #x1 = xmin+Wbox*1.0
    #y0 = ymin+Hbox*-0.05
    #y1 = ymin+Hbox*0.9
    #ix0_sub2 = np.argmin(np.abs(x-x0))
    #ix1_sub2 = np.argmin(np.abs(x-x1))
    #iy0_sub2 = np.argmin(np.abs(y-y0))
    #iy1_sub2 = np.argmin(np.abs(y-y1))
    
    
    
    
    
    strainRate_profile_sub0 = np.zeros((ix1_sub0-ix0_sub0,nt))
    strainRate_profile_sub1 = np.zeros((ix1_sub0-ix0_sub0,nt))
    #VxTopo = np.zeros((nx-1,nt))
    
    
    
    
    
    
    
    
    
    

    it=-1
    
    for iStep in range(i0,nSteps,jump):
    #for iSim in range(0,nSim):
    
        it+=1
        rootFolder = superRootFolder + superDirList[iSim] + "/"
        outFolder = "Out_%05d" % iStep #DirList[iStep]
        print("iStep = %i/%i" % (iStep, nSteps-1))
        # index of the first node that register the minimum of khi (where khi is active)
        # Set file
        # =====================
        
        
        #Setup = Output.readInput(rootFolder +  'Input/input.json')sa
        #Char = Setup.Char
        #CharExtra = Input.CharExtra(Char)
    
        dataFolder = rootFolder + outFolder + "/"
        State       = Output.readState(dataFolder + "modelState.json")
        timeSim   = (State.time+ State.dt) * Setup.Char.time
        
        phase = Output.getData(dataFolder + 'phase.bin').data
        mask = phase == 0
        
        strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
        
        
#        Vx = Output.getData(dataFolder + 'Vx.bin',True).data
        
        
        # Extract Topo
    #    yvMasked = yv#np.ma.masked_array(yv, mask)
    #    yvMasked[mask] = 0.0
    #    yvMasked = yvMasked 
    ##    Topo = np.max(yvMasked,1)
    #    Topo_iy = np.argmax(yvMasked,1)
    #    Topo_iy -= 2# to be sure that I don't sample the air, because phase and Vx don't have the same number of points
    #    
    #    for ix in range(nx-1):
    #       VxTopo[ix,it] = Vx[ix,Topo_iy[ix]] 
    #    
    #    
        
    #    strain = Output.getData(dataFolder + 'strain.bin',True,mask).data
    #    Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True,mask).data
    #    Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True,mask).data
    #    SII = Output.getData(dataFolder + 'sigma_II.bin',True,mask).data
    #    Pressure   = Output.getData(dataFolder + 'P.bin',True,mask).data
    #    strainRate = thisData.data * 1.0/Char.time
    #    strainRate = np.ma.masked_array(strainRate, mask) # Hide the air, note: should be done with phase but I forgot to save it
    
    
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
    #    # Plot image
    #    plt.figure(1)
    #    plt.clf()
    ##    plt.subplot(nSim,1,iSim+1)
    #    plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    ###    plt.pcolor(xv - dx/2.0,yv - dy/2.0,strain/strainUnit,vmin=cAx_sMin,vmax=cAx_sMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    ###    plt.pcolor(xv - dx/2.0,yv - dy/2.0,diffStrain/strainUnit,vmin=cAx_dsMin,vmax=cAx_dsMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    #    plt.colorbar()
    #    plt.axis('equal')
    #    plt.set_cmap("Pressure")
    #    plt.title("tstep = %i" % (it))
    #    
    ##    refP = np.mean(Pressure)*10.0
    #    
        
    #    
    #    plt.plot(x[[ix0_sub0, ix1_sub0, ix1_sub0, ix0_sub0,ix0_sub0]],y[[iy0_sub0,iy0_sub0,iy1_sub0,iy1_sub0,iy0_sub0]])
    #    plt.plot(x[[ix0_sub1, ix1_sub1, ix1_sub1, ix0_sub1,ix0_sub1]],y[[iy0_sub1,iy0_sub1,iy1_sub1,iy1_sub1,iy0_sub1]])
    #    plt.pause(0.01)
        
        
        # Extract fault activation information
        
        maxStrainRate_sub0[it] = np.max(strainRate[ix0_sub0:ix1_sub0,iy0_sub0:iy1_sub0])
        maxStrainRate_sub1[it] = np.max(strainRate[ix0_sub1:ix1_sub1,iy0_sub1:iy1_sub1])
    #    maxStrainRate_sub2[it] = np.max(strainRate[ix0_sub2:ix1_sub2,iy0_sub2:iy1_sub2])
        
        avStrainRate_sub0[it] = np.sum(strainRate[ix0_sub0:ix1_sub0,iy0_sub0:iy1_sub0])/((ix1_sub0-ix0_sub0)*(iy1_sub0-iy0_sub0))
        avStrainRate_sub1[it] = np.sum(strainRate[ix0_sub1:ix1_sub1,iy0_sub1:iy1_sub1])/((ix1_sub1-ix0_sub1)*(iy1_sub1-iy0_sub1))
    #    avStrainRate_sub2[it] = np.sum(strainRate[ix0_sub2:ix1_sub2,iy0_sub2:iy1_sub2])/((ix1_sub2-ix0_sub2)*(iy1_sub2-iy0_sub2))
        
        posMaxStrainRate_sub0[it] = np.argmax(strainRate[ix0_sub0:ix1_sub0,iy0_sub0:iy1_sub0])
    #    maxStrainRate_sub0 = np.max(strainRate)
        
        strainRate_profile_sub0[:,it] = strainRate[ix0_sub0:ix1_sub0,iy0_sub0]
        strainRate_profile_sub1[:,it] = strainRate[ix0_sub0:ix1_sub0,iy0_sub1]
    #    
    ## end of time loop
    #plt.figure(2)
    #plt.clf()
    #plt.subplot(2,1,1)
    #plt.plot(np.log10(maxStrainRate_sub0))
    #plt.plot(np.log10(maxStrainRate_sub1))
    ##plt.plot(np.log10(maxStrainRate_sub2))
    #plt.subplot(2,1,2)
    #plt.plot(np.log10(avStrainRate_sub0))
    #plt.plot(np.log10(avStrainRate_sub1))
    ##plt.plot(np.log10(avStrainRate_sub2))
    #



loadedData = np.load("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/FigData/FaultActivation.npz");
strainRate_profile_sub0     = loadedData["strainRate_profile_sub0"][()]
strainRate_profile_sub1 = loadedData["strainRate_profile_sub1"][()]


plt.figure(3)
plt.clf()

vmin = -14
vmax = -10
plt.set_cmap("wcyrk")

# sample rates (for tests)
sx = 1
st = 1
dt = Setup.Numerics.dtMin/yr

ax1 = plt.subplot(2,1,1)
#plt.pcolor(np.log10(strainRate_profile_sub0),vmin=vmin,vmax=vmax)
#plt.pcolor(np.log10(strainRate_profile_sub0[0::sx,0::st]),vmin=vmin,vmax=vmax)
plt.pcolor(np.arange(0,nt,st)*dt,np.arange(0,nx,sx)*dx,np.log10(strainRate_profile_sub0[0::sx,0::st]),vmin=vmin,vmax=vmax)
#plt.colorbar()
Pos = ax1.get_position()

plt.ylabel("x [m]",fontname="Times New Roman",size=16)


A = np.array([Pos.x0,Pos.y0, Pos.width, Pos.height])+np.array([0.0,-Pos.height-0.01,0.0,0.0])
ax2 = plt.axes(A)
#ax2 = plt.subplot(2,1,2)
#plt.pcolor(np.log10(strainRate_profile_sub1),vmin=vmin,vmax=vmax)
P = plt.pcolor(np.arange(0,nt,st)*dt,np.arange(0,nx,sx)*dx,np.log10(strainRate_profile_sub1[0::sx,0::st]),vmin=vmin,vmax=vmax)


Setup.Numerics.dtMin/yr


plt.ylabel("x [m]",fontname="Times New Roman",size=16)
plt.xlabel("time [yr]",fontname="Times New Roman",size=16)

#plt.colorbar()
#plt.subplot(3,1,3)
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#divider = make_axes_locatable(ax2)
Pos = ax2.get_position()
A = np.array([Pos.x0+Pos.width/4.0,Pos.y0, Pos.width/2.0, Pos.height*.05])+np.array([0.0,+0.035,0.0,0.0])
cax = plt.axes(A)
#text.set_font_properties(font)
#font = matplotlib.font_manager.FontProperties(family='times new roman', style='italic', size=16)
cbar = plt.colorbar(P,cax=cax,orientation="horizontal",ticks=np.arange(vmin,vmax+0.001))

#plt.rc('text', usetex=False)
#plt.rc('font', family='serif')
plt.title("log10 strain rate [1/s]",fontname="Times New Roman",size=16)


#ax1.get_xaxis().set_visible(False)

for l in ax1.yaxis.get_ticklabels():
    l.set_family("Times New Roman")
    l.set_size(12)

    
for l in ax2.yaxis.get_ticklabels():
    l.set_family("Times New Roman")
    l.set_size(12)
    
for l in ax2.xaxis.get_ticklabels():
    l.set_family("Times New Roman")
    l.set_size(12)
    
for l in cbar.ax.xaxis.get_ticklabels():
    l.set_family("Times New Roman")
    l.set_size(12)


ax1.tick_params(direction='in',top=True,right=True)
ax2.tick_params(direction='in',top=True,right=True)


for l in ax1.xaxis.get_ticklabels():
    l.set_visible(False)
    
ax1.xaxis.get_ticklabels()[0].set_text('100')
 

ax1.text(0.01*nt*dt,0.02*nx*dx,"a",fontname="Times New Roman",size=20)
ax2.text(0.01*nt*dt,0.02*nx*dx,"b",fontname="Times New Roman",size=20)
   
plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/FaultActivity.png",r=500)
#
#plt.subplot(3,1,3)
#Vback = abs(Setup.BC.Stokes.backStrainRate*Wbox)
#plt.pcolor(1.0-VxTopo/(Vback),vmin=0.0,vmax=1.0)
#plt.colorbar(orientation='horizontal')










#np.savez("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/FigData/FaultActivation.npz",
#             strainRate_profile_sub0 = strainRate_profile_sub0,
#             strainRate_profile_sub1 = strainRate_profile_sub1,
#             )
#else:


#
#
### Plot the last time step
#iStep = nSteps-1
#rootFolder = superRootFolder + superDirList[iSim] + "/"
#outFolder = "Out_%05d" % iStep #DirList[iStep]
#print("iStep = %i/%i" % (iStep, nSteps-1))
## index of the first node that register the minimum of khi (where khi is active)
## Set file
## =====================
#
#
##Setup = Output.readInput(rootFolder +  'Input/input.json')sa
##Char = Setup.Char
##CharExtra = Input.CharExtra(Char)
#
#dataFolder = rootFolder + outFolder + "/"
#State       = Output.readState(dataFolder + "modelState.json")
#timeSim   = (State.time+ State.dt) * Setup.Char.time
#
#phase = Output.getData(dataFolder + 'phase.bin').data
#mask = phase == 0
#
#strainRate = Output.getData(dataFolder + 'strainRate.bin',True,mask).data
#strain      = Output.getData(dataFolder + 'strain.bin',True,mask).data
#
## Strain Rate    
#cAx_srMin = -14
#cAx_srMax = -10
## Strain
#cAx_sMin = 0.0
#cAx_sMax = 3.0
#
#
## Plot image
#plt.figure(1)
#plt.clf()
#plt.subplot(2,1,1)
##    plt.subplot(nSim,1,iSim+1)
#plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#plt.subplot(2,1,2)
#plt.pcolor(xv - dx/2.0,yv - dy/2.0,strain/strainUnit,vmin=cAx_sMin,vmax=cAx_sMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
###    plt.pcolor(xv - dx/2.0,yv - dy/2.0,diffStrain/strainUnit,vmin=cAx_dsMin,vmax=cAx_dsMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#plt.colorbar()
#plt.axis('equal')
##plt.set_cmap("Pressure")
##plt.title("tstep = %i" % (it))
#
##    refP = np.mean(Pressure)*10.0
##    

    