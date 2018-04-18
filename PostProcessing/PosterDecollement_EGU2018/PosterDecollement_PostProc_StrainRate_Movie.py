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


#superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNew/FixedDt_Method1/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNew/FixedDt_Method0/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/Test00/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/Corotational/FixedDt_Method0/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/StrucStyle/NoOutFlow_phib25_cohW_50_fricW75_cStrain1/Output/"
#superRootFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/StrucStyle/NoOutFlow_phib29_cohW_50_fricW75_cStrain1/Output/"

superRootFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/StrucStyle_WeakOrNotBase/NoOutFlow_phib30_cohW_50_fricW85_cStrain1_HFac1_G5e8_30it_NotWeakBase/Output/"

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



#time_t = np.zeros(len(range(i0,nSteps,jump)))
#Pfault_t    = np.zeros(len(range(i0,nSteps,jump)))
#Pfar_t      = np.zeros(len(range(i0,nSteps,jump)))
#maxStrain_Sim = np.zeros(nSim)


#
#if sys.platform == "linux" or sys.platform == "linux2":
#    # linux
#    os.system("mkdir /home/abauvill/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
#elif sys.platform == "darwin":
#    # OS X
#    os.system("mkdir /Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Movies/" + superDirList[iSim])
    

iStep = 10
nSim = len(superDirList)


i0 = 0
#nSteps = 3701
nSubPlots = 1
jump=1# round(nSteps/nSubPlots)

nt = len(range(i0,nSteps,jump))
maxStrainRate_sub0 = np.zeros(nt)
maxStrainRate_sub1 = np.zeros(nt)
maxStrainRate_sub2 = np.zeros(nt)
avStrainRate_sub0 = np.zeros(nt)
avStrainRate_sub1 = np.zeros(nt)
avStrainRate_sub2 = np.zeros(nt)
posMaxStrainRate_sub0 = np.zeros(nt)

xminP = xmin#*0.75
xmaxP = xmax
yminP = ymin
ymaxP = ymax#/2.0

# ColorScales
# Strain Rate    
cAx_srMin = -13
cAx_srMax = -11
# Strain
cAx_sMin = 0.0
cAx_sMax = 10.0

axDict = dict()
nSim = 1
nSubPlots = 1
iSub = 0
plt.figure(2)
plt.clf()

for iVert in range(0,nSubPlots):
    for iHor in range(0,nSim):
        iSub += 1
        xPad = 0.05
        yPad = 0.05
        
        yPadTop = 0.05
        yPadBot = 0.15
        
        xPadLeft = 0.05
        xPadRight = 0.05
        
        WFig = 1.0 - xPadLeft - xPadRight
        WSub = (WFig-(nSim-1)*xPad)/nSim
        HFig = 1.0 - yPadTop - yPadBot
        HSub = (HFig-(nSubPlots-1)*yPad)/nSubPlots
#            axDict["ax%i" % iSub] = plt.axes([xPad*(iHor+1)+iHor*WFig  ,yPad*(iVert+1)+iVert*HSub,WFig/nSim,HFig/nSubPlots])
        axDict["ax%i" % iSub] = plt.axes([xPadLeft+xPad*(iHor)+iHor*WSub,1.0-yPadTop-yPad*iVert-(iVert+1)*HSub, WSub, HSub])
        



# Custom colorbar
WcBar = 0.5
HcBar = 0.025
axDict["axColorBar"] = plt.axes([0.5-WcBar/2.0, 0.15, WcBar, HcBar],xticks=np.arange(cAx_srMin,cAx_srMax+0.0000001,0.5))

img = np.ones([2,128,4])
dataPlot = np.ones([2,128])
dataPlot[0,:] = np.linspace(0,1,128)
dataPlot[1,:] = dataPlot[0,:]
cmap = plt.get_cmap("wcyrk")
#dataPlotNorm = Normalize(cAx_srMin, cAx_srMax, clip=True)(dataPlot)
img = cmap(dataPlot)
# alpha:
dataPlot = Normalize(0.0, 0.5, clip=True)(dataPlot)
img[:,:,-1] = dataPlot

imgplot = plt.imshow(img,extent=[cAx_srMin,cAx_srMax,0,1],aspect="auto",cmap=cmap)
plt.title("log10 second strainrate invariant [1/s]")





plt.sca(axDict["ax%i" % 1])


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
#            if iSim == 1:
#            rootFolder = rootFolder.replace("/NoOutFlow_","/OutFlow_")
#            rootFolder = rootFolder.replace("_phib30","_phib20")
            print(rootFolder)
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
#            strainRate = np.ma.masked_array(strainRate, strainRate<1e-12)
#            strain  = Output.getData(dataFolder + 'strain.bin',True,mask).data
#    #        strain = np.ma.masked_array(strain, strain<3.0)
#
            
            plt.cla()
            plt.axis("off")
            
#            plt.subplot(nSubPlots,2,iSub)

            dataFolderImage = rootFolder.replace("/Output/dx","/Visu/dx")
            img = mpimg.imread(dataFolderImage + 'Frame_%.5d.png' % (iStep))
    #        imgplot = ax1.imshow(img,extent=[xmin-Wbox*0.05,xmax+Wbox*0.05,ymin-Hbox*0.05,ymax+Hbox*0.05])
            scalePad = 0.025
            padY = -Hbox*0.125
            imgplot = plt.imshow(img,extent=[xmin-Wbox*scalePad,xmax+Wbox*scalePad,ymin-Hbox*scalePad+padY,ymax+Hbox*scalePad-padY])
            
#            plt.pcolormesh(xv,yv+dy,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax,shading="gouraud") # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center            
            plt.axis([xminP,xmaxP,yminP,ymaxP])
#            
#            if (iSim==1):
#                plt.xlim(xmax, xmin)
            
            
            #### Transparent plot
            img = np.ones([ny,nx,4])
            dataPlot = np.flipud(np.log10(strainRate.T))
#            imgStrain[:,:,0] = strainPlot
#            imgStrain[:,:,1] = strainPlot
#            imgStrain[:,:,2] = strainPlot
#    #        imgStrain[:,:,3] = strainPlot
#            cmap = plt.get_cmap("Pressure")
            cmap = plt.get_cmap("wcyrk")

    #        cmap = plt.cm.RdYlBu

            dataPlotNorm = Normalize(cAx_srMin, cAx_srMax, clip=True)(dataPlot)
            img = cmap(dataPlotNorm)
            # alpha:
            dataPlotNorm = Normalize(cAx_srMin, (cAx_srMin+cAx_srMax)/2.0, clip=True)(dataPlot)
            dataPlotNorm[np.flipud(phase.T)==0] = 0.0

            img[:,:,-1] = dataPlotNorm
    
            imgplot = plt.imshow(img,extent=[xmin,xmax,ymin,ymax],cmap=cmap)

#            imgplot2 = plt.imshow(imgStrainRate,extent=[xmin,xmax,ymin,ymax],cmap=cmap)
            
             #### Transparent plot
             
             
            # Scale
#            if (iSub == 1):
            plt.fill(xminP+np.array([0.0,1000.0+200.0,1000.0+200.0,0.0]),yminP+np.array([0.0,0.0,75.0+350.0,75.0+350.0]),"w")
            plt.fill(xminP+100.0+np.array([0.0,1000.0,1000.0,0.0]),yminP+100.0+np.array([0.0,0.0,75.0,75.0]),"k")
            plt.text(xminP+100.0+500.0 , yminP+75.0+100.0+50.0,"1 km",horizontalAlignment="center")
#            plt.pause(0.0001)
             
            
    #    #    plt.subplot(nSim,1,iSim+1)
    #        plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(strainRate/strainRateUnit),vmin=cAx_srMin,vmax=cAx_srMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    ##        plt.pcolor(xv - dx/2.0,yv - dy/2.0,strain/strainUnit,vmin=cAx_sMin,vmax=cAx_sMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    #    ##    plt.pcolor(xv - dx/2.0,yv - dy/2.0,diffStrain/strainUnit,vmin=cAx_dsMin,vmax=cAx_dsMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#            plt.colorbar()
    #        plt.axis('equal')
    #        plt.set_cmap("Pressure")
    #        plt.title("tstep = %i" % (it))
            
        #    refP = np.mean(Pressure)*10.0
            

#        plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/phiB20/StrainRateEvo/Frame%05d.png" % (iSub-1),dpi=400)
        plt.savefig("/Users/abauville/Dropbox/00_ConferencesAndSeminars/EGU2018/Figz/phiB30/StrainRateEvo/Frame%05d.png" % (iSub-1),dpi=400)