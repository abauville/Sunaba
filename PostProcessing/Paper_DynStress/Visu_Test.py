# Output reading test

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
#from pprint import pprint

#import scipy
import json

import OutputDef as Output

import InputDef as Input

from matplotlib.colors import LinearSegmentedColormap

from math import pi, tan



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

Pa      = kg/m/s/s
MPa     = 1e6 * Pa
GPa     = 1e9 * Pa

degree = pi/180


# Colormap
# =====================
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.75 , 1.0, 1.0),
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


CMAP = LinearSegmentedColormap('StokesFD', cdict1)
plt.register_cmap(cmap=CMAP)
plt.set_cmap('StokesFD')
#plt.set_cmap("gray")



# Set file
# =====================
RFac = 1
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_Test/"
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
simFolder  = ""
inFolder  = "Input/"
nSteps = Output.getNumberOfOutFolders(rootFolder);
#nSteps = 150
outFolder = "Out_%05d/" % (nSteps-1)


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

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'khi.bin')
khi = dataSet.data

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
P = dataSet.data

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
sigmaII = dataSet.data

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

# Define grid
# =====================
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)

Hmatrix = 1000.0;



# Choose a cell to monitor
# =====================
halfSpan = 76
ixCellMin = 0#76 - halfSpan
ixCellMax = nx#ixCellMin + 2*halfSpan
iyCell = np.argmin( np.abs(y-Hmatrix/2.0) )
ixCell = int(round((ixCellMin+ixCellMax)/2.0))

#ixCellMin = 30
#ixCellMax = ixCellMin#+60
#iyCell = np.argmin( np.abs(y-Hmatrix/4.5) )
#ixCell = round((ixCellMin+ixCellMax)/2.0)




# Plotting 2D
# =====================
plt.figure(1)
plt.clf()
plt.subplot(2,1,1)
#plt.pcolor(xv,yv,sigmaII*CharExtra.stress/MPa,vmin=0.0, vmax=4.0*Setup.Physics.Pback/MPa)
plt.pcolor(xv - dx/2.0,yv - dy/2.0,P*CharExtra.stress/MPa,vmin=1.0*Setup.Physics.Pback/MPa, vmax=3.0*Setup.Physics.Pback/MPa) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
plt.plot(xv[ixCell,iyCell], yv[ixCell,iyCell],'og')
plt.axis('equal')
plt.colorbar(orientation='horizontal')


# Extracting data
# =====================
sigmaIIEvo = np.zeros(nSteps)
#sigmaII_altEvo = np.zeros(nSteps)
PEvo = np.zeros(nSteps)
khiEvo = np.zeros(nSteps)
ZEvo = np.zeros(nSteps)
timeEvo = np.zeros(nSteps)
dtEvo = np.zeros(nSteps)
sigmaYieldEvo = np.zeros(nSteps)
IEvo = np.zeros(nSteps)
for it in range(0,nSteps) :
    outFolder   = "Out_%05d/" % (it)
    State       = Output.readState(rootFolder + simFolder + outFolder + "modelState.json")
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'khi.bin')
    khi = dataSet.data
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'Z.bin')
    Z = dataSet.data

    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
    sigmaII = dataSet.data * CharExtra.stress
    subset_sigmaII  = sigmaII[ixCellMin:ixCellMax+1,iyCell]
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
    P = dataSet.data * CharExtra.stress
    subset_P  = P[ixCellMin:ixCellMax+1,iyCell]
    
    I = np.argmin(subset_P) # find the minimum khi
    sigmaIIEvo[it] = subset_sigmaII[I] # get the stress corresponding to the minimum value
    
    PEvo[it] = subset_P[I] # get the stress corresponding to the minimum value
    sigmaYieldEvo[it] = C*np.cos(phi) + PEvo[it]*np.sin(phi)
    
    khiEvo[it] = khi[ixCell, iyCell]
    ZEvo[it]   = Z  [ixCell, iyCell]
    
    timeEvo[it] = State.time * Setup.Char.time
    dtEvo[it] = State.dt #* Char.time
    
    IEvo[it] = I
    
# Analytical yield stress
# =====================  
P = Setup.Physics.Pback
sigmaYield_back = C*np.cos(phi) + P*np.sin(phi)
sigmaYield_back /= MPa


# Find interesting values
# =====================
I_sigmaMax = np.argmax(sigmaIIEvo)
I_sigmaMin = I_sigmaMax + np.argmin(sigmaIIEvo[I_sigmaMax:])

#Delta_timeSoft = timeEvo[I_sigmaMin] - timeEvo[I_sigmaMax]
    
# Rate of change of sigma
# =====================
timeEvo_centered = (timeEvo[1:]+timeEvo[:-1])/2.0
timeEvo_diff = np.diff(timeEvo)
sigmaIIRateEvo = np.diff(sigmaIIEvo)/timeEvo_diff
I_sigmaIIRateMax = np.argmax(sigmaIIRateEvo)
I_sigmaIIRateMin = np.argmin(sigmaIIRateEvo)
sigmaIIRate_LimitEndOfSoftening = 0.1 * sigmaIIRateEvo[I_sigmaIIRateMin] # arbitrary limit to determine the end of the softening period and the beginning of the stable one
I_EndOfSoftening = I_sigmaIIRateMin + np.argmax(sigmaIIRateEvo[I_sigmaIIRateMin:]>sigmaIIRate_LimitEndOfSoftening)

# Plotting Graph
# =====================
#plt.figure(1)
plt.subplot(2,1,2)
#plt.clf()
timePlotUnit = 1000*yr
sigmaPlotUnit = MPa

plt.plot(timeEvo/timePlotUnit,sigmaYieldEvo/sigmaPlotUnit,'-r')
plt.plot(timeEvo/timePlotUnit,sigmaIIEvo/sigmaPlotUnit,'.k')
#plt.plot(timeEvo/timePlotUnit,sigmaII_altEvo/MPa,'.r')
plt.plot(timeEvo/timePlotUnit,PEvo/sigmaPlotUnit,'.b')

plt.plot((0,timeEvo[-1]/timePlotUnit), (sigmaYield_back,sigmaYield_back),'--k')

plt.plot(timeEvo[I_sigmaMax]/timePlotUnit,sigmaIIEvo[I_sigmaMax]/sigmaPlotUnit,'or',markerSize=8,markerFaceColor='none')
#plt.plot(timeEvo[I_sigmaMin]/timePlotUnit,sigmaIIEvo[I_sigmaMin]/MPa,'ob',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo[I_EndOfSoftening]/timePlotUnit,sigmaIIEvo[I_EndOfSoftening]/sigmaPlotUnit,'ob',markerSize=8,markerFaceColor='none')

plt.title("$P_{back}$ = %.f MPa, C = %.f MPa, G = %.f GPa" % (Setup.Physics.Pback/MPa, Setup.MatProps['1'].cohesion/MPa, Setup.MatProps['1'].G/GPa))
plt.legend(["$\\tau_{y}$","$\\tau_{II}$","P","$\\tau_{y}$ at $P_{back}$"])
plt.xlabel("time [kyr]")
plt.ylabel("Stress [MPa]")

plt.figure(2)
plt.clf()

plt.plot(timeEvo_centered/timePlotUnit, sigmaIIRateEvo/(MPa/timePlotUnit),'-ok')
plt.plot(timeEvo_centered[I_sigmaIIRateMax]/timePlotUnit,sigmaIIRateEvo[I_sigmaIIRateMax]/(MPa/timePlotUnit),'or',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo_centered[I_sigmaIIRateMin]/timePlotUnit,sigmaIIRateEvo[I_sigmaIIRateMin]/(MPa/timePlotUnit),'ob',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo_centered[I_EndOfSoftening]/timePlotUnit,sigmaIIRateEvo[I_EndOfSoftening]/(MPa/timePlotUnit),'og',markerSize=8,markerFaceColor='none')

Delta_timeSoft = timeEvo_centered[I_EndOfSoftening] - timeEvo[I_sigmaMax]

print("Delta_timeSoft = %.f Kyr" % (Delta_timeSoft/Kyr))
print("sigmaMax = %.f MPa, sigmaMin = %.f MPa" % (sigmaIIEvo[I_sigmaMax]/MPa, sigmaIIEvo[I_sigmaMin]/MPa))
print("sigmaMaxTime = %.f Kyr, sigmaMinTime = %.f Kyr" % (timeEvo[I_sigmaMax]/Kyr, timeEvo[I_sigmaMin]/Kyr))



plt.figure(3)
plt.clf()
plt.subplot(2,1,1)
plt.plot(timeEvo/timePlotUnit, np.log10(ZEvo) ,'or')
plt.subplot(2,1,2)
plt.plot(timeEvo/timePlotUnit, np.log10(khiEvo) ,'or')

plt.figure(4)
plt.clf()
plt.subplot(2,1,1)
plt.plot(timeEvo/timePlotUnit,dtEvo,'ok')
plt.subplot(2,1,2)
plt.plot(timeEvo/timePlotUnit,IEvo,'ok')