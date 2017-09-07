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
Fac = 1
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+08_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+08_C_1.00e+07_fric_1.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+09_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+10_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+20_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/nx_183_ny_128_G_5.00e+10_C_4.00e+07_fric_3.00e+01_Pref_5.00e+07/"
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/nx_183_ny_128_G_5.00e+10_C_2.00e+06_fric_3.00e+01_Pref_5.00e+07/"
rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/Test/"
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
ixCellMin = 75-5
ixCellMax = ixCellMin+10
iyCell = np.argmin( np.abs(y-Hmatrix/2.0) )
ixCell = round((ixCellMin+ixCellMax)/2.0)

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
timeEvo = np.zeros(nSteps)
sigmaYieldEvo = np.zeros(nSteps)
for it in range(0,nSteps) :
    outFolder   = "Out_%05d/" % (it)
    State       = Output.readState(rootFolder + simFolder + outFolder + "modelState.json")
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'khi.bin')
    khi = dataSet.data

    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
    sigmaII = dataSet.data * CharExtra.stress
    subset_sigmaII  = sigmaII[ixCellMin:ixCellMax+1,iyCell]
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
    P = dataSet.data * CharExtra.stress
    subset_P  = P[ixCellMin:ixCellMax+1,iyCell]
    
#    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_xx.bin')
#    sxx = dataSet.data * CharExtra.stress
#    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_xy.bin')
#    sxy = dataSet.data * CharExtra.stress
#    sigmaII_alt = np.sqrt(sxx**2 + sxy**2)
#    subset_sigmaII_alt  = sigmaII_alt[ixCellMin:ixCellMax+1,iyCell]
    
    #subset_khi      = khi    [ixCellMin:ixCellMax+1,iyCell]
    #I = np.argmin(subset_khi) # find the minimum khi
    I = np.argmin(subset_sigmaII) # find the minimum khi
    sigmaIIEvo[it] = subset_sigmaII[I] # get the stress corresponding to the minimum value
    
    #I = np.argmin(subset_P) # find the minimum khi
    PEvo[it] = subset_P[I] # get the stress corresponding to the minimum value
    sigmaYieldEvo[it] = C*np.cos(phi) + PEvo[it]*np.sin(phi)
    
#    I = np.argmin(subset_sigmaII_alt) # find the minimum khi
#    sigmaII_altEvo[it] = subset_sigmaII_alt[I] # get the stress corresponding to the minimum value
    
    timeEvo[it] = State.time * Setup.Char.time
    
    
# Analytical yield stress
# =====================  
P = Setup.Physics.Pback
sigmaYield_back = C*np.cos(phi) + P*np.sin(phi)
sigmaYield_back /= MPa


# Find interesting values
# =====================
I_sigmaMax = np.argmax(sigmaIIEvo)
I_sigmaMin = I_sigmaMax + np.argmin(sigmaIIEvo[I_sigmaMax:])

Delta_timeSoft = timeEvo[I_sigmaMin] - timeEvo[I_sigmaMax]
    
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
plt.plot(timeEvo/1000/yr,sigmaYieldEvo/MPa,'-r')
plt.plot(timeEvo/1000/yr,sigmaIIEvo/MPa,'.k')
#plt.plot(timeEvo/1000/yr,sigmaII_altEvo/MPa,'.r')
plt.plot(timeEvo/1000/yr,PEvo/MPa,'.b')

plt.plot((0,timeEvo[-1]/1000/yr), (sigmaYield_back,sigmaYield_back),'--k')

plt.plot(timeEvo[I_sigmaMax]/1000/yr,sigmaIIEvo[I_sigmaMax]/MPa,'or',markerSize=8,markerFaceColor='none')
#plt.plot(timeEvo[I_sigmaMin]/1000/yr,sigmaIIEvo[I_sigmaMin]/MPa,'ob',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo[I_EndOfSoftening]/1000/yr,sigmaIIEvo[I_EndOfSoftening]/MPa,'ob',markerSize=8,markerFaceColor='none')

plt.title("$P_{back}$ = %.f MPa, C = %.f MPa, G = %.f GPa" % (Setup.Physics.Pback/MPa, Setup.MatProps['1'].cohesion/MPa, Setup.MatProps['1'].G/GPa))
plt.legend(["$\\tau_{y}$","$\\tau_{II}$","P","$\\tau_{y}$ at $P_{back}$"])
plt.xlabel("time [kyr]")
plt.ylabel("Stress [MPa]")

plt.figure(2)
plt.clf()

plt.plot(timeEvo_centered/1000/yr, sigmaIIRateEvo/(MPa/1000/yr),'-ok')
plt.plot(timeEvo_centered[I_sigmaIIRateMax]/1000/yr,sigmaIIRateEvo[I_sigmaIIRateMax]/(MPa/1000/yr),'or',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo_centered[I_sigmaIIRateMin]/1000/yr,sigmaIIRateEvo[I_sigmaIIRateMin]/(MPa/1000/yr),'ob',markerSize=8,markerFaceColor='none')
plt.plot(timeEvo_centered[I_EndOfSoftening]/1000/yr,sigmaIIRateEvo[I_EndOfSoftening]/(MPa/1000/yr),'og',markerSize=8,markerFaceColor='none')

Delta_timeSoft = timeEvo_centered[I_EndOfSoftening] - timeEvo[I_sigmaMax]