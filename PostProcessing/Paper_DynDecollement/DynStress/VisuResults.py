# Output reading test

import sys
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
rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+08_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+08_C_1.00e+07_fric_1.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+09_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+10_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
#rootFolder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_5.00e+20_C_1.00e+07_fric_3.00e+01_Hsed_1.00e+03/" % (128*Fac, 64*Fac)
simFolder  = ""
inFolder  = "Input/"
outFolder = "Out_00030/"


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

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
sigmaII = dataSet.data


xmin    = dataSet.xmin
xmax    = dataSet.xmax
ymin    = dataSet.ymin
ymax    = dataSet.ymax
W       = xmax-xmin
H       = ymax-ymin

# Define grid
# =====================
x = np.linspace(xmin,xmax,dataSet.nx)
y = np.linspace(ymin,ymax,dataSet.ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)



# Choose a cell to monitor
# =====================
ixCellMin = 83*Fac
ixCellMax = 83*Fac
ixCell = round((ixCellMin+ixCellMax)/2.0)
iyCell = 32*Fac


# Plotting
# =====================
plt.figure(1)
plt.clf()
plt.pcolor(xv,yv,sigmaII*CharExtra.stress/MPa,vmin=0.0, vmax=100)
plt.plot(xv[ixCell,iyCell], yv[ixCell,iyCell],'og')
plt.axis('equal')
plt.colorbar()

nSteps = 31;
dataEvo = np.zeros(nSteps)
timeEvo = np.zeros(nSteps)
for it in range(0,nSteps) :
    outFolder   = "Out_%05d/" % (it)
    State       = Output.readState(rootFolder + simFolder + outFolder + "modelState.json")
    
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'khi.bin')
    khi = dataSet.data

    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
    sigmaII = dataSet.data * CharExtra.stress
    
    subset_sigmaII  = sigmaII[ixCellMin:ixCellMax+1,iyCell]
    #subset_khi      = khi    [ixCellMin:ixCellMax+1,iyCell]
    #I = np.argmin(subset_khi) # find the minimum khi
    I = np.argmin(subset_sigmaII) # find the minimum khi
    print(I)
    dataEvo[it] = subset_sigmaII[I] # get the stress corresponding to the minimum value
    timeEvo[it] = State.time
    
plt.figure(2)
plt.clf()
plt.plot(timeEvo,dataEvo/MPa)
