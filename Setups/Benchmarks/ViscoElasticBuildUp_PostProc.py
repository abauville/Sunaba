# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
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
Myr     = 1e6       * yr



rootFolder = "/Users/abauville/StokesFD_Output/ViscoElasticBuildUp/"
DirList = os.listdir(rootFolder)
DirList.remove('.DS_Store')
DirList.remove('Input')

nSteps = len(DirList)
time = np.zeros(nSteps)
sxx = np.zeros(nSteps)
iStep = 0
for outFolder in DirList:

    # Set file
    # =====================
    #outFolder = "Out_00009"
    dataFolder = rootFolder + outFolder + "/"
    
    # Read data
    # =====================
    state = Output.readState(dataFolder + '/modelState.json')
    thisData = Output.getData(dataFolder + 'Sigma_xx.bin')
    
    # Write data
    # =====================
    time[iStep] = state.time
    sxx[iStep] = thisData.data[round(thisData.nx/2),round(thisData.ny/2)]
    
    iStep += 1;




## Define grid
## =====================
#x = np.linspace(thisData.xmin,thisData.xmax,thisData.nx)
#y = np.linspace(thisData.ymin,thisData.ymax,thisData.ny)
#
#xv, yv = np.meshgrid(x,y)
#xv = np.transpose(xv)
#yv = np.transpose(yv)
#
## Plot
## =====================
#data = np.transpose(thisData.data)
#plt.pcolor(xv,yv,np.log10(data))
#plt.colorbar()
#plt.axis('equal')
#plt.show()


# Analytical solution
# =====================
Setup = Output.readInput(rootFolder +  'Input/input.json')

Char = Setup.Char
CharExtra = Input.CharExtra(Char)

time = time * Char.time
sxx = sxx * CharExtra.stress


Exx = Setup.BC.Stokes.backStrainRate
MatProps_Matrix = Setup.MatProps['0']
eta = MatProps_Matrix.getRefVisc(1.0,1.0,Exx)
G   = MatProps_Matrix.G

sxx = -sxx
Exx = -Exx


timeAna = np.linspace(0.0,np.max(time),100)
sxxAna = 2.0*Exx*eta * (1.0 - np.exp(-timeAna * G/eta))


## Plot
## =====================
plt.plot(timeAna/yr,sxxAna,'k')
plt.plot(time/yr,sxx,'ok')
plt.show()
