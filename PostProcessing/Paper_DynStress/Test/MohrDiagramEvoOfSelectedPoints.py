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

Setup = Output.readInput(rootFolder +  'Input/input.json')
nSteps = len(DirList)
time = np.zeros(nSteps)
sxx = np.zeros(nSteps)
P = np.zeros(nSteps)
eII = np.zeros(nSteps)
iStep = 0
Char = Setup.Char
CharExtra = Input.CharExtra(Char)
Ix_t = np.zeros(nSteps)
Iy_t = np.zeros(nSteps)


for outFolder in DirList:

    # Set file
    # =====================
    #outFolder = "Out_00009"
    dataFolder = rootFolder + outFolder + "/"
    
    # Read data
    # =====================
    state = Output.readState(dataFolder + 'modelState.json')
    
#    thisData = Output.getData(dataFolder + 'sigma_II.bin')
#    thisData = Output.getData(dataFolder + 'P.bin')
    thisData = Output.getData(dataFolder + 'khi.bin')
#    iyCell = int(thisData.ny/2)
#    iyCell = 103#thisData.ny-2
#    iyCell = 65#thisData.ny-2
#    iyCell = 20#thisData.ny-2
    iyCell = 85#thisData.ny-2
    Ix = int(np.argmin(thisData.data[:,iyCell]) )
    
#    ixCell = 89#thisData.nx-2
#    Iy = np.argmin(thisData.data[ixCell,0:103]) 
    
#    I = np.argmax(thisData.data[:,:]) 
    Ix_t[iStep] = Ix
#    Iy_t[iStep] = Iy
    iStep += 1;



#I = int(max(I_t)-4)


iStep = 0
for outFolder in DirList:
#    Ix = Ix_t[iStep]
#    Iy = Iy_t[iStep]
#    Ix = 1
#    Ix = 91
#    Ix = 27
    Ix = Ix_t[int(np.argmax(Ix_t>0))] # index of the first node that register the minimum of khi (where khi is active)
    # Set file
    # =====================
    #outFolder = "Out_00009"
    dataFolder = rootFolder + outFolder + "/"
    
    # Read data
    # =====================
    state = Output.readState(dataFolder + 'modelState.json')
    # Write data
    # =====================
    time[iStep] = state.time
    thisData = Output.getData(dataFolder + 'P.bin')
    P[iStep] = thisData.data[Ix,iyCell]
#    P[iStep] = thisData.data[ixCell,Iy]
#    P[iStep] = np.max(thisData.data)
    thisData = Output.getData(dataFolder + 'Sigma_II.bin')
    sxx[iStep] = thisData.data[Ix,iyCell]
#    sxx[iStep] = thisData.data[ixCell,Iy]
#    sxx[iStep] = np.max(thisData.data)
    thisData = Output.getData(dataFolder + 'strainRate.bin')
    eII[iStep] = thisData.data[Ix,iyCell]
#    eII[iStep] = thisData.data[ixCell,Iy]
#    eII[iStep] = np.max(thisData.data)
    
    
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


time = time * Char.time
sxx = sxx * CharExtra.stress
P = P * CharExtra.stress
eII = eII*1.0/Char.time

Exx = np.abs(Setup.BC.Stokes.backStrainRate)
MatProps_Matrix = Setup.MatProps['1']
eta = MatProps_Matrix.getRefVisc(1.0,1.0,Exx)
G   = MatProps_Matrix.G

#sxx = -sxx
#Exx = -Exx

timeAna = np.linspace(0.0,np.max(time),100)
sxxAna = 2.0*Exx*eta * (1.0 - np.exp(-timeAna * G/eta))

Pback= 100.0*MPa
C = 50.0*MPa
phi = 30/180*np.pi
Sy_back = C*np.cos(phi) + Pback*np.sin(phi)

Sy_DP = C*np.cos(phi) + P*np.sin(phi) # SY drucker prager

dt = Setup.Numerics.dtMin
nSteps = int(np.floor(np.max(time)/dt)) + 1


## Plot
## =====================
sigmaPlotUnit = Pback
plt.figure(1)
plt.clf()

plt.plot(timeAna/Kyr,sxxAna/sigmaPlotUnit,'k')
plt.plot(time/Kyr,sxx/sigmaPlotUnit,'.r')
plt.plot(time/Kyr,P/sigmaPlotUnit,'.b')
plt.plot(time/Kyr,Sy_DP/sigmaPlotUnit,'-g')
plt.plot(timeAna[[0,-1]]/Kyr,np.array([Sy_back,Sy_back])/sigmaPlotUnit,'k')

plt.show()





## In terms of principal stresses:
## =====================
#TauII = C*cos(phi) + P*sin(phi)
#with, TauII = (S1-S3)/2
#        P   = (S1+S3)/2


# Substistuting and rearranging yields:
# S1 = 1/(1-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )

# Using S3 = Pback
S3 = Pback
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
TauII_Lim = (S1-S3)/2.0
P_Lim = (S1+S3)/2.0

#S3 = Pback/2.0
#S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
#TauII_LimHalfS3 = (S1-S3)/2.0


plt.plot(timeAna[[0,-1]]/Kyr,np.array([TauII_Lim,TauII_Lim])/sigmaPlotUnit,'--r')
#plt.plot(timeAna[[0,-1]]/Kyr,np.array([TauII_LimHalfS3,TauII_LimHalfS3])/sigmaPlotUnit,'-r')
plt.plot(timeAna[[0,-1]]/Kyr,np.array([P_Lim,P_Lim])/sigmaPlotUnit,'--b')
plt.plot(timeAna[[0,-1]]/Kyr,np.array([S1,S1])/sigmaPlotUnit,'--g')


S1_num = sxx+P
S3_num = -sxx+P
plt.plot(time/Kyr,S1_num/sigmaPlotUnit,'.g')
plt.plot(time/Kyr,S3_num/sigmaPlotUnit,'.y')


plt.figure(2)
plt.clf()
plt.plot(time/Kyr,np.log10(eII),'.k')
