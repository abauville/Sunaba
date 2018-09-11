#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 14:04:05 2018

@author: abauville
"""
import os
import sys
import numpy as np
from numpy import array as arr
Machine = 2
if Machine==0:
    outFile = "/Users/abauville/Output/Paper_Decollement/Figz/Data/slopes.npz"
    superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/wWater/Beta00/"
    sys.path.insert(0, '../../src/UserInput')
elif Machine==2:
    outFile = "/work/G10501/abauville/Paper_Decollement/PostProcessing/Data/slopes.npz"
    superRootFolder = "/work/G10501/abauville/Paper_Decollement/wWater/Beta00/"
    sys.path.insert(0, '/work/G10501/abauville/Software/StokesFD/src/UserInput')
    
import InputDef as Input
import OutputDef as Output

weakList = arr([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])*100.0
#weakList = arr([1, 5, 10, 20, 40, 60, 80])
#weakList = arr([1, 5])
nW = len(weakList)
LambdaList = arr([00,40,60,80])
nL = len(LambdaList)


IL = np.zeros(nL,dtype=np.int)
IW = np.zeros(nW,dtype=np.int)

#i = 0
#for Lambda in LambdaList/100.0:
#    IL[i] = np.argmin(abs(LambdaRef_list-Lambda))
#    i+=1
#
#i = 0
#for weak in weakList/100.0:
#    IW[i] = np.argmin(abs(chi_list-weak))
#    i+=1
    
    
beta = 0.0
alphas_diff_bar  = np.zeros([nL,nW])
IRefColorMap_Big = np.zeros([nL,nW],dtype=np.int)
IRefColorMap = np.zeros(nL*nW,dtype=np.int)
Lambdas = np.zeros(nL*nW)
chis = np.zeros(nL*nW)
#taper_angles = np.zeros([nL,nW])

#
## find colors in the colormap
#cmin = -1.0
#cmax = 1.0
#
#nC = CMAP_ref.shape[0]
#refCmapValues = np.linspace(cmin,cmax,nC)
#
#iiL = 0
#for iL in IL:
#    iiW = 0
#    for iW in IW:
#        iB = np.argmin(np.abs(betas_all[iL,iW,:]-beta))
##        alphas_diff[iL,iW] = alphas_diff_all[iL,iW,iB]
##        alphas_width[iL,iW] = alphas_WB_up_all[iL,iW,iB] - alphas_WB_low_all[iL,iW,iB]
##        alphas_WB_up[iL,iW] = alphas_WB_up_all[iL,iW,iB]
##        taper_angles[iL,iW] = betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB]
#        alphas_diff_bar[iiL,iiW] = alphas_diff_all[iL,iW,iB]/(betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB])
#
##        if (Lambdas[iL,iW]+chis[iL,iW]<1.0):
##            alphas_diff_bar[iiL,iiW] = (Lambdas_Ref_all[iL,iW,iB]-1.0+chis_all[iL,iW,iB])
##        else:
##            alphas_diff_bar[iiL,iiW] = Lambdas_Ref_all[iL,iW,iB]*chis_all[iL,iW,iB]
#        IRefColorMap_Big[iiL,iiW] = np.argmin(np.abs(refCmapValues-alphas_diff_bar[iiL,iiW]))
#        iiW+=1
#    iiL+=1




## Create the folder tree
# ================================

#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Static2/Beta0/"

#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/NoTopo/Beta0/Weak%i%s/" % (Weak, wWater)
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"

superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    
ProductionMode = False




superDirList = []
i = 0
iL = 0
Lambdas = np.zeros(nW*nL)
chis = np.zeros(nW*nL)
IL = np.zeros(nW*nL,dtype=np.int)
IW = np.zeros(nW*nL,dtype=np.int)
for Lambda in LambdaList:
    iW = 0
    for weak in weakList:
        superDirList.append("Weak%02d/Lambda%02d" % (weak, Lambda))
#        IRefColorMap[i] = IRefColorMap_Big[iL,iW]
        Lambdas[i] = Lambda/100.0
        chis[i] = weak/100.0
        IL[i] = iL
        IW[i] = iW
        i+=1
        iW+=1
    iL+=1




rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]


rootFolders = [''] * len(superDirList) # initialize to the right size
nStepsList = np.zeros(len(superDirList),dtype='int16');
for i in range(len(superDirList)):
    rootFolders[i] = superRootFolder + superDirList[i] + "/Output/"

    
    stepFolderList = os.listdir(rootFolders[i])
    try:
        stepFolderList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char

#pushVel = 10.0*cm/yr



#print("koko")

#Ax = []
#nrows = nW
#ncols = nL
#fig    = Figz_Utils.Figure(2,height=29.7,mode='draft')
##fig    = Figz_Utils.Figure(6,height=13.0,mode='draft')
##AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
##    AxesDum['12'].axis('off')
#Axes   = Figz_Utils.makeAxes(fig,nrows,ncols,aspectRatio=0.8,leftMarginPad=1.5,rightMarginPad=1.0,topMarginPad = 0.0)
#for i in range(nrows):
#    for j in range(ncols):
#        Ax.append(Axes['%i%i' % (i+1,j+1)])

#print("soko")
  
Hgrid = 64 # Thickness H in terms of grid points

Setup = Output.readInput(rootFolders[0] +  'Input/input.json')
ix0_surf = 1
ix1_surf = Setup.Grid.nxC+1
iy0_surf = Hgrid-5
iy1_surf = Hgrid+5


ix0_base = ix0_surf
ix1_base = ix1_surf
iy0_base = 0
iy1_base = 5

ix0_mid = ix0_surf
ix1_mid = ix1_surf
iy0_mid = Hgrid-13
iy1_mid = Hgrid-8


iSub = 0
Colors = ([1.0,0.2,0.4],[0.2,0.4,1.0])
xFront = []
xBase = []
timeList = []

strainFront = []
strainBase = []



nSim = len(superDirList)

#    nSim = 7
iSim0 = 0#nSim-1
#Hsed = 2.0 * km
nSteps = 747#nStepsList[-1]
i0 = 746#nSteps-1#jump-2
#    iStep = i0
jump = 1000000
frame = 0

iCol = 0
iRow = 0
slopes = []#np.zeros(nSim)
timeLists = []
xFronts = []
xBases = []
xMids = []
outFolders = []

alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)

for iSim in range(iSim0,nSim):
    i0 = 0
    jump = 1
    nSteps = nStepsList[iSim]
    slope = np.zeros(len(np.arange(i0,nSteps,jump)))
    timeList = np.zeros(len(np.arange(i0,nSteps,jump)))
    xFront = np.zeros(len(np.arange(i0,nSteps,jump)))
    xBase = np.zeros(len(np.arange(i0,nSteps,jump)))
    xMid = np.zeros(len(np.arange(i0,nSteps,jump)))
    print(rootFolders[iSim])
    
#    list_OutFolder = os.listdir(rootFolders[iSim])
#    list_OutFolder.sort()
    
    
    
    iStep=-1
    for iTimeStep in range(i0,nSteps,jump):
        iStep+=1
        ## Set Folder
        # ================================
        
        outFolder = 'Out_%05d' % iTimeStep
        dataFolder = rootFolders[iSim] + outFolder + "/"
        
        if ((iTimeStep-1)%20==0):
            print("iStep = %i/%i" % (iTimeStep, nSteps-1))
#            print(outFolder)
        
        ## Get Data and Pattern
        # ================================
        Char = Output.readInput(rootFolders[iSim] +  'Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
#            Setup = Output.readInput(rootFolders[iSim] +  'Input/input.json')
#            Char = Setup.Char
        
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        #phase = Output.getData(dataFolder + 'phase.bin').data
        rawData  = Output.getData(dataFolder + 'phase.bin',True)
        phase = rawData.data
        H = 1.0
        ny = rawData.ny
        ymin = rawData.ymin/H
        ymax = rawData.ymax/H
        nx = rawData.nx
        xmin = rawData.xmin/H
        xmax = rawData.xmax/H
        dx = (xmax-xmin)/(nx-1)
        ymin = rawData.ymin/H
        ymax = rawData.ymax/H
        dy = (ymax-ymin)/(ny-1)
        IHsed = 64
        
        strain  = Output.getData(dataFolder + 'strain.bin',True).data
        strain[phase==0] = 0.0
        A = np.where(strain[ix0_surf:ix1_surf,iy0_surf:iy1_surf]>0.5)
        B = np.where(strain[ix0_base:ix1_base,iy0_base:iy1_base]>0.5)
        C = np.where(strain[ix0_mid:ix1_mid,iy0_mid:iy1_mid]>0.5)

        
        iFront = 0
        try:
            iFront           = A[0][0]
            xFront[iStep]    = A[0][0]
        except IndexError:
            xFront[iStep]    = 0

                
        try:
            xBase[iStep]     = B[0][0]
        except IndexError:
            xBase[iStep]    = 0
            
        try:
            xBase[iStep]     = B[0][0]
        except IndexError:
            xBase[iStep]    = 0
            
        try:
            xMid[iStep]     = C[0][0]
        except IndexError:
            xMid[iStep]    = 0
        
        
        Y = np.ones(phase.shape) * dy
        Y = np.cumsum(Y,1) + ymin
        Y[phase==0] = 0.0
        Topo = np.max(Y,1)
        iMaxTopo = np.argmax(Topo)
        if iMaxTopo == 0:
            iMaxTopo = 1
        if iMaxTopo<=iFront:
            iMaxTopo = nx
            
        IBack = nx-int(1.5*IHsed)
        if iFront<IBack-int(IHsed/2.0):
            TopoPrism = Topo[iFront:IBack]
            x = np.arange(IBack-iFront) * dx
        else:
#            TopoPrism = Topo[iFront:iMaxTopo]
#            x = np.arange(iMaxTopo-iFront) * dx
            TopoPrism = Topo[iFront::]
            x = np.arange(nx-iFront) * dx
        
        p = np.polyfit(x,TopoPrism,1)
        slope[iStep] = np.arctan(p[0])
        timeList[iStep] = timeSim
        
        
        
    # end iStep
    
    slopes.append(slope)
    timeLists.append(timeList)
    xFronts.append(xFront)
    xBases.append(xBase)
    xMids.append(xMid)
    outFolders.append(outFolder)
    
    
# end iSim
    



    
np.savez(outFile,
         slopes = slopes,
         timeLists = timeLists,
         xFronts=xFronts,
         xBases=xBases,
         xMids=xMids,
         outFolders=outFolders)