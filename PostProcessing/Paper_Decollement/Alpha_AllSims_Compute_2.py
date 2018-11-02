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
    sys.path.insert(0, '../../src/UserInput')
    sys.path.insert(0, './CriticalTaper')
elif Machine==2:
    sys.path.insert(0, '/work/G10501/abauville/Software/StokesFD/src/UserInput')
import InputDef as Input
import OutputDef as Output
#import matplotlib.pyplot as plt
#from CritTaper_utils import Taper
#from PaperDecollement_Utils import getColormap, get_XYandPattern

#weakList = arr([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])*100.0
#weakList = arr([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80])

#weakList = arr([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80])
#weakList = arr([1, 10, 20, 30, 40, 50, 60, 70, 80])
weakList = arr([1, 30, 60])

#weakList = arr([35])
winSize = 32

sampleRate = 50
pointSize = sampleRate/1.0

jump = 1


#weakList = arr([1, 5])
nW = len(weakList)
#LambdaList = arr([00,40,60,80])
LambdaList = arr([60])
nL = len(LambdaList)

beta = 0.0

if nL==1:
    add = "_Beta%02d_Lambda%02d_WeakDense" % (beta*10.0,LambdaList[0])
else:
    add = "_Beta%02d_all" % (beta*10.0)

if Machine==0:
    outFile = "/Users/abauville/Output/Paper_Decollement/Figz/Data/locSlopes%s.npy" % add
#    superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/wWater/Beta%02d/" % (beta*10.0)
    superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta%02d/" % (beta*10.0)
    
elif Machine==2:
    outFile = "/work/G10501/abauville/Paper_Decollement/PostProcessing/Scripts/Data/locSlopes%s.npy" % add
    superRootFolder = "/work/G10501/abauville/Paper_Decollement/wWater/Beta%02d/" % (beta*10.0)
    
    



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
altSuperDirList = []
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
        altSuperDirList.append("Weak%02d/Done_Lambda%02d" % (weak, Lambda))
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

    try:
        stepFolderList = os.listdir(rootFolders[i])
    except FileNotFoundError:
        rootFolders[i] = superRootFolder + superDirList[i][0:7] +  "Done_" + superDirList[i][7:] + "/Output/"
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
#    iStep = i0
#jump = 1000000
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


Data = {}


for iSim in range(iSim0,nSim):
    
    
    iStart = 0
    nSteps = nStepsList[iSim]
    slope = np.zeros(len(np.arange(iStart,nSteps,jump)))
    timeList = np.zeros(len(np.arange(iStart,nSteps,jump)))
    xFront = np.zeros(len(np.arange(iStart,nSteps,jump)))
    xBase = np.zeros(len(np.arange(iStart,nSteps,jump)))
    xMid = np.zeros(len(np.arange(iStart,nSteps,jump)))
    print(rootFolders[iSim])
    
    Lambda = Lambdas[iSim]
    chi = chis[iSim]
    
    Data['Lambda%02d_chi%02d' % (Lambda*100, chi*100)] = {}
    thisData = Data['Lambda%02d_chi%02d' % (Lambda*100, chi*100)]
#    list_OutFolder = os.listdir(rootFolders[iSim])
#    list_OutFolder.sort()
    
    thisData['tSteps'] = np.arange(iStart,nSteps,jump)
#    thisData['time_list'] = []
    iStep=-1
    
    locSlopes = []
    lenPrisms = []
    
    for iTimeStep in range(iStart,nSteps,jump):
        iStep+=1
        ## Set Folder
        # ================================
        
        outFolder = 'Out_%05d' % iTimeStep
#        try:
#            outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
#        except FileNotFoundError:
#            outFolder = os.listdir(superRootFolder + altSuperDirList[iSim] + "/Output/")[-1]
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
        xGlob = np.linspace(xmin,xmax,nx)
        
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
        
        
#        plt.plot(np.diff(TopoPrism),'.')
        

        locSlope = np.zeros(Topo.shape)  
        for i in range(winSize,len(Topo)-winSize):
            locSlope[i] = np.arctan(np.polyfit(xGlob[i-winSize:i+winSize],Topo[i-winSize:i+winSize],1)[0])
            
        
        lenPrism = len(TopoPrism)-2*winSize
        if lenPrism>0:
            locSlope = np.zeros(lenPrism)  
            i0 = 0
            for i in range(winSize,len(TopoPrism)-winSize):
                locSlope[i0] = np.arctan(np.polyfit(x[i-winSize:i+winSize],TopoPrism[i-winSize:i+winSize],1)[0])
                i0+=1
        else:
            locSlope = np.array([])
                    
        
        
        
        locSlopes.append(locSlope)
        lenPrisms.append(lenPrism)
#        diffTopo = np.diff(Topo)
#        locSlope = np.zeros(diffTopo.shape)
#        for i in range(winSize,len(diffTopo)-winSize):
#            locSlope[i] = 1.0/(2.0*winSize*dx)*np.sum(diffTopo[i-winSize:i+winSize])

#            locSlope[i] = 1.0/(2.0*winSize)*np.sum(Topo[i-winSize:i*winSize])
            
      
        
        
        
#        
### Plotting test
## =============================================================================
#        
#        # =============================================================================
#        #                       Create taper and get data
#        
#        rho_w = 1000.0
#        rho = 2500.0
#        phiRef   = 30.0*np.pi/180.0
#        
##        chi = 1e-7
#        
#        LambdaRef = LambdaList[0]/100.0
#        #chi_list = [.05,0.7,0.7]
##        beta_list = np.linspace(35.0,-5.0,9)/deg
#        beta = 0.0
#        chi = weakList[iSim]/100.0
#        LambdaWeak = (1.0-chi) * LambdaRef   + chi
#        
#        ## ============= RefTaper =================    
#        tprRef = Taper(phi=phiRef, phi_b=phiRef,
#                    Lambda=LambdaRef, Lambda_b=LambdaRef+1e-6,
#                    rho_w=rho_w, rho=rho)
#        tprBasalWeak = Taper(phi=phiRef, phi_b=phiRef,
#                    Lambda=LambdaRef, Lambda_b=LambdaWeak,
#                    rho_w=rho_w, rho=rho)
#        tprTotalWeak = Taper(phi=phiRef, phi_b=phiRef,
#                    Lambda=LambdaRef, Lambda_b=LambdaWeak,
#                    rho_w=rho_w, rho=rho)
#        tprRef.computeAlphaVsBeta(n=2010)
#        tprBasalWeak.computeAlphaVsBeta(n=2010)
#        tprTotalWeak.computeAlphaVsBeta(n=2010)
#        
#        alpha_Ref = tprRef.findAlpha(beta,"average")
#        alpha_WB_up = tprBasalWeak.findAlpha(beta,"lower")
#        alpha_WB_low = tprBasalWeak.findAlpha(beta,"upper")
#        alpha_WF = tprTotalWeak.findAlpha(beta,"average")
#        
#        
#        
#        
##        tpr_list.append(tpr)
#        
#                
#        #                       Create taper and get data
#        # =============================================================================
#        
#        
#        # Plot
#        # ======================
#        
#        plt.figure(205)
#        plt.clf()
#        plt.subplot(311)
#        x0 = 0
#        x1 = nx
#        y0 = 0
#        y1 = 25
#        plt.xlim(x0,x1)
#        plt.ylim(y0,y1)
#        
#        deg = 180.0/np.pi
#        transparency = 0.2
#        Color  = arr([[.85,.15,.25],[.25,.5,.5],[.85,.15,.25]])
#        lineWidth = 0.75
#        plt.plot([x0,x1],[alpha_WF *deg, alpha_WF*deg],'b',linewidth=lineWidth)
#        plt.fill([x0,x1,x1,x0],arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,color=Color[1],alpha=transparency)
#        plt.plot([x0,x1],[alpha_WB_up*deg, alpha_WB_up*deg],color=Color[1],linewidth=lineWidth)
#        plt.plot([x0,x1],[alpha_WB_low*deg, alpha_WB_low*deg],color=Color[1],linewidth=lineWidth)
#        plt.plot([x0,x1],[alpha_Ref*deg, alpha_Ref*deg],color=Color[0],linewidth=lineWidth)
#        plt.plot([x0,x1],arr([slope,slope])*deg,'--k')
##        plt.plot(timeLists[iSim][I]/kyr,slopes[iSim][I]*deg,'-k',linewidth=.5,markersize=.5)
#        
##        plt.plot([0,nx],alphaRef*deg,)
#        
#        plt.plot(np.arctan(locSlope)*180.0/np.pi,'.k')
##        plt.plot(diffTopo)
#        
#        
##        plt.subplot(312)
##        plt.plot(Topo)
###        plt.contourf(plt.contour(phase.T))
##        plt.xlim(x0,x1)
##        
#        
#        plt.subplot(312)
#        plt.hist(locSlope*deg,50)
#        
#        
#        
#        plt.subplot(313)
##        plt.contourf(strain.T,np.linspace(0,25,20))
##        plt.colorbar()
#        
#        
##        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=1, nLayersY=4.00,maxStrain=50.0)
#        lc=2.0e3
#        PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
#        PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
#        
#    
#        PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sampleRate]
##        
#        plt.scatter(PartX,PartY,c=PartStrain,s=pointSize,vmin=0.0,vmax=20.0,edgecolors='None')     
#
##        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
##        
#        ymax = 4.
##        plt.axis([-1.0/aspectRatioDrawing*ymax,0.0,0.0,ymax])
#    #    plt.axis("off")
##        
##        CMAP = arr([[0.5,0.5,.8,1.]])
##        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=True)
##        plt.register_cmap(cmap=CMAP)
##        plt.set_cmap("myColorMap")
##        plt.set_cmap('viridis')
#        plt.set_cmap('inferno')
#
#        plt.axis('off')
#        plt.axis('equal')
#        
#        
#        
##        plt.xlim(0,1200)
#        
#        
#        
#        
#        
#        
#        
#        
        
        
        
        
        
        
        
        
    # end iStep
    
    thisData['time_list'] = timeList
    thisData['xFront'] = xFront
    thisData['xBase'] = xBase
    thisData['xMid'] = xMid
    thisData['locSlopes'] = locSlopes
    thisData['lenPrisms'] = lenPrisms
    
#    slopes.append(slope)
#    timeLists.append(timeList)
#    xFronts.append(xFront)
#    xBases.append(xBase)
#    xMids.append(xMid)
#    outFolders.append(outFolder)
    
    
# end iSim
    
np.save(outFile,Data)

#    
#np.savez(outFile,
#         slopes = slopes,
#         timeLists = timeLists,
#         xFronts=xFronts,
#         xBases=xBases,
#         xMids=xMids,
#         outFolders=outFolders,
#         weakList=weakList,
#         LambadList=LambdaList)