#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 14:04:05 2018

@author: abauville
"""


# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan, arctan
from PaperDecollement_Utils import getColormap, get_XYandPattern
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
from matplotlib.colors import Normalize
import CritTaper_Style
#from scipy import interpolate
#from scipy import ndimage
import CritTaper_dataMaker
from Units import *
from numpy import array as arr
import Figz_Utils


# Load Data
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/slopes.npz")
slopes = loadedData["slopes"][()]
timeLists = loadedData["timeLists"][()]
xFronts = loadedData["xFronts"][()]
xBases = loadedData["xBases"][()]
xMids = loadedData["xMids"][()]
weakList = arr([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])*100.0
weakList_short = arr([1, 5, 10, 20, 40, 60, 80])
#weakList = arr([1, 5])
nW = len(weakList)
LambdaList = arr([00,40,60,80])
nL = len(LambdaList)

## Get Data from CritTaper theory analalysis
#Style = CritTaper_Style.Style()
#CMAP_ref = plt.get_cmap(Style.colormap)._lut

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False,nChi=21, nBeta=21, nLambda = 21)
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all 



## Create the folder tree
# ================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/wWater/Beta00/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    
ProductionMode = False


if ProductionMode:
    sampleRate = 1
    pointSize = 0.0002
else:
    sampleRate = 100
    pointSize = 0.1

superDirList = []
i = 0
iL = 0
Lambdas = np.zeros(nW*nL)
chis = np.zeros(nW*nL)
for Lambda in LambdaList:
    iW = 0
    for weak in weakList:
        superDirList.append("Weak%02d/Lambda%02d" % (weak, Lambda))
#        IRefColorMap[i] = IRefColorMap_Big[iL,iW]
        Lambdas[i] = Lambda/100.0
        chis[i] = weak/100.0
        i+=1
        iW+=1
    iL+=1




rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]


#rootFolders = [''] * len(superDirList) # initialize to the right size
#nStepsList = np.zeros(len(superDirList),dtype='int16');
#for i in range(len(superDirList)):
#    rootFolders[i] = superRootFolder + superDirList[i] + "/Output/"
#
#    
#    stepFolderList = os.listdir(rootFolders[i])
#    try:
#        stepFolderList.remove('.DS_Store')
#    except ValueError:
#        Dummy=0
#    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char

pushVel = 10.0*cm/yr



plot = 0 # 0 angle, 1 xFront
Ax = []
nrows = nL
ncols = nW
if plot==0:
#    fig    = Figz_Utils.Figure(2,width=29.7,height=21.0,mode='draft')
    fig    = Figz_Utils.Figure(2,width=29.7,height=21.0,mode='production')
elif plot==1:
#    fig    = Figz_Utils.Figure(3,width=29.7,height=21.0,mode='draft')
    fig    = Figz_Utils.Figure(3,width=29.7,height=21.0,mode='production')
else:
    raise ValueError('unknown plot type')
    
#fig2   = Figz_Utils.Figure(3,height=29.7,mode='draft')
#fig    = Figz_Utils.Figure(6,height=13.0,mode='draft')
#AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
#    AxesDum['12'].axis('off')
Axes   = Figz_Utils.makeAxes(fig,nrows,ncols,aspectRatio=1.5,leftMarginPad=1.5,rightMarginPad=1.0,topMarginPad = 0.0,xPad = 0.1,yPad = 1.0)
#Axes_xFront = Figz_Utils.makeAxes(fig,nrows,ncols,aspectRatio=0.8,leftMarginPad=1.5,rightMarginPad=1.0,topMarginPad = 0.0)
for i in range(nrows):
    for j in range(ncols):
        Ax.append(Axes['%i%i' % (i+1,j+1)])

 
  


nSim = len(superDirList)

#    nSim = 7
iSim0 = 0#nSim-1
Hsed = 2.0 * km
nSteps = 747#nStepsList[-1]
i0 = 746#nSteps-1#jump-2
#    iStep = i0
jump = 1000000
frame = 0

iCol = 0
iRow = 0
#slopes = []#np.zeros(nSim)
#timeLists = []

alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)

for iSim in range(iSim0,nSim):

    
    Lambda = Lambdas[iSim]
    chi = chis[iSim]
    
#    if chi in weakList_short:
    
    IC = np.argmin(np.abs(chi_list-chi))
    IL = np.argmin(np.abs(LambdaRef_list-Lambda))
    deg = 180.0/np.pi
    beta = 0.0
    alphas_Ref[iSim] = Taper_Ref[IL].findAlpha(beta,"average")
    alphas_WF[iSim] = Taper_WF[IL*nChi+IC].findAlpha(beta,"average")
    alphas_WB_up[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"upper")
    alphas_WB_low[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"lower")

    
    
    
    
for iSim in range(iSim0,nSim):  
    ## Plot stuff
    I = np.all([xBases[iSim]>0,xFronts[iSim]>0,xMids[iSim]>0],axis=0)
    Lambda = Lambdas[iSim]
    chi = chis[iSim]
    plt.sca(Ax[iSim])
    ax = plt.gca()
    plt.xticks([])
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if iCol == 0:
        plt.ylabel('$\\lambda=%i$'%(Lambda*100))
    else: 
        plt.yticks([])
        
    if iRow == 0:
        plt.xlabel('$\\chi=%i$'%(chi*100))
        ax.xaxis.set_label_position('top')
    

    
    if plot==0:
        alpha_Ref = alphas_Ref[iSim]
        alpha_WF = alphas_WF[iSim]
        alpha_WB_up = alphas_WB_up[iSim]
        alpha_WB_low = alphas_WB_low[iSim]
        
        if iCol == 0:
            y0 = 0.0#alpha_WB_low*deg*.25
            y1 = np.max(np.concatenate([alphas_WB_up[iSim:iSim+ncols],[alpha_Ref]]))*1.1*deg
#            plt.ylim([y0,y1])
        
        x0 = timeLists[iSim][0]/kyr
        x1 = timeLists[iSim][-1]/kyr

        
        plt.plot([x0,x1],[alpha_Ref*deg, alpha_Ref*deg],'r')
        plt.plot([x0,x1],[alpha_WF *deg, alpha_WF*deg],'b')
        plt.fill([x0,x1,x1,x0],arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,'g',alpha=0.1)
        plt.plot([x0,x1],[alpha_WB_up*deg, alpha_WB_up*deg],'g')
        plt.plot([x0,x1],[alpha_WB_low*deg, alpha_WB_low*deg],'g')
        
        
    #    plt.plot(timeLists[iSim]/kyr,slopes[iSim]*deg,'ok',markersize=1.0)
        plt.plot(timeLists[iSim][I]/kyr,slopes[iSim][I]*deg,'-k',linewidth=.5,markersize=.5)
        
        
        plt.xlim([x0,x1])
#        y0 = alpha_WB_low*deg*.75
#        y1 = np.max([alpha_WB_up,alpha_Ref])*1.1*deg
##        y0 = 0.0
##        y1 = 35.0
        plt.ylim([y0,y1])

    elif plot==1:
#        xFronts[iSim][xFronts[iSim]==0] = Setup.Grid.nxC
#        xMids[iSim][xMids[iSim]==0] = Setup.Grid.nxC
#        xBases[iSim][xBases[iSim]==0] = Setup.Grid.nxC
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xFronts[iSim][I],'-k',linewidth=1.0,markersize=.5)
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xMids[iSim][I]  ,'-r',linewidth=1.0,markersize=.5)
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xBases[iSim][I] ,'-b',linewidth=1.0,markersize=.5)


    iCol += 1
    if iCol>=ncols:
        iCol = 0
        iRow+=1

# end iSim

        
dpi = 300
if plot==0:
    plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/Angles",dpi=dpi)
elif plot==1:
    plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/FrontPos",dpi=dpi)












