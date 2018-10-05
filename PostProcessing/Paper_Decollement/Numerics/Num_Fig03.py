#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:47 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../')
import OutputDef as Output
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern


# Misc
# =========================================
Style = CritTaper_Style.Style()

#   Define chi_list
# =========================================
thisFile_chi_list = [1, 10, 20, 30, 40, 50, 60, 70]
#thisFile_chi_list = [1, 5, 10, 15, 20, 25, 30, 40]
#chi_list = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
nC = len(thisFile_chi_list)
nSim = nC



#  Load Data of Type and dAlpha
# =========================================
beta= 0.0 * np.pi/180.0
Lambda = 0.6
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0))
Lambdas = loadedData["Lambdas"][()]
chis    = loadedData["chis"][()]
floatType = loadedData["floatType"][()]
alphas_diff = loadedData["alphas_diff"][()]
taper_angles = loadedData["taper_angles"][()]



#  Figure
# =========================================
aspectRatio = 1.0/3.0
#fig             = Figz_Utils.Figure(104,height=29.7,mode='draft')
fig             = Figz_Utils.Figure(103,height=29.7)

bottomMarginPad = 2.0

#Axes_back       = Figz_Utils.makeAxes(fig,nSim,1,leftMarginPad=0.0,yPad=0.1)
#AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.075,leftMarginPad=10.0,rightMarginPad=10.0+10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,setAspectRatioBasedOn='y')
#AxesDrawing     = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=1.0/aspectRatio,leftMarginPad=0.1,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesAlpha       = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.75,leftMarginPad=13.0,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 0.25,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesxFault      = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.75,leftMarginPad=AxesAlpha['info']['leftMarginPad']+AxesAlpha['info']['plotsWidth']+AxesAlpha['info']['xPad'],rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')

Axes_back       = Figz_Utils.makeAxes(fig,nSim,1,leftMarginPad=0.0,yPad=0.1,bottomMarginPad=bottomMarginPad)
#AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.075,leftMarginPad=9.5,rightMarginPad=9.5+10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,setAspectRatioBasedOn='y')
AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.07,leftMarginPad=0.75,rightMarginPad=0.75+10.5,topMarginPad = 2.0,bottomMarginPad = bottomMarginPad,xPad = 0.25,setAspectRatioBasedOn='y')
AxesDrawing     = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=1.0/aspectRatio,leftMarginPad=2.5,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')
AxesAlpha       = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.9,leftMarginPad=12.0,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 0.25,yPad = 0.1,setAspectRatioBasedOn='y')
AxesxFault      = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.9,leftMarginPad=AxesAlpha['info']['leftMarginPad']+AxesAlpha['info']['plotsWidth']+AxesAlpha['info']['xPad'],rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')



for iSim in range(0,nSim,2):
    plt.sca(Axes_back['%i1' % (iSim+1)])
    ax = plt.gca()
    ax.patch.set_facecolor([.94,.94,.99])
    
for iSim in range(nSim):
    plt.sca(Axes_back['%i1' % (iSim+1)])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks([])
    plt.yticks([])
    

#   File system
# =========================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta%02d/" % round(beta*180.0/np.pi*10.0)
superDirList = []
i = 0

for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda60" % (thisFile_chi_list[iSim]))




#  Production mode
# =========================================
ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/30.0
else:
    sampleRate = 60
    pointSize = sampleRate/30.0




#  Plot wedge drawings
# =========================================
for iSim in range(nSim):
    outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
    dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
    Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
    timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
    
    PartX = []
    PartY = []
    PartPattern = []
    
    ax = plt.sca(AxesDrawing["%i1" % (iSim+1)])
    
    PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
    plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
    
    ymax = 3.5
    plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
#    plt.axis("off")
    
    CMAP = arr([])
    CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False)
    plt.register_cmap(cmap=CMAP)
    plt.set_cmap("myColorMap")

    plt.axis('off')



#  Plot dAlpha
# =========================================
plt.sca(AxesDAlpha['11'])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xlabel('$\\mathbf{\\bar{\\Delta \\alpha}}$')
plt.ylabel('$\\mathbf{\\chi}$ $\\mathbf{[\\%]}$')

alpha_diff_list = np.zeros(nSim)
alphas_diff_bar = alphas_diff/taper_angles

for iSim in range(nSim):
    chi = thisFile_chi_list[iSim]/100.0
    Lambda = 0.6
    iC = np.argmin(np.abs(chis-chi))
    iL = np.argmin(np.abs(Lambdas-Lambda))
    alpha_diff_list[iSim] = alphas_diff_bar[iL,iC]

    
max_alphaDiff = np.max(alphas_diff_bar[iL,:])
chi_max_alphaDiff = chis[np.argmax(alphas_diff_bar[iL,:])]


plt.plot(alphas_diff_bar[iL,:],chis*100.0,'-b')
#plt.plot(max_alphaDiff,chi_max_alphaDiff*100.0,'ok')
plt.plot(alpha_diff_list,thisFile_chi_list,'ob')

plt.plot([0.0,0.0],thisFile_chi_list[slice(0,nC,nC-1)],'--k')

x0 = -0.5
x1 = 0.5
plt.xlim([x0,x1])
plt.xticks([-.5,0.0,.5])

#plt.contourf(Lambdas*100.0,chis*100.0,floatType,np.linspace(-1.0001,2.0001,130),vmin=-1.00,vmax=2.00)
#plt.plot([60.0,60.0],[0.0,100.0],'--k')
plt.contourf(arr([x0,x1]),chis*100.0,arr([floatType[iL,:],floatType[iL,:]]).T,np.linspace(-1.0001,2.0001,130),vmin=-1.00,vmax=2.00)
#plt.plot(floatType[iL,:],chis*100.0)
plt.ylim(arr([72.0,-2.0]))
#   Colormap
# ============================================
CMAP, color_list = Style.getCmap_Type()
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")
#plt.colorbar(ticks=[-1,0,1,2])
    
    
    
    
#   Text
# ============================================
plt.text(x0+(x1-x0)*0.15 ,12,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(x0+(x1-x0)*0.15 ,40,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(x0+(x1-x0)*0.15 ,65,'III',family='Times New Roman',color='w',size=28,weight='bold')
    
    
    
    
    
    
    
    
    
    
## Figure Alpha
# ============================================
(nChi, nBeta, nLambda, LambdaRef_list, 
 chi_list, betas_all, alphas_Ref_all, 
 alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
 Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=np.linspace(0.0,30.0,13.0)*np.pi/180.0, nChi=61, nLambda=61,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)

Setup = Output.readInput(superRootFolder + superDirList[0] +  '/Output/Input/input.json')

yr = 365.25*24.0*3600.0
kyr = 1000.0 * yr

loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/slopes_Beta%02d_Lambda%02d.npz" % (beta*10.0*180.0/np.pi,Lambda*100.0))
slopes = loadedData["slopes"][()]
timeLists = loadedData["timeLists"][()]
xFronts = loadedData["xFronts"][()]
xBases = loadedData["xBases"][()]
xMids = loadedData["xMids"][()]


alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)
iSim0 = 0




for iSim in range(iSim0,nSim):

    
    Lambda = 0.6#Lambdas[iSim]
    chi = thisFile_chi_list[iSim]/100.0
    
    IC = np.argmin(np.abs(chi_list-chi))
    IL = np.argmin(np.abs(LambdaRef_list-Lambda))
    deg = 180.0/np.pi
    beta = 0.0
    alphas_Ref[iSim] = Taper_Ref[IL].findAlpha(beta,"average")
    alphas_WF[iSim] = Taper_WF[IL*nChi+IC].findAlpha(beta,"average")
    alphas_WB_up[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"upper")
    alphas_WB_low[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"lower")


plot=0
for iSim in range(iSim0,nSim):  
    ## Plot stuff
    I = np.all([xBases[iSim]>0,xFronts[iSim]>0,xMids[iSim]>0],axis=0)
    I = np.arange(len(timeLists[iSim]))
    Lambda = Lambdas[iSim]
    chi = chis[iSim]
    plt.sca(AxesAlpha['%i1' % (iSim+1)])
    ax = plt.gca()
    ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
    plt.xticks([])
    plt.yticks([])
    
    

        
    

    
#    if plot==0:
    alpha_Ref = alphas_Ref[iSim]
    alpha_WF = alphas_WF[iSim]
    alpha_WB_up = alphas_WB_up[iSim]
    alpha_WB_low = alphas_WB_low[iSim]
    
    flip = False
    
    if not flip:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
#        if iCol == 0:
#            y0 = 0.0#alpha_WB_low*deg*.25
#            y1 = np.max(np.concatenate([alphas_WB_up[iSim:iSim+ncols],[alpha_Ref]]))*1.1*deg
##            plt.ylim([y0,y1])
        y0 = 0.0
        y1 = 22.0
        x0 = timeLists[iSim][0]/kyr
        x1 = timeLists[iSim][-1]/kyr

    
        plt.plot([x0,x1],[alpha_Ref*deg, alpha_Ref*deg],'r')
        plt.plot([x0,x1],[alpha_WF *deg, alpha_WF*deg],'b')
        plt.fill([x0,x1,x1,x0],arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,'g',alpha=0.1)
        plt.plot([x0,x1],[alpha_WB_up*deg, alpha_WB_up*deg],'g')
        plt.plot([x0,x1],[alpha_WB_low*deg, alpha_WB_low*deg],'g')
        plt.plot(timeLists[iSim][I]/kyr,slopes[iSim][I]*deg,'-k',linewidth=.5,markersize=.5)
    
        plt.xlim([x0,x1])
        plt.ylim([y0,y1])
        
    else:
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        y0 = 0.0
        y1 = 22.0
        x0 = timeLists[iSim][0]/kyr
        x1 = timeLists[iSim][-1]/kyr
        plt.plot([alpha_Ref*deg, alpha_Ref*deg],[x0,x1],'r')
        plt.plot([alpha_WF *deg, alpha_WF*deg],[x0,x1],'b')
        plt.fill(arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,[x0,x1,x1,x0],'g',alpha=0.1)
        plt.plot([alpha_WB_up*deg, alpha_WB_up*deg],[x0,x1],'g')
        plt.plot([alpha_WB_low*deg, alpha_WB_low*deg],[x0,x1],'g')
        plt.plot(slopes[iSim][I]*deg,timeLists[iSim][I]/kyr,'-k',linewidth=.5,markersize=.5)
        
    
        plt.ylim([x1,x0])
        plt.xlim([y0,y1])

#    elif plot==1:
        
    if not flip:
        plt.sca(AxesxFault['%i1' % (iSim+1)])
        ax = plt.gca()
        ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
        plt.xticks([])
        plt.yticks([])
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xFronts[iSim][I],'-k',linewidth=1.0,markersize=.5)
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xMids[iSim][I]  ,'-r',linewidth=1.0,markersize=.5)
        plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xBases[iSim][I] ,'-b',linewidth=1.0,markersize=.5)
    
#        plt.plot(Setup.Grid.nxC-xFronts[iSim][I],timeLists[iSim][I]/kyr,'-k',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xMids[iSim][I]  ,timeLists[iSim][I]/kyr,'-r',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xBases[iSim][I] ,timeLists[iSim][I]/kyr,'-b',linewidth=1.0,markersize=.5)
    
        plt.xlim([x0,x1])
        plt.ylim([0,Setup.Grid.nxC/2.0])
    else:
        plt.sca(AxesxFault['%i1' % (iSim+1)])
        ax = plt.gca()
        ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
        plt.xticks([])
        plt.yticks([])
        
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    #    plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xFronts[iSim][I],'-k',linewidth=1.0,markersize=.5)
    #    plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xMids[iSim][I]  ,'-r',linewidth=1.0,markersize=.5)
    #    plt.plot(timeLists[iSim][I]/kyr,Setup.Grid.nxC-xBases[iSim][I] ,'-b',linewidth=1.0,markersize=.5)
    
        plt.plot(Setup.Grid.nxC-xFronts[iSim][I],timeLists[iSim][I]/kyr,'-k',linewidth=1.0,markersize=.5)
        plt.plot(Setup.Grid.nxC-xMids[iSim][I]  ,timeLists[iSim][I]/kyr,'-r',linewidth=1.0,markersize=.5)
        plt.plot(Setup.Grid.nxC-xBases[iSim][I] ,timeLists[iSim][I]/kyr,'-b',linewidth=1.0,markersize=.5)
    
        plt.ylim([x1,x0])
        plt.xlim([0,Setup.Grid.nxC/2.0])


# end iSim
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    