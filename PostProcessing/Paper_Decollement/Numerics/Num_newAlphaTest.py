#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:47 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '..//CriticalTaper')
sys.path.insert(0, '../')
import OutputDef as Output
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern
from CritTaper_utils import Taper

# Misc
# =========================================
Style = CritTaper_Style.Style()

#   Define chi_list
# =========================================
#thisFile_chi_list = [1, 10, 20, 30, 40, 50, 60, 70]
#thisFile_chi_list = [1, 5, 10, 15, 20, 25, 30, 40]
thisFile_chi_list = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65,70,80]
#thisFile_chi_list = [30]
nC = len(thisFile_chi_list)
nSim = nC


#
##  Load Data of Type and dAlpha
## =========================================
#beta= 0.0 * np.pi/180.0
#Lambda = 0.6
#loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0))
#Lambdas = loadedData["Lambdas"][()]
#chis    = loadedData["chis"][()]
#floatType = loadedData["floatType"][()]
#alphas_diff = loadedData["alphas_diff"][()]
#taper_angles = loadedData["taper_angles"][()]
#


#  Figure
# =========================================
#aspectRatio = 1.0/5.0
#fig             = Figz_Utils.Figure(103,height=29.7,mode='draft')
fig             = Figz_Utils.Figure(103,height=21.0,width=29.7,mode='draft')
#fig             = Figz_Utils.Figure(103,height=29.7)

bottomMarginPad = 2.0

#Axes_back       = Figz_Utils.makeAxes(fig,nSim,1,leftMarginPad=0.0,yPad=0.1)
#AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.075,leftMarginPad=10.0,rightMarginPad=10.0+10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,setAspectRatioBasedOn='y')
#AxesDrawing     = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=1.0/aspectRatio,leftMarginPad=0.1,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesAlpha       = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.75,leftMarginPad=13.0,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 0.25,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesxFault      = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=0.75,leftMarginPad=AxesAlpha['info']['leftMarginPad']+AxesAlpha['info']['plotsWidth']+AxesAlpha['info']['xPad'],rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')

#Axes_back       = Figz_Utils.makeAxes(fig,nSim,1,leftMarginPad=0.0,yPad=0.1,bottomMarginPad=bottomMarginPad)
#AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.075,leftMarginPad=9.5,rightMarginPad=9.5+10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,setAspectRatioBasedOn='y')
#AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.07,leftMarginPad=0.75,rightMarginPad=0.75+10.5,topMarginPad = 2.0,bottomMarginPad = bottomMarginPad,xPad = 0.25,setAspectRatioBasedOn='y')
#AxesDrawing     = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=1.0/aspectRatio,leftMarginPad=2.5,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')

#aspectRatio = 2.5
#AxesAlpha       = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=aspectRatio,leftMarginPad=1.0,rightMarginPad=0.5,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 0.25,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesxFault      = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=aspectRatio,leftMarginPad=AxesAlpha['info']['leftMarginPad']+AxesAlpha['info']['plotsWidth']+AxesAlpha['info']['xPad'],rightMarginPad=0.0,topMarginPad = 0.0,bottomMarginPad = bottomMarginPad,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')

aspectRatio = 0.3
aspectRatioDrawing = .27


yShiftTop = 1.0
#yShift0 = AxesAlpha['info']['plotsHeight']+AxesAlpha['info']['yPad']
yShift = yShiftTop
xPad = 0.25
yPad = 1.5

# =====

yShift0 = AxesDrawing['info']['plotsHeight']+AxesDrawing['info']['yPad']
yShift += yShift0
AxesAlpha       = Figz_Utils.makeAxes(fig,1,int(nSim),aspectRatio=aspectRatio,
                                      leftMarginPad=1.0,rightMarginPad=0.0,
                                      topMarginPad = yShift,bottomMarginPad = bottomMarginPad,
                                      xPad = 0.25,yPad=yPad,
                                      setAspectRatioBasedOn='x')
yShift0 = AxesAlpha['info']['plotsHeight']+AxesAlpha['info']['yPad']
yShift += yShift0
AxesxFault      = Figz_Utils.makeAxes(fig,1,int(nSim),aspectRatio=aspectRatio,
                                      leftMarginPad=1.0,rightMarginPad=0.0,
                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                      xPad=xPad,
                                      setAspectRatioBasedOn='x')


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
    sampleRate = 120
    pointSize = sampleRate/30.0



    
    
    
    
## Figure Alpha
# ============================================
#(nChi, nBeta, nLambda, LambdaRef_list, 
# chi_list, betas_all, alphas_Ref_all, 
# alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=np.linspace(0.0,30.0,13.0)*np.pi/180.0, nChi=61, nLambda=61,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)

Setup = Output.readInput(superRootFolder + superDirList[0] +  '/Output/Input/input.json')

yr = 365.25*24.0*3600.0
kyr = 1000.0 * yr
beta = 0.0
Lambda = 0.6
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/locSlopes_Beta%02d_Lambda%02d_WeakDense.npy" % (beta*10.0*180.0/np.pi,Lambda*100.0)).item(0)



alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)
iSim0 = 0

transparency = 0.2
Color  = arr([[.85,.15,.25],[.25,.5,.5],[.85,.15,.25]])
Color_w_transparency = arr([[.25,.5,.5,transparency],
                            [.85,.15,.25,transparency],
                            [.85,.15,.25,transparency]])


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
    Lambda = 0.6#Lambdas[iSim]
    chi = thisFile_chi_list[iSim]/100.0
#    plt.sca(AxesAlpha['%i1' % (iSim+1)])
    plt.sca(AxesAlpha['1%i' % (iSim+1)])    
    ax = plt.gca()
    ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
    plt.xticks([])
    plt.yticks([])
    
    thisData = loadedData["Lambda%02d_chi%02d" % (Lambda*100,chi*100)]
    locSlopes = thisData["locSlopes"]
    lenPrisms = thisData["lenPrisms"]
    time_list = thisData["time_list"]
    tSteps = thisData["tSteps"]
    xFront = thisData["xFront"]
    xBase = thisData["xBase"]
    xMid = thisData["xMid"]
    
    nSteps = len(tSteps)

        
    
#     =============================================================================
    #                       Create taper and get data
    
    rho_w = 1000.0
    rho = 2500.0
    phiRef   = 30.0*np.pi/180.0
    
#        chi = 1e-7
    
    LambdaRef = Lambda
    #chi_list = [.05,0.7,0.7]
#        beta_list = np.linspace(35.0,-5.0,9)/deg
    beta = 0.0
#        chi = chi
    LambdaWeak = (1.0-chi) * LambdaRef   + chi
    
    ## ============= RefTaper =================    
    tprRef = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaRef+1e-6,
                rho_w=rho_w, rho=rho)
    tprBasalWeak = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaWeak,
                rho_w=rho_w, rho=rho)
    tprTotalWeak = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaWeak,
                rho_w=rho_w, rho=rho)
    tprRef.computeAlphaVsBeta(n=2010)
    tprBasalWeak.computeAlphaVsBeta(n=2010)
    tprTotalWeak.computeAlphaVsBeta(n=2010)
    
    alpha_Ref = tprRef.findAlpha(beta,"average")
    alpha_WB_up = tprBasalWeak.findAlpha(beta,"lower")
    alpha_WB_low = tprBasalWeak.findAlpha(beta,"upper")
    alpha_WF = tprTotalWeak.findAlpha(beta,"average")
    
    
    
    
#        tpr_list.append(tpr)
    
            
    #                       Create taper and get data
    # =============================================================================
    
    

    
#    if plot==0:
#    alpha_Ref = alphas_Ref[iSim]
#    alpha_WF = alphas_WF[iSim]
#    alpha_WB_up = alphas_WB_up[iSim]
#    alpha_WB_low = alphas_WB_low[iSim]
    
    
    ## Compute the colormaps from Histograms
    plt.cla()
    recompute = True
    if recompute:
        n = 45 
        res = 0.25 # in degrees
        N = int(n/res)
        transformedSlopes = np.zeros((len(time_list),N))
        bins_in=np.linspace(0.5*res,n-.5*res,N)
        
        for iStep in range(len(time_list)):
            Hist = np.zeros(N)
            if (iStep%50==0):
                print("iStep = %i/%i" % (iStep,nSteps))
    #        Hist = plt.hist(locSlopes[iStep]*deg,bins=np.linspace(0.0,n,N+1))    
            for i in range(len(locSlopes[iStep])):
                I = np.argmin(np.abs(locSlopes[iStep][i]*deg-bins_in))
                Hist[I] += 1
            transformedSlopes[iStep] = Hist/lenPrisms[iStep]
        plt.cla()
    
    
    flip = False
    
    if not flip:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
#        if iCol == 0:
#            y0 = 0.0#alpha_WB_low*deg*.25
#            y1 = np.max(np.concatenate([alphas_WB_up[iSim:iSim+ncols],[alpha_Ref]]))*1.1*deg
##            plt.ylim([y0,y1])
        y0 = 0.0
        y1 = 30.0
        x0 = time_list[0]/kyr
        x1 = time_list[-1]/kyr

    
        lineWidth = 0.75
        plt.plot([x0,x1],[alpha_WF *deg, alpha_WF*deg],'b',linewidth=lineWidth)
        plt.fill([x0,x1,x1,x0],arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,color=Color[1],alpha=transparency)
        plt.plot([x0,x1],[alpha_WB_up*deg, alpha_WB_up*deg],color=Color[1],linewidth=lineWidth)
        plt.plot([x0,x1],[alpha_WB_low*deg, alpha_WB_low*deg],color=Color[1],linewidth=lineWidth)
        plt.plot([x0,x1],[alpha_Ref*deg, alpha_Ref*deg],color=Color[0],linewidth=lineWidth)
#        plt.contour(timeLists[iSim]/kyr,locSlopes*deg,'-k',linewidth=.5,markersize=.5)
#        plt.plot(timeLists[iSim][I]/kyr,slopes[iSim][I]*deg,'-k',linewidth=.5,markersize=.5)
        
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
#        plt.plot(slopes[iSim][I]*deg,timeLists[iSim][I]/kyr,'-k',linewidth=.5,markersize=.5)
        
    
        plt.ylim([x1,x0])
        plt.xlim([y0,y1])


    TT,BB = np.meshgrid(time_list/kyr,bins_in)    
    TT = TT.T
    BB = BB.T
#    plt.pcolor(TT,BB,transformedSlopes)

    plt.pcolor(TT,BB,transformedSlopes,vmin=0,vmax=0.2)

    
    #   Colormap
    # ============================================
    from matplotlib.colors import LinearSegmentedColormap
    n = 256
    nBeg = int(np.floor(n/2))
    nEnd = n-nBeg
#    Colors = np.array([np.linspace(1.0,0.0,n), # Red
#                       np.linspace(1.0,0.0,n), # Green
#                       np.linspace(1.0,0.0,n), # Blue
#                       np.linspace(1.0,1.0,n)]).T # Alpha
    
#    Colors = np.array([np.concatenate([np.linspace(1.0,0.0,nBeg),np.ones(nEnd)]), # Red
#                       np.concatenate([np.linspace(1.0,0.0,nBeg),np.ones(nEnd)]), # Green
#                       np.concatenate([np.linspace(1.0,0.0,nBeg),np.ones(nEnd)]), # Blue
#                       np.linspace(1.0,1.0,n)]).T # Alpha
    
#    
#    Colors = np.array([np.concatenate([np.linspace(1.0,0.0,nBeg),np.linspace(0.0,0.0,nBeg)]), # Red
#                       np.concatenate([np.linspace(0.0,1.0,nBeg),np.linspace(1.0,0.0,nBeg)]), # Green
#                       np.concatenate([np.linspace(0.0,0.0,nBeg),np.linspace(0.0,1.0,nBeg)]), # Blue
#                       np.linspace(1.0,1.0,n)]).T # Alpha

#    Colors[0,3] = 0.0
    
    Colors = [[1.0,1.0,1.0,1.0],
              [.5,0.5,1.0,1.0],
              [0.0,0.0,1.0,1.0],
              [1.0,.25,.5,1.0],
              [1.0,0.0,0.0,1.0]]
    CMAP = LinearSegmentedColormap.from_list('custom',Colors,N=n)        
    plt.register_cmap(cmap=CMAP)
    plt.set_cmap("custom")
    
#    plt.set_cmap("hot")
#    CMAP = plt.get_cmap()
#    CMAP._lut[0,:] = [1.0,1.0,1.0,1.0]

#    elif plot==1:
        
        
        
        
        
        
        
        
        
        
        
    if not flip:
#        plt.sca(AxesxFault['%i1' % (iSim+1)])
        plt.sca(AxesxFault['1%i' % (iSim+1)])
        ax = plt.gca()
        ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
        plt.xticks([])
        plt.yticks([])
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.plot(time_list/kyr,Setup.Grid.nxC-xFront,'-',color=[.5,.5,.5],linewidth=1.0,markersize=.5)
        plt.plot(time_list/kyr,Setup.Grid.nxC-xMid  ,'-',color=[.6,.1,.0],linewidth=1.0,markersize=.5)
        plt.plot(time_list/kyr,Setup.Grid.nxC-xBase ,'-',color=[.0,.6,.1],linewidth=1.0,markersize=.5)
        
        
    
#        plt.plot(Setup.Grid.nxC-xFronts[iSim][I],timeLists[iSim][I]/kyr,'-k',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xMids[iSim][I]  ,timeLists[iSim][I]/kyr,'-r',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xBases[iSim][I] ,timeLists[iSim][I]/kyr,'-b',linewidth=1.0,markersize=.5)
    
        plt.xlim([x0,x1])
        plt.ylim([0,Setup.Grid.nxC*0.75])
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
    
#        plt.plot(Setup.Grid.nxC-xFront,time_list/kyr,'-k',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xMid  ,time_list/kyr,'-r',linewidth=1.0,markersize=.5)
#        plt.plot(Setup.Grid.nxC-xBase ,time_list/kyr,'-b',linewidth=1.0,markersize=.5)
    
        plt.ylim([x1,x0])
        plt.xlim([0,Setup.Grid.nxC/2.0])


# end iSim
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    