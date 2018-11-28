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
chi_list = [1,20,60]
Y1 = [19,24,19]
Y0 = [0, 0, 0]
Y1Ref = 24



#chi_list = [10,40,80]
#chi_list = [1,10,20,40,60,80]
#Y1 = [20,24,20,24,24,24]
#Y0 = [0, 0, 0, 0, 0, 0]
nC = len(chi_list)
nSim = nC

Production = True
recompute = False
if Production:
    res = 0.125/4.0 # in degrees
else:
    res = 0.5 # in degrees


#tSteps_list_Fig02 = arr([[214, 428, 641, 855, 1069, 1282, 1496, 1710],
#                         [248, 495, 743, 990, 1238, 1486, 1733, 1981],
#                         [134, 269, 403, 538,  672,  806,  941, 1075]]) 
tSteps_list_Fig02 = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
                     [249, 497, 746, 994, 1243, 1492, 1740, 1989],
                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 
xticks = np.zeros(tSteps_list_Fig02.shape)




#  Get colors of Type
# =========================================
beta= 0.0 * np.pi/180.0
Lambda = 0.6
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0))
Lambdas = loadedData["Lambdas"][()]
chis    = loadedData["chis"][()]
floatType = loadedData["floatType"][()]
alphas_diff = loadedData["alphas_diff"][()]
taper_angles = loadedData["taper_angles"][()]

CMAP, colorList_Type = Style.getCmap_Type()
Type_list = np.linspace(1.0,4.0,colorList_Type.shape[0])

ITs = np.zeros(nSim,np.int)
Types = np.zeros(nSim)
for iSim in range(nSim):
    Lambda = 0.6#Lambdas[iSim]
    chi = chi_list[iSim]/100.0
    IL = np.argmin(np.abs(Lambdas-Lambda))
    IC = np.argmin(np.abs(chis-chi))
    
    Type = floatType[IL,IC]
    Types[iSim] = Type
    ITs[iSim] = np.argmin(np.abs(Type_list-Type))
    
    


#  Figure
# =========================================
#fig             = Figz_Utils.Figure(104,height=29.7,width=21.0,mode='draft')
fig             = Figz_Utils.Figure(104,height=29.7,width=21.0,mode='production')

bottomMarginPad = 2.0

aspectRatio = 0.3
nCol = 1

#aspectRatio = 0.6
#nCol = 2



yShiftTop = 0.5
#yShift0 = AxesAlpha['info']['plotsHeight']+AxesAlpha['info']['yPad']
yShift = yShiftTop
xPad = 0.25
yPad = 0.5

# =====

yShift0 = 0.0
yShift += yShift0
AxesAlpha = {}
AxesxFault = {}
AxesNotes = {}
nRow = int(np.ceil(nSim/nCol))
yShiftList = [yShift]
for iRow in range(nRow):
    aspectRatio = 0.25*Y1[iRow]/Y1Ref

    AxesTemp1        = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio,
                                          leftMarginPad=0.75,rightMarginPad=0.0,
                                          topMarginPad = yShift,bottomMarginPad = bottomMarginPad,
                                          xPad = 0.25,yPad=yPad,
                                          setAspectRatioBasedOn='x')
    
    AxesTemp2       = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio,
                                          leftMarginPad=0.75,rightMarginPad=0.0,
                                          topMarginPad=yShift-.5,bottomMarginPad = bottomMarginPad,
                                          xPad=xPad,yPad=yPad,
                                          setAspectRatioBasedOn='x')
    AxesTemp3       = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio,
                                      leftMarginPad=0.75,rightMarginPad=0.0,
                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                      xPad=xPad,yPad=yPad,
                                      setAspectRatioBasedOn='x')
    
    
    yShift += AxesTemp1['info']['plotsHeight']+yPad#AxesAlpha['info']['yPad']
    for iCol in range(nCol):
        AxesAlpha  ['%i%i' % (iRow+1,iCol+1)] = AxesTemp1['1%i' % (iCol+1)]
        AxesxFault ['%i%i' % (iRow+1,iCol+1)] = AxesTemp2['1%i' % (iCol+1)]
        AxesNotes  ['%i%i' % (iRow+1,iCol+1)] = AxesTemp3['1%i' % (iCol+1)]



yShift += yPad # extra y padding
AxesLegend       = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio*.4,
                                      leftMarginPad=0.75,rightMarginPad=0.0,
                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                      xPad=xPad,yPad=yPad,
                                      setAspectRatioBasedOn='x')
yShift += yPad*.5+.15
AxesColorbar       = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio*.45,
                                  leftMarginPad=12.5,rightMarginPad=0.5,
                                  topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                  xPad=xPad,yPad=yPad,
                                  setAspectRatioBasedOn='x')



#yShift0 = AxesAlpha['info']['plotsHeight']+0.5#AxesAlpha['info']['yPad']
#yShift += yShift0
yShift=yShiftTop
#yShift += -0.4
#AxesxFault      = Figz_Utils.makeAxes(fig,nRow,nCol,aspectRatio=aspectRatio,
#                                      leftMarginPad=0.75,rightMarginPad=0.0,
#                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
#                                      xPad=xPad,yPad=yPad,
#                                      setAspectRatioBasedOn='x')

#yShift = 0.0
#yShift -= -0.4






for iRow in range(nRow):
    
    
    AxesAlpha['%i1' % (iRow+1)].spines['right'].set_visible(False)
    AxesAlpha['%i1' % (iRow+1)].spines['top'].set_visible(False)
    
    if nCol==2:
        AxesAlpha['%i2' % (iRow+1)].set_yticklabels([])
        AxesAlpha['%i2' % (iRow+1)].spines['right'].set_visible(False)
        AxesAlpha['%i2' % (iRow+1)].spines['top'].set_visible(False)
for iRow in range(nRow-1):
    AxesAlpha['%i1' % (iRow+1)].set_xticklabels([])
    if nCol==2:
        AxesAlpha['%i2' % (iRow+1)].set_xticklabels([])


#   File system
# =========================================
beta = 0.0
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta%02d/" % round(beta*180.0/np.pi*10.0)
superDirList = []
i = 0

for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda60" % (chi_list[iSim]))




    
    
    
    
## Figure Alpha
# ============================================

Setup = Output.readInput(superRootFolder + superDirList[0] +  '/Output/Input/input.json')

yr = 365.25*24.0*3600.0
kyr = 1000.0 * yr
beta = 0.0
Lambda = 0.6
DataFolder = "/Users/abauville/Output/Paper_Decollement/Figz/Data/"
loadedData = np.load(DataFolder + "locSlopes_Beta%02d_Lambda%02d_WeakDense_winSize64.npy" % (beta*10.0*180.0/np.pi,Lambda*100.0)).item(0)



alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)
iSim0 = 0

transparency = 0.2
#Color  = arr([[.85,.15,.25],[.25,.5,.5],[.85,.15,.25]])
#Color_w_transparency = arr([[.25,.5,.5,transparency],
#                            [.85,.15,.25,transparency],
#                            [.85,.15,.25,transparency]])
#    
Color  = arr([Style.colorRef,Style.colorBW,Style.colorFW])
Color_w_transparency = arr([np.concatenate([Style.colorRef,[Style.alphaBW]]),
                            np.concatenate([Style.colorBW,[Style.alphaBW]]),
                            np.concatenate([Style.colorFW,[Style.alphaBW]])])
    
#
#
#for iSim in range(iSim0,nSim):
#
#    
#    Lambda = 0.6#Lambdas[iSim]
#    chi = chi_list[iSim]/100.0
#    
#    IC = np.argmin(np.abs(chi_list-chi))
#    IL = np.argmin(np.abs(LambdaRef_list-Lambda))
#    deg = 180.0/np.pi
#    alphas_Ref[iSim] = Taper_Ref[IL].findAlpha(beta,"average")
#    alphas_WF[iSim] = Taper_WF[IL*nChi+IC].findAlpha(beta,"average")
#    alphas_WB_up[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"upper")
#    alphas_WB_low[iSim] = Taper_WB[IL*nChi+IC].findAlpha(beta,"lower")

deg = 180.0/np.pi
plot=0
for iSim in range(iSim0,nSim):  
    print("iSim = %i/%i" % (iSim,nSim))
    ## Plot stuff
    Lambda = 0.6#Lambdas[iSim]
    chi = chi_list[iSim]/100.0
    iRow = int(np.ceil((iSim+1)/nCol))
    iCol = iSim%nCol + 1
    plt.sca(AxesAlpha['%i%i' % (iRow,iCol)])    
    
    thisData = loadedData["Lambda%02d_chi%02d" % (Lambda*100,chi*100)]
    locSlopes = thisData["locSlopes"]
    lenPrisms = thisData["lenPrisms"]
    time_list = thisData["time_list"]
    tSteps = thisData["tSteps"]
    xFront = thisData["xFront"]
    xBase = thisData["xBase"]
    xMid = thisData["xMid"]
    
    nSteps = len(tSteps)

    xticks[iSim,:] = time_list[tSteps_list_Fig02[iSim,:]]/kyr
        
    
# =============================================================================
#                       Create taper and get data    
    rho_w = 1000.0
    rho = 2500.0
    phiRef   = 30.0*np.pi/180.0
    
    LambdaRef = Lambda
    beta = 0.0
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
#                       Create taper and get data
# =============================================================================




    
# =============================================================================
#                           Compute Intensity
    n = 45 
    th = 0.005
    N = int(n/res)
    bins_in=np.linspace(0.5*res,n-.5*res,N)    
    if not recompute:
        try: 
            Intensity3 = np.load(DataFolder + "Alpha_Lambda%02d_chi%02d_Intensity_res%5f.npy" % (Lambda*100,chi*100,res))
        except FileNotFoundError:
            print("File not found: " + DataFolder + "Alpha_Lambda%02d_chi%02d_Intensity_res%5f.npy" % (Lambda*100,chi*100,res) + ". Recomputing the values")
            recompute = True
    
    if recompute:
        Intensity = np.zeros((len(time_list),N))
        Intensity2 = np.zeros((len(time_list),N))
        Intensity3 = np.zeros((len(time_list),N))
        meanSlopes = np.zeros(len(time_list))
        meanSlopes2 = np.zeros(len(time_list))
        IBack = 64
        
        for iStep in range(len(time_list)):
            Hist = np.zeros(N)
            if (iStep%50==0):
                print("iStep = %i/%i" % (iStep,nSteps))
    #        Hist = plt.hist(locSlopes[iStep]*deg,bins=np.linspace(0.0,n,N+1))    
            for i in range(len(locSlopes[iStep])):
                if iStep>100 and i>len(locSlopes[iStep])-IBack and locSlopes[iStep][i]*deg<2.0:
                    doNothing=1.0
                else:
                    I = np.argmin(np.abs(locSlopes[iStep][i]*deg-bins_in))
                    Hist[I] += 1
            Intensity[iStep] = Hist/lenPrisms[iStep]
            
            I = Intensity[iStep]>th
            meanSlopes2[iStep] = np.mean(bins_in[I])

            th2=0.1
            I = Intensity[iStep]>th2
            
            
            
            Intensity2[iStep] = Intensity[iStep].copy()
            winSize = 5
            for i in range(winSize, N-winSize):
                Intensity2[iStep,i] = np.mean(Intensity2[iStep,i-winSize:i+winSize])
                
            
            I = Intensity2[iStep]<th
            Intensity3[iStep] = Intensity2[iStep].copy()
#            Intensity3[iStep,I] = 0.0
            
        #end iStep            
        np.save(DataFolder + "Alpha_Lambda%02d_chi%02d_Intensity_res%5f.npy" % (Lambda*100,chi*100, res),Intensity3)        
#                           Compute Intensity
# =============================================================================
    



# =============================================================================
#                           Plot
    y0 = Y0[iSim]
    y1 = Y1[iSim]
    x0 = time_list[0]/kyr
#    x1 = time_list[-1]/kyr
    x1 = time_list[tSteps_list_Fig02[iSim,-1]]/kyr


    lineWidth = 0.75
    plt.fill([x0,x1,x1,x0],arr([alpha_WB_up, alpha_WB_up, alpha_WB_low, alpha_WB_low])*deg,color=Color[1],alpha=transparency)
    plt.plot([x0,x1],[alpha_WB_up*deg, alpha_WB_up*deg],color=Color[1],linewidth=lineWidth)
    plt.plot([x0,x1],[alpha_WB_low*deg, alpha_WB_low*deg],color=Color[1],linewidth=lineWidth)
    plt.plot([x0,x1],[alpha_Ref*deg, alpha_Ref*deg],color=Color[0],linewidth=lineWidth)
    plt.plot([x0,x1],[alpha_WF *deg, alpha_WF*deg],color=Color[2,:],linewidth=lineWidth)
    
    plt.xlim([x0,x1])
    plt.ylim([y0,y1])

    TT,BB = np.meshgrid(time_list/kyr,bins_in)    
    TT = TT.T
    BB = BB.T


    vmax = res/3.0
    Intensity3[Intensity3>vmax] = vmax
    plt.contourf(TT,BB,Intensity3,np.linspace(0.0,vmax,32),vmin=0,vmax=vmax)
    
    
    
#                           Plot
# =============================================================================

    

# ============================================
#                 Colormap
    from matplotlib.colors import LinearSegmentedColormap
    n = 256
    nBeg = int(np.floor(n/2))
    nEnd = n-nBeg


    Colors = [[1.0,1.0,1.0,0.0],
              [0.5,0.5,0.5,0.3],
              [0.2,0.2,0.2,0.6],
              [1.0,0.0,0.5,1.0],
              [1.0,1.0,0.2,1.0]]
    

    CMAP = LinearSegmentedColormap.from_list('custom',Colors,N=n)        
    plt.register_cmap(cmap=CMAP)
    plt.set_cmap("custom")

#                 Colormap
# ============================================


# =============================================================================
#                           xFaults

    plt.sca(AxesxFault['%i%i' % (iRow,iCol)])    
    ax = plt.gca()
    ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
    plt.xticks([])
    plt.yticks([])
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.plot(time_list/kyr,Setup.Grid.nxC-xMid  ,'-',color=[.0,.0,.0],linewidth=1.0,markersize=.5)


    plt.xlim([x0,x1])
    plt.ylim([0,Setup.Grid.nxC*0.75*3.0])
    
    
    plt.axis('off')
    
#                          xFaults
# =============================================================================

# end iSim
    
    

    
    
# =============================================================================
#                           Annotations  
Letters='ABCDEF'
iSim = 0
YR = (arr([.0,.0,.0,4.0,3.0,3.0])-y0)/(y1-y0)
for iRow in range(nRow):
    for iCol in range(nCol):

        plt.sca(AxesNotes['%i%i' % (iRow+1,iCol+1)])
        plt.axis([.0,1.0,.0,1.0])
        plt.axis('off')
#        yr = .92
#        xr = .003
        yr = .0
        xr = .003
#        yr = YR[iSim]
        
#        plt.fill(xr+arr([.0,.18,.18,.0]),yr+arr([.0,.0,.1,.1]),color=colorList_Type[ITs[iSim],:],linestyle='None')
        if Types[iSim]>-0.05:
            plt.text(xr+0.005,yr+0.03,Letters[iSim] + '. Type=%i' % np.abs(np.round(Types[iSim])),weight='bold')
        else:
            plt.text(xr+0.005,yr+0.03,Letters[iSim] + '. Type=%.1f' %       (Types[iSim]),weight='bold')
    
        if iRow==0:
            plt.text(.75,.325,'fault@y=0.8')
        
        iSim+=1
  
# x, y labels
plt.sca(AxesNotes['11'])
plt.text(-.04,.9,'$\mathbf{\\alpha [°]}$',fontsize=12,rotation=90)
plt.sca(AxesNotes['31'])
plt.text(0.01,-0.15,'time []',weight='bold',fontsize=12)      

iSim=0
for iRow in range(nRow):
    for iCol in range(nCol):
        plt.sca(AxesAlpha['%i%i' % (iRow+1,iCol+1)])
        plt.xticks(xticks[iSim,:])
        if iRow==0 and iCol==0:
            plt.yticks([0,5,10,15],[0,5,10,''])
            
        iSim+=1

xticklabels = []
for iStep in range(len(tSteps_list_Fig02[0,:])):
    I = tSteps_list_Fig02[2,iStep]
    if iStep<len(tSteps_list_Fig02[0,:])-1:
        xticklabels.append('$t_{%i}$'%(iStep+1))
    else:
        xticklabels.append('$t_{ref}$')
plt.gca().set_xticklabels(xticklabels)

#                           Annotations     
# =============================================================================
 
    

# =============================================================================
#                           Legend  

plt.sca(AxesLegend['11'])
plt.axis([.0,1.0,.0,1.0])
plt.axis('off')
plt.text(.65/2.0,.87,'Wedge stability domains',horizontalAlignment='center',weight='bold')

# stability domain of the basally weakened wedge
fx0 = .015
fx1 = fx0+.18
fy0 = .15
fy1 = .7
fyText = .375
plt.fill([fx0,fx1,fx1,fx0],[fy0,fy0,fy1,fy1],color=Color[1],alpha=transparency,lineWidth=0.0)
plt.plot([fx0,fx1],[fy0,fy0],color=Color[1],linewidth=lineWidth)
plt.plot([fx0,fx1],[fy1,fy1],color=Color[1],linewidth=lineWidth)
plt.text(fx0+0.5*(fx1-fx0),fyText,'basally weakened',horizontalAlignment='center',verticalAlignment='center')


# stability domain of the intact wedge
fx0 = fx1+.015
fx1 = fx0+.18
fy0 = .25
#fy1 = .75
plt.plot([fx0,fx1],[fy0,fy0],color=Color[0],linewidth=lineWidth)
plt.text(fx0+0.5*(fx1-fx0),fyText,'intact',horizontalAlignment='center')

# stability domain of the fully weakened wedge
fx0 = fx1+.015
fx1 = fx0+.18
#fy0 = .5
plt.plot([fx0,fx1],[fy0,fy0],color=Color[2],linewidth=lineWidth)
plt.text(fx0+0.5*(fx1-fx0),fyText,'fully weakened',horizontalAlignment='center')



#                           Legend    
# =============================================================================    
    
    
    
    
# =============================================================================
#                           ColorMap  

plt.sca(AxesColorbar['11'])
cbar = plt.colorbar(ax=AxesAlpha['11'],cax=AxesColorbar['11'],orientation='horizontal')
plt.text(0.5,1.45,'Intensity [%]',horizontalAlignment='center')
cbar.set_ticks([0.0,vmax])
cbar.set_ticklabels(['0','%.0f' % np.round(vmax*100.0)])


#                           ColorMap    
# =============================================================================
 
    
    
    