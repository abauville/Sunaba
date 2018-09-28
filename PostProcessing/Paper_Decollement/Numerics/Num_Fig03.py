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
chi_list = [1, 10, 20, 30, 40, 50, 60, 70]
#chi_list = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
nC = len(chi_list)
nSim = nC



#  Load Data of Type and dAlpha
# =========================================
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType.npz")
Lambdas = loadedData["Lambdas"][()]
chis    = loadedData["chis"][()]
floatType = loadedData["floatType"][()]
alphas_diff = loadedData["alphas_diff"][()]
taper_angles = loadedData["taper_angles"][()]



#  Figure
# =========================================
aspectRatio = 1.0/3.0
fig             = Figz_Utils.Figure(104,height=29.7,mode='draft')
AxesDAlpha      = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.1,leftMarginPad=1.0,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,setAspectRatioBasedOn='y')
#AxesAlpha       = Figz_Utils.makeAxes(fig,nSubs,1,aspectRatio=2.0,leftMarginPad=4.5,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')
#AxesxFault      = Figz_Utils.makeAxes(fig,nSubs,1,aspectRatio=2.0,leftMarginPad=7.5,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')
AxesDrawing     = Figz_Utils.makeAxes(fig,nSim,1,aspectRatio=1.0/aspectRatio,leftMarginPad=8.5,rightMarginPad=10.5,topMarginPad = 0.0,bottomMarginPad = 0.0,xPad = 1.0,yPad = 0.1,setAspectRatioBasedOn='y')



#   File system
# =========================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta50/"
superDirList = []
i = 0

for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda60" % (chi_list[iSim]))




#  Production mode
# =========================================
ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 60
    pointSize = sampleRate/100.0




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
    
    ymax = 5.0
    plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
#    plt.axis("off")
    
    CMAP = arr([])
    CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False)
    plt.register_cmap(cmap=CMAP)
    plt.set_cmap("myColorMap")





#  Plot dAlpha
# =========================================
alpha_diff_list = np.zeros(nSim)
alphas_diff_bar = alphas_diff/taper_angles

for iSim in range(nSim):
    chi = chi_list[iSim]/100.0
    Lambda = 0.6
    iC = np.argmin(np.abs(chis-chi))
    iL = np.argmin(np.abs(Lambdas-Lambda))
    alpha_diff_list[iSim] = alphas_diff_bar[iL,iC]

    
max_alphaDiff = np.max(alphas_diff_bar[iL,:])
chi_max_alphaDiff = chis[np.argmax(alphas_diff_bar[iL,:])]
plt.sca(AxesDAlpha['11'])
plt.plot(alphas_diff_bar[iL,:],chis*100.0,'-b')
#plt.plot(max_alphaDiff,chi_max_alphaDiff*100.0,'ok')
plt.plot(alpha_diff_list,chi_list,'ob')

plt.plot([0.0,0.0],chi_list[slice(0,nC,nC-1)],'--k')

x0 = -0.45
x1 = 0.45
plt.xlim([x0,x1])

#plt.contourf(Lambdas*100.0,chis*100.0,floatType,np.linspace(-1.0001,2.0001,130),vmin=-1.00,vmax=2.00)
#plt.plot([60.0,60.0],[0.0,100.0],'--k')
plt.contourf(arr([x0,x1]),chis*100.0,arr([floatType[iL,:],floatType[iL,:]]).T,np.linspace(-1.0001,2.0001,130),vmin=-1.00,vmax=2.00)
#plt.plot(floatType[iL,:],chis*100.0)
plt.ylim([75.0,-5.0])
#   Colormap
# ============================================
CMAP = Style.getCmap_Type()
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")
#plt.colorbar(ticks=[-1,0,1,2])
    
    
    
    
#   Text
# ============================================
plt.text(x0+(x1-x0)*0.15 ,12,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(x0+(x1-x0)*0.15 ,40,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(x0+(x1-x0)*0.15 ,65,'III',family='Times New Roman',color='w',size=28,weight='bold')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    