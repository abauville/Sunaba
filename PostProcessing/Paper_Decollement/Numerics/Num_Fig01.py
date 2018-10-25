#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '../CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
import OutputDef as Output
from PaperDecollement_Utils import getColormap, get_XYandPattern


#   Misc
# ============================================
Style = CritTaper_Style.Style()


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


#   Chi, Lambda
# ============================================
chi_list = [ 1, 10, 20, 30, 40, 50, 60, 70, 80]
#chi_list = [ 1, 20, 40, 60, 80]
#chi_list = [ 1, 10, 20, 20, 40, 40, 60, 60, 80]
Lambda_list = [0, 40, 60, 80]
#Lambda_list = [60]
nC = len(chi_list)
nL = len(Lambda_list)
nHor = len(Lambda_list)
nVer = len(chi_list)

#   Fig, Axes
# ============================================
aspectRatio = 0.27
fig  = Figz_Utils.Figure(101,height=21.0,width=29.7,mode='draft')
bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.62,leftMarginPad=1.25,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad=1.5,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')

Axes = Figz_Utils.makeAxes(fig,nVer,nHor,aspectRatio=aspectRatio,leftMarginPad=1.5,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')
axInfo = Axes['info']
cBarAxes = Figz_Utils.makeAxes(fig,1,1,topMarginPad=axInfo['topMarginPad']+axInfo['plotsHeight']*nVer+axInfo['yPad']*(nVer-1)+0.5,bottomMarginPad=1.25,
                               leftMarginPad=axInfo['leftMarginPad']+5.0, rightMarginPad=axInfo['rightMarginPad']+5.0)

#   File stuff
# ============================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"
superDirList = []
i = 0
for iC in range(nC):
    for iL in range(nL):
        superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iC], Lambda_list[iL]))


#   Production Mode
# ============================================
ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 100
    pointSize = sampleRate/100.0



#   Colormap
# ============================================
CMAP, colorList_Type = Style.getCmap_Type()
Type_list = np.linspace(-1.0,2.0,colorList_Type.shape[0])
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")
#plt.colorbar(ticks=[-1,0,1,2])




#   BigAxes
# ============================================
plt.sca(bigAxes['11'])
ax = plt.gca()
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.xlim(-10,70)
plt.xticks([0,20,40,60],[0,40,60,80])

plt.ylim(89.5,-8.5)
plt.yticks([0,10,20,30,40,50,60,70,80])
ax.tick_params(which='both',direction='in')

plt.text(-10,-9.5,'$\\mathbf{\\lambda}$ $\\mathbf{[\\%]}$',size=12)
plt.text(-11.6,-7,'$\\mathbf{\\chi}$ $\\mathbf{[\\%]}$',rotation=90,size=12,verticalAlignment='baseline')


#   Plotting loop
# ============================================
nSim = nC*nL
iSim = 0
for iC in range(nC):
    for iL in range(nL):
        
        IL = np.argmin(np.abs(Lambdas-Lambda_list[iL]/100.0))
        IC = np.argmin(np.abs(chis-chi_list[iC]/100.0))
        
        Type = floatType[IL,IC]
        IT = np.argmin(np.abs(Type_list-Type))
#        print(IC)
        outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
        Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
        
        ymax = 4.25
        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
        plt.axis("off")
        
        CMAP = arr([.0,.0,.0,.0])
        CMAP[:-1] = colorList_Type[IT,:]
        
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False)
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")
        iSim+=1



#   Plotting loop
# ============================================
plt.sca(cBarAxes['11'])
plt.contourf(Type_list,[0.0,1.0],arr([Type_list,Type_list]),Type_list,cmap='custom')
plt.xlim(2,-1)
plt.xticks([-1,0,1,2])
plt.yticks([])
plt.text(0.5,-.25,"Type",weight='bold',verticalAlignment='top',horizontalAlignment='center')
plt.text(1.5,0.45 ,'I'  ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
plt.text(0.5,0.45 ,'II' ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
plt.text(-0.5,0.45,'III',horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)