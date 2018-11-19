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

#chi_list = [ 5, 15, 25, 35, 45, 55, 65]
#chi_list = [ 1, 5, 10, 15, 20, 25, 30, 35, 40]
#
#chi_list = [ 1, 20, 40, 60, 80]
#chi_list = [ 1, 10, 20, 20, 40, 40, 60, 60, 80]
Lambda_list = [40, 60, 80]

tSteps_list = [426, 1728, 264, 574, 1946, 260, 610, 1989, 241, 528, 1720, 215, 525, 1413, 189, 444, 1296, 165, 372, 1075, 142, 280, 883, 118, 239, 696, 96]

aspectRatio_list = arr([0.31984994, 0.34588737, 0.35917899, 0.34664289, 0.33786318,
       0.28453353, 0.26216297, 0.23113716, 0.20000444])


#Lambda_list = [60,60,60]
#Lambda_list = [60]
nC = len(chi_list)
nL = len(Lambda_list)
nCol = len(Lambda_list)
nRow = len(chi_list)



#   Fig, Axes
# ============================================
aspectRatio = 0.29
#fig  = Figz_Utils.Figure(101,height=21.0,width=29.7,mode='draft')
fig  = Figz_Utils.Figure(101,height=20.5,width=21.0,mode='production')
bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.82,leftMarginPad=.75,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad=1.5,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')

yPad = 0.0
Axes = {}
topMarginPad = 1.5
for iRow in range(nRow):
#    yPad = yPad_list[iStep]
    aspectRatio = aspectRatio_list[iRow]
 
#    tempAxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=0,setAspectRatioBasedOn='x')
    tempAxes = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio,leftMarginPad=1.,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')
    topMarginPad += tempAxes['info']['plotsHeight']+yPad
    
    for iCol in range(nCol):
        Axes['%i%i' % (iRow+1,iCol+1)] = tempAxes['1%i' % (iCol+1)]


axInfo = tempAxes['info']
topMarginPad += 0.5
cBarAxes = Figz_Utils.makeAxes(fig,1,1,topMarginPad=topMarginPad,bottomMarginPad=.9,
                               leftMarginPad=axInfo['leftMarginPad']+5.0, rightMarginPad=axInfo['rightMarginPad']+5.0)

#   File stuff
# ============================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select2/Beta00/"
superDirList = []
i = 0
for iC in range(nC):
    for iL in range(nL):
        superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iC], Lambda_list[iL]))


#   Production Mode
# ============================================
ProductionMode = True
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 300
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
#plt.xlim(-10,70)
#plt.xticks([0,20,40,60],[0,40,60,80])
plt.xlim(-10,50)
plt.xticks([0,20,40],[40,60,80])

#plt.ylim(89.5,-8.5)
#plt.yticks([0,10,20,30,40,50,60,70,80])
dum=-.03
plt.ylim(np.sum(aspectRatio_list)-dum,0.0-dum)
plt.yticks(np.cumsum(aspectRatio_list),chi_list)
ax.tick_params(which='both',direction='in')

plt.text(-10.0,-.01,'$\\mathbf{\\lambda}$ $\\mathbf{[\\%]}$',size=12)
plt.text(-12.0,+.09,'$\\mathbf{\\chi}$ $\\mathbf{[\\%]}$',rotation=90,size=12,verticalAlignment='baseline')



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
#        outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
        outFolder = 'Out_%05d' % tSteps_list[iSim]
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
        Char = Output.readInput(superRootFolder + superDirList[iSim] + '/Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,minStrain=1.0,maxStrain=5.0)
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
        
        
        aspectRatio = aspectRatio_list[iC]
        xmin = -15.75
        x0 = xmin; x1 = 0.0
        y0 = 0.0 ; y1 = -1.0*aspectRatio*xmin

        plt.axis([xmin,0.0,0.0,-1.0*aspectRatio*xmin])
        plt.axis('off')
        
        CMAP = arr([.0,.0,.0,.0])
        CMAP[:-1] = colorList_Type[IT,:]
        
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,0.0,1.0,0.0])
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")
        
        
        if (iL==0):
            pyMax = np.max(PartY)
            pyPad = 0.3
            aspectRatio_list[iC] = ((pyMax+pyPad)/(-xmin))
            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
        
        
        iSim+=1



#   Plotting loop
# ============================================
plt.sca(cBarAxes['11'])
plt.contourf(Type_list,[0.0,1.0],arr([Type_list,Type_list]),Type_list,cmap='custom')
plt.xlim(2,-1)
plt.xticks([-1,0,1,2])
plt.yticks([])
plt.text(0.5,-.25,"Type",weight='bold',verticalAlignment='top',horizontalAlignment='center')
plt.text(1.5,0.4 ,'I'  ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
plt.text(0.5,0.4 ,'II' ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
plt.text(-0.5,0.4,'III',horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)