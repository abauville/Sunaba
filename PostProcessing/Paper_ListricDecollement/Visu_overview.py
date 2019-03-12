#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
#sys.path.insert(0, '../CriticalTaper')
sys.path.insert(0, '../Utils/')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
import OutputDef as Output
from PaperDecollement_Utils import getColormap, get_XYandPattern





#   Chi, Lambda
# ============================================
Hc_nd_list = [1.0/16.0, 1.0/4.0, 1.0/2.0, 1.0, 2.0]
#Hc_nd = 1.0/1.0
#Hc_nd = 1.0/8.0


Lambda = 0.9
#weakFac = 0.4
PfWeakFac_list = [0.05, 0.1, 0.25, 0.5]
frictionWeakFac = 0.0
cohesionWeakFac_list = [0.1, 0.5, 0.9]
Lambda_b_Fac = 0.0

#
#nC = len(chi_list)
#nL = len(Lambda_list)

column_list = cohesionWeakFac_list
row_list = PfWeakFac_list
fixed_list = Hc_nd_list
iFixed = 2
nCol = len(column_list)
nRow = len(row_list)



#   Fig, Axes
# ============================================
aspectRatio = 0.29
#fig  = Figz_Utils.Figure(101,height=21.0,width=29.7,mode='draft')
fig  = Figz_Utils.Figure(100+iFixed,height=25.0,width=40.7,mode='crop')
#bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.85,leftMarginPad=.75,rightMarginPad=0.25,topMarginPad=1.2,bottomMarginPad=1.5,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')

yPad = 0.0
Axes = {}
topMarginPad = 1.2
for iRow in range(nRow):
#    yPad = yPad_list[iStep]
#    aspectRatio = aspectRatio_list[iRow]
 
#    tempAxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=0,setAspectRatioBasedOn='x')
    tempAxes = Figz_Utils.makeAxes(fig,1,nCol,aspectRatio=aspectRatio,leftMarginPad=1.,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')
    topMarginPad += tempAxes['info']['plotsHeight']+yPad
    
    for iCol in range(nCol):
        Axes['%i%i' % (iRow+1,iCol+1)] = tempAxes['1%i' % (iCol+1)]


#axInfo = tempAxes['info']
#topMarginPad += 0.5
#cBarAxes = Figz_Utils.makeAxes(fig,1,1,topMarginPad=topMarginPad,bottomMarginPad=.9,
#                               leftMarginPad=axInfo['leftMarginPad']+5.0, rightMarginPad=axInfo['rightMarginPad']+5.0)

#   File stuff
# ============================================
superRootFolder = "/Users/abauville/Output/ListricDecollement/LightOutput/Output_Test_Dilation/"
superDirList = []
i = 0

if fixed_list == cohesionWeakFac_list:
    cohesionWeakFac = cohesionWeakFac_list[iFixed]
elif fixed_list == PfWeakFac_list:
    PfWeakFac = PfWeakFac_list[iFixed]
elif fixed_list == Hc_nd_list:
    Hc_nd = Hc_nd_list[iFixed]
    
for iC in range(nCol):
    for iR in range(nRow):
        if column_list == cohesionWeakFac_list:
            cohesionWeakFac = cohesionWeakFac_list[iC]
        elif column_list == PfWeakFac_list:
            PfWeakFac = PfWeakFac_list[iC]
        elif column_list == Hc_nd_list:
            Hc_nd = Hc_nd_list[iC]
            
        if row_list == cohesionWeakFac_list:
            cohesionWeakFac = cohesionWeakFac_list[iR]
        elif row_list == PfWeakFac_list:
            PfWeakFac = PfWeakFac_list[iR]
        elif row_list == Hc_nd_list:
            Hc_nd = Hc_nd_list[iR]
            
        
        superDirList.append( "Lambda%02d_Hc%03d_CW%02d_PfW%02d_GFac005/" % (Lambda*100, Hc_nd*100, cohesionWeakFac*100, PfWeakFac*100) )


#   Production Mode
# ============================================
ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/100.0
else:
    sampleRate = 5
    pointSize = sampleRate/15.0



##   Colormap
## ============================================
#CMAP, colorList_Type = Style.getCmap_Type()
#Type_list = np.linspace(1.0,4.0,colorList_Type.shape[0])
#plt.register_cmap(cmap=CMAP)
#plt.set_cmap("custom")
##plt.colorbar(ticks=[-1,0,1,2])



#
##   BigAxes
## ============================================
#plt.sca(bigAxes['11'])
#ax = plt.gca()
#ax.xaxis.tick_top()
#ax.xaxis.set_label_position('top')
#ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
##plt.xlim(-10,70)
##plt.xticks([0,20,40,60],[0,40,60,80])
#plt.xlim(-10,50)
#plt.xticks([0,20,40],[0,33,66])
#
##plt.ylim(89.5,-8.5)
##plt.yticks([0,10,20,30,40,50,60,70,80])
#dum=-.03
#plt.ylim(np.sum(aspectRatio_list)-dum,0.0-dum)
#plt.yticks(np.cumsum(aspectRatio_list),chi_list)
#ax.tick_params(which='both',direction='in')
#
#plt.text(-10.0,-.01,'$\\mathbf{\\lambda^*}$ $\\mathbf{[\\%]}$',size=12)
#plt.text(-12.0,+.09,'$\\mathbf{\\chi}$ $\\mathbf{[\\%]}$',rotation=90,size=12,verticalAlignment='baseline')
#


#   Plotting loop
# ============================================
nSim = nCol*nRow
iSim = 0
for iC in range(nCol):
    for iR in range(nRow):
        
#        IL = np.argmin(np.abs(Lambdas-Lambda_ov_list[iR]/100.0))
#        IC = np.argmin(np.abs(chis-chi_list[iC]/100.0))
        
#        Type = floatType[IL,IC]
#        IT = np.argmin(np.abs(Type_list-Type))
#        print(IC)
        folderList = os.listdir(superRootFolder + superDirList[iSim] + "Output/")
        folderList.sort()
        outFolder = folderList[-1]
#        outFolder = 'Out_%05d' % tSteps_list[iSim]
        dataFolder = superRootFolder + superDirList[iSim] + "Output/" + outFolder + "/"
        Char = Output.readInput(superRootFolder + superDirList[iSim] + 'Output/Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iR+1,iC+1)])
        
#        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,minStrain=1.0,maxStrain=5.0)
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=12, xmin=-24.0, xmax=0.0, ymin=0.0,ymax=1.00,nLayersY=5,minStrain=2.0,maxStrain=10.0,mainDir='x')
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
        
        
#        aspectRatio = aspectRatio_list[iC]
        xmin = -20.0
        x0 = xmin; x1 = 0.0
        y0 = 0.0 ; y1 = -1.0*aspectRatio*xmin

        plt.axis([xmin,0.0,0.0,-1.0*aspectRatio*xmin])
        plt.axis('off')
        
#        CMAP = arr([.4,.6,.8,.0])
#        CMAP[:-1] = colorList_Type[IT,:]
#        
#        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,0.0,1.0,0.0])
#        plt.register_cmap(cmap=CMAP)
#        plt.set_cmap("myColorMap")
        CMAP=([[]])
#        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,.0,.95,.0],
#                           RGBShift=[[0.0, 0.0, 0.0], 
#                                     [0.0, 0.0, 0.0], 
#                                     [-.2,-.2,0.2], 
#                                     [0.0, 0.0, 0.0]])
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,1.0,1.0,1.0],
                           RGBShift=[[0.0, 0.0, 0.0], 
                                     [0.0, 0.0, 0.0], 
                                     [-.2,-.2,0.2], 
                                     [0.0, 0.0, 0.0]])
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")
        
#        if (iR==0):
#            pyMax = np.max(PartY)
#            pyPad = 0.3
#            aspectRatio_list[iC] = ((pyMax+pyPad)/(-xmin))
#            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
        
        
        iSim+=1


#
##   Plotting loop
## ============================================
#plt.sca(cBarAxes['11'])
#plt.contourf(Type_list,[0.0,1.0],arr([Type_list,Type_list]),Type_list,cmap='custom')
#plt.xlim(1,4)
#plt.xticks([1,2,3,4])
#plt.yticks([])
#plt.text(2.5,-.4,"M",weight='bold',verticalAlignment='top',horizontalAlignment='center')
#plt.text(1.5,0.4 ,'I'  ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
#plt.text(2.5,0.4 ,'II' ,horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)
#plt.text(3.5,0.4,'III',horizontalAlignment='center',verticalAlignment='center',family='Times New Roman',color='white',size=16)