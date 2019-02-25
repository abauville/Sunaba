#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '../CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
#import CritTaper_dataMaker
import Figz_Utils
#import CritTaper_Style
import OutputDef as Output
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern
#

#tSteps_list = arr([[214, 428, 641, 855, 1069, 1282, 1496, 1710],
#                   [248, 495, 743, 990, 1238, 1486, 1733, 1981],
#                   [134, 269, 403, 538,  672,  806,  941, 1075]])


nWeak = 2
chi_list = arr([20,20,20,20,60,60,60,60])
tSteps_list = arr([766,766,766,766,638,638,638,638]) 

res_list = [25,50,100,200,25,50,100,200]
#nSim = len(res_list)
nSim = 4
#nSteps = tSteps_list.shape[1]
aspectRatio = 0.5
fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='production')
#fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='draft')
#fig  = Figz_Utils.Figure(1040,height=7.0,width=29.7)

yPad_list = [.0, .0, .0, .0, .0, .0, .0]
#aspectRatio_list = arr([.2, .25, .3, .3, .4, .4, .45, .45])*.75

aspectRatio_list = arr([0.25,.25,.25,0.25,.25,.25,.25,.25])

#aspectRatio_list =arr([0.36334528])

topMarginPad = 0.2
Axes = {}
yPad = 0


backgroundBoxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=2.09,leftMarginPad=0.25/2.0,rightMarginPad=0.25/2.0,topMarginPad=.0,bottomMarginPad = 0.0,xPad = 0.5/2.0,yPad=0,setAspectRatioBasedOn='x')

plt.sca(backgroundBoxes['11'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
plt.axis([0,1,0,1])
plt.sca(backgroundBoxes['12'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.8,.8,.88,1.0])
plt.axis([0,1,0,1])
#plt.sca(backgroundBoxes['13'])
#plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
#plt.axis([0,1,0,1])

backgroundBoxes['11'].axis('off')
backgroundBoxes['12'].axis('off')
#backgroundBoxes['13'].axis('off')

for iSim in range(nSim):
#    yPad = yPad_list[iStep]
    aspectRatio = aspectRatio_list[iSim]
 
    tempAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=0,setAspectRatioBasedOn='x')
    topMarginPad += tempAxes['info']['plotsHeight']+yPad
    
    Axes['%i1' % (iSim+1)] = tempAxes['11']
    Axes['%i2' % (iSim+1)] = tempAxes['12']
#    Axes['%i3' % (iSim+1)] = tempAxes['13']
    
    






superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/ResTest/"
superDirList = []
i = 0
iW = 0
for iWeak in range(nWeak):
    for iSim in range(nSim):
        superDirList.append("Weak%02d/sameDispLim/Res_%03d" % (chi_list[iWeak*nSim+iSim], res_list[iWeak*nSim+iSim]))



ProductionMode = True
#ProductionMode = True
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/60.0
    pointSize_list  = [sampleRate/1.0,sampleRate/3.0,sampleRate/14.0,sampleRate/60.0]
else:
    sampleRate = 100
    pointSize = sampleRate/60.0
    pointSize_list  = [sampleRate/2.0,sampleRate/5.0,sampleRate/12.0,sampleRate/60.0]
for iWeak in range(nWeak):
#    if iSim == 0:
#        pointSize = sampleRate/12.0

    for iSim in range(nSim):
        dataFolder = superRootFolder + superDirList[iWeak*nSim+iSim] + "/Output/" + "Out_%05d/" % tSteps_list[iWeak*nSim+iSim]
        Char = Output.readInput(superRootFolder + superDirList[iWeak*nSim+iSim] + '/Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        print(dataFolder)
        PartX = []
        PartY = []
        PartPattern = []
        
        pointSize = pointSize_list[iSim]
        ax = plt.sca(Axes["%i%i" % (iSim+1,iWeak+1)])
        xmin = -16.0
        x0 = xmin; x1 = 0.0
        aspectRatio = aspectRatio_list[iWeak*nSim+iSim]
        y0 = 0.0 ; y1 = -1.0*aspectRatio*xmin



#        # BW Version        
#        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=1, xmin=-4.0, xmax=0.0, ymin=0.0,ymax=1.15,nLayersY=6,minStrain=2.0,maxStrain=10.0,mainDir='y')
#        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
#        if (iSim==1):
#            pyMax = np.max(PartY)
#            pyPad = 0.56
#            aspectRatio_list[iStep] = ((pyMax+pyPad)/(-xmin))
#            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
#            
#        CMAP = arr([[0.6,0.6,0.6,1.0],
#                    [0.0,0.0,0.0,1.0],
#                    [0.8,0.8,0.8,1.0],
#                    [0.0,0.0,0.0,1.0],
#                    [1.0,1.0,1.0,1.0],
#                    [0.0,0.0,0.0,1.0]])
#        
#        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,1.0,1.0,1.0])
#        plt.register_cmap(cmap=CMAP)
#        plt.set_cmap("myColorMap")
        
        # Color Version        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=12, xmin=-24.0, xmax=0.0, ymin=0.0,ymax=1.00,nLayersY=5,minStrain=2.0,maxStrain=10.0,mainDir='x')
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
#        if (iSim==1):
#            pyMax = np.max(PartY)
#            pyPad = 0.56
#            aspectRatio_list[iStep] = ((pyMax+pyPad)/(-xmin))
#            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
#            
#        CMAP = arr([[0.6,0.6,0.6,1.0],
#                    [0.0,0.0,0.0,1.0],
#                    [0.8,0.8,0.8,1.0],
#                    [0.0,0.0,0.0,1.0],
#                    [1.0,1.0,1.0,1.0],
#                    [0.0,0.0,0.0,1.0]])
        CMAP=([[]])
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,.0,.95,.0],
                           RGBShift=[[0.0, 0.0, 0.0], 
                                     [0.0, 0.0, 0.0], 
                                     [-.2,-.2,0.2], 
                                     [0.0, 0.0, 0.0]])
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")
        
        
        
        plt.axis([xmin,0.0,0.0,-1.0*aspectRatio*xmin])
        plt.axis('off')

        Letters = 'ABCDEFGH'
        res = [16,32,64,128]
#        if iStep<nSteps-1:
        plt.fill([x0,x0,x0+1.8,x0+1.8],[y0,y0+.65,y0+.65,y0],'k')
        plt.text(x0+0.03,y0+.14,'%s. H/%i'% (Letters[iWeak*nSim+iSim],res[i]),color='w')

        if iSim==0:
            titleText1 = ['Type=1\n',
                          'Type=2\n',
                          'Type=3\n']
            titleText2 = ['Pure accretion',
                         'Accretion and underplating',
                         'Pure underplating']
            titleText2 = ['Style 2',
                          'Style 3']
    #            plt.text(x0+(x1-x0)/2.0,y0+1.8,titleText1[iSim],horizontalAlignment='center',weight='bold')
            plt.text(x0+(x1-x0)*.00,y0+2.5,titleText2[iWeak],horizontalAlignment='left',weight='bold')