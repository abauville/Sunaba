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



chi_list = [ 1, 20, 60]
tSteps_list = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
                     [249, 497, 746, 994, 1243, 1492, 1740, 1989],
                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 



#
#chi_list = [ 1, 15, 60]
#tSteps_list = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
#                     [256, 513, 769, 1026, 1282, 1538, 1795, 2051],
#                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 

#
#tSteps_list = arr([[ 1544],
#                   [ 1981],
#                   [ 1075]])

nSteps = tSteps_list.shape[1]
aspectRatio = 0.5
fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='crop')
#fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='draft')
#fig  = Figz_Utils.Figure(1040,height=7.0,width=29.7)

yPad_list = [.0, .0, .0, .0, .0, .0, .0]
#aspectRatio_list = arr([.2, .25, .3, .3, .4, .4, .45, .45])*.75
aspectRatio_list = arr([0.14767239, 0.19295576, 0.21157659, 0.23101366, 0.27265011,
       0.28095144, 0.31765885, 0.35418436])
aspectRatio_list = arr([0.18007018, 0.21782023, 0.2332407 , 0.25981279, 0.28025861,
       0.30212507, 0.30531708, 0.33758197])

#aspectRatio_list =arr([0.36334528])

topMarginPad = 0.2
Axes = {}
yPad = 0


backgroundBoxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=2.09,leftMarginPad=0.25/2.0,rightMarginPad=0.25/2.0,topMarginPad=.0,bottomMarginPad = 0.0,xPad = 0.5/2.0,yPad=0,setAspectRatioBasedOn='x')

plt.sca(backgroundBoxes['11'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
plt.axis([0,1,0,1])
plt.sca(backgroundBoxes['12'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.8,.8,.88,1.0])
plt.axis([0,1,0,1])
plt.sca(backgroundBoxes['13'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
plt.axis([0,1,0,1])

backgroundBoxes['11'].axis('off')
backgroundBoxes['12'].axis('off')
backgroundBoxes['13'].axis('off')

for iStep in range(nSteps):
#    yPad = yPad_list[iStep]
    aspectRatio = aspectRatio_list[iStep]
 
    tempAxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=0,setAspectRatioBasedOn='x')
    topMarginPad += tempAxes['info']['plotsHeight']+yPad
    
    Axes['%i1' % (iStep+1)] = tempAxes['11']
    Axes['%i2' % (iStep+1)] = tempAxes['12']
    Axes['%i3' % (iStep+1)] = tempAxes['13']
    
    

Lambda_list = [60,60,60]

nSim = len(chi_list)



superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select2/Beta00/"
superDirList = []
i = 0
iW = 0
for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iSim], Lambda_list[iSim]))



ProductionMode = True
#ProductionMode = True
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/65.0
else:
    sampleRate = 500
    pointSize = sampleRate/60.0

for iSim in range(nSim):
    if iSim == 0:
        pointSize = sampleRate/12.0

    for iStep in range(nSteps):
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + "Out_%05d/" % tSteps_list[iSim,iStep]
        Char = Output.readInput(superRootFolder + superDirList[iSim] + '/Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iStep+1,iSim+1)])
        xmin = -16.0
        x0 = xmin; x1 = 0.0
        aspectRatio = aspectRatio_list[iStep]
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
        if (iSim==1):
            pyMax = np.max(PartY)
            pyPad = 0.56
            aspectRatio_list[iStep] = ((pyMax+pyPad)/(-xmin))
            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
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

        Letters = 'ABC'
        if iStep<nSteps-1:
            plt.fill([x0,x0,x0+1.4,x0+1.4],[y0,y0+.75,y0+.75,y0],'k')
            plt.text(x0+0.03,y0+.14,'%s.$t_{%i}$' % (Letters[iSim],(iStep+1)),color='w')
        else:
            plt.fill([x0,x0,x0+1.9,x0+1.9],[y0,y0+.75,y0+.75,y0],'k')
            plt.text(x0+0.03,y0+.14,'$%s.t_{ref}$' % (Letters[iSim]),color='w')

        if iStep==0:
            titleText1 = ['Type=1\n',
                          'Type=2\n',
                          'Type=3\n']
#            titleText2 = ['Pure accretion',
#                         'Accretion and underplating',
#                         'Pure underplating']
            titleText2 = ['Style 1',
                          'Style 2',
                          'Style 3']
#            plt.text(x0+(x1-x0)/2.0,y0+1.8,titleText1[iSim],horizontalAlignment='center',weight='bold')
            plt.text(x0+(x1-x0)*.00,y0+2.25,titleText2[iSim],horizontalAlignment='left',weight='bold',size=13)