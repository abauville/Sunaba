#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
import OutputDef as Output
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern


nSubs = 5
aspectRatio = 0.5
fig  = Figz_Utils.Figure(104,height=21.0,width=29.7,mode='draft')
Axes = Figz_Utils.makeAxes(fig,nSubs,3,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=2.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=-1.5,setAspectRatioBasedOn='x')

chi_list = [ 1, 25, 80]
Lambda_list = [60,60,60]

nSim = len(chi_list)

tSteps_list = arr([[200,400,600,800,1000],
                   [200,400,600,800,1000],
                   [200,400,600,800,1000]])
nSteps = nSubs

superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select/Beta00/"
superDirList = []
i = 0
iW = 0
for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iSim], Lambda_list[iSim]))



ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 60
    pointSize = sampleRate/100.0



for iSim in range(nSim):
    for iStep in range(nSteps):
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + "Out_%05d/" % tSteps_list[iSim,iStep]
        Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iStep+1,iSim+1)])
        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=16, nLayersY=5.00,maxStrain=5.0)
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
        
        ymax = 9.0
        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
        plt.axis("off")
        
        CMAP = arr([])
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=True)
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")


