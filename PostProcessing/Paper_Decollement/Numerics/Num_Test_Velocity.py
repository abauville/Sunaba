#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:22:22 2018

@author: abauville
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
import OutputDef as Output
from PaperDecollement_Utils import getColormap, get_XYandPattern



#chi_list = [ 1, 10, 20, 30, 40, 50, 60, 70, 80]
#chi_list = [ 1, 20, 40, 60, 80]
chi_list = [80]
Lambda_list = [60]
nC = len(chi_list)
nL = len(Lambda_list)
nHor = len(Lambda_list)
nVer = len(chi_list)


aspectRatio = 0.3
fig  = Figz_Utils.Figure(101,height=21.0,width=29.7,mode='draft')
bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.66,leftMarginPad=1.25,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')
ax = plt.gca()
#plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
#plt.axis([.0,100.0,100.0,.0])

Axes = Figz_Utils.makeAxes(fig,nVer,nHor,aspectRatio=aspectRatio,leftMarginPad=1.5,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')


superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"
superDirList = []
i = 0
for iC in range(nC):
    for iL in range(nL):
        superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iC], Lambda_list[iL]))



ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 1
    pointSize = sampleRate/3.0


nSim = nC*nL
iSim = 0
for iC in range(nC):
    for iL in range(nL):
        outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
        Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
        
        CMAP = arr([ [0.4,0.5,0.8,0.0] ])
#        CMAP = arr([ [0.7,0.7,0.9,0.0] ])
        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False,strainDarknessFactor=0.0)
        plt.register_cmap(cmap=CMAP)
        plt.set_cmap("myColorMap")
        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')
        if iSim == 0:
            raw = Output.getData(dataFolder + 'Vx.bin',True)
            xmin = raw.xmin; ymin = raw.ymin
            xmax = raw.xmax; ymax = raw.ymax
            x = np.linspace(raw.xmin,raw.xmax,raw.nx)
            y = np.linspace(raw.ymin,raw.ymax,raw.ny-1)
            dx = x[1]-x[0]
            X, Y = np.meshgrid(x,y)
            X = X.T
            Y = Y.T
        Vx = Output.getData(dataFolder + 'Vx.bin',True).data
        Vx = 0.5*(Vx[:,0:-1] + Vx[:,1:])
        Vy = Output.getData(dataFolder + 'Vy.bin',True).data
        Vy = 0.5*(Vy[0:-1,:] + Vy[1:,:])
        phase = Output.getData(dataFolder + 'phase.bin',True).data
        phase = 0.5*(phase[:,0:-1] + phase[:,1:])
        phase = 0.5*(phase[0:-1,:] + phase[1:,:])
        
        Vx = Vx*phase
        Vy = Vy*phase
        
        Vi = np.sqrt(Vx**2+Vy**2)
        
        vmin = -1.0
        vmax = 1.0
        Vc = Vx.copy()
        Vc[Vx<vmin] = vmin
        Vc[Vx>vmax] = vmax
        
#        start_points_y = np.linspace(0.0,0.9,50)
#        start_points=arr([ (xmin+dx)*np.ones(len(start_points_y)) , start_points_y]).T        
#        plt.streamplot(x[0::rx],y[0::ry],Vx[0::rx,0::ry].T,Vy[0::rx,0::ry].T,color='r',density=2,arrowsize=1,start_points=start_points,minlength=0.9)
#        start_points_y = arr([0.6,0.7,0.8])#np.linspace(0.6,0.9,100)
#        start_points=arr([ (xmin+dx)*np.ones(len(start_points_y)) , start_points_y]).T
        n = 250
        I = (np.random.rand(n)*PartX.size).astype(np.int)
        start_points = arr([  PartX[I], PartY[I] ]).T
        
        colors = arr([[255,255,  0],
                      [255,255,255],
                      [  0,255,255]]) / 255.0
        CMAP = LinearSegmentedColormap.from_list('custom2',colors,N=nTot)       
        plt.register_cmap(cmap=CMAP)
        plt.streamplot(x,y,Vx.T,Vy.T,color=Vc.T, density=2.5,arrowsize=1,maxlength=0.15,start_points=start_points,linewidth=0.75,cmap='custom2')        
#        plt.streamplot(x,y,Vx.T,Vy.T,color=Vc.T, density=1.5,arrowsize=1)     
        
#        plt.set_cmap('custom2')
        plt.colorbar()
#        plt.set_cmap("seismic")
#        plt.colorbar()
        
        ymax = 5.0
        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
        plt.axis("off")
        
        
        
        
        iSim+=1
