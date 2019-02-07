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



#chi_list = [ 1, 20, 60]
#tSteps_list = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
#                     [249, 497, 746, 994, 1243, 1492, 1740, 1989],
#                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 




chi_list = [ 1, 20, 60]

Lambda_list = [60,60,60]

nSim = len(chi_list)
#tSteps_list = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
#                     [256, 513, 769, 1026, 1282, 1538, 1795, 2051],
#                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 

tSteps_list = arr([[1728,1728],
                   [1989,1989],
                   [1075,1075]]) 

#
#tSteps_list = arr([[ 1544],
#                   [ 1981],
#                   [ 1075]])

nSteps = tSteps_list.shape[1]
aspectRatio = 0.5
#fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='production')
fig  = Figz_Utils.Figure(102,height=29.7,width=21.0,mode='production')
#fig  = Figz_Utils.Figure(102,height=21.0,width=29.7,mode='draft')
#fig  = Figz_Utils.Figure(1040,height=7.0,width=29.7)

yPad_list = [.0, .0, .0, .0, .0, .0, .0]
#aspectRatio_list = arr([.2, .25, .3, .3, .4, .4, .45, .45])*.75
aspectRatio_list = arr([0.33758197,0.33758197,0.33758197])

#aspectRatio_list =arr([0.36334528])

topMarginPad = 0.2
Axes = {}
yPad = 0


backgroundBoxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=2.09,leftMarginPad=0.25/2.0,rightMarginPad=0.25/2.0,topMarginPad=.0,bottomMarginPad = 0.0,xPad = 0.5/2.0,yPad=0,setAspectRatioBasedOn='x')

plt.sca(backgroundBoxes['11'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
plt.axis([0,1,0,1])
#plt.sca(backgroundBoxes['12'])
#plt.fill([0,0,1,1],[0,1,1,0],color=[.8,.8,.88,1.0])
#plt.axis([0,1,0,1])
plt.sca(backgroundBoxes['12'])
plt.fill([0,0,1,1],[0,1,1,0],color=[.9,.9,.97,1.0])
plt.axis([0,1,0,1])

backgroundBoxes['11'].axis('off')
backgroundBoxes['12'].axis('off')
#backgroundBoxes['13'].axis('off')

#for iStep in range(nSteps):
for iSim in range(nSim):
#    yPad = yPad_list[iStep]
    aspectRatio = aspectRatio_list[iSim]
 
    tempAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=aspectRatio,leftMarginPad=0.25,rightMarginPad=0.25,topMarginPad=topMarginPad,bottomMarginPad = 0.0,xPad = 0.5,yPad=0,setAspectRatioBasedOn='x')
    topMarginPad += tempAxes['info']['plotsHeight']+yPad
    
    Axes['%i1' % (iSim+1)] = tempAxes['11']
    Axes['%i2' % (iSim+1)] = tempAxes['12']
#    Axes['%i3' % (iStep+1)] = tempAxes['13']
    
    




superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select2/Beta00/"
superDirList = []
i = 0
iW = 0
for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iSim], Lambda_list[iSim]))



ProductionMode = False
#ProductionMode = True
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/60.0
else:
    sampleRate = 250#50
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
        
        ax = plt.sca(Axes["%i%i" % (iSim+1,iStep+1)])
        xmin = -16.0
        x0 = xmin; x1 = 0.0
        aspectRatio = aspectRatio_list[iStep]
        y0 = 0.0 ; y1 = -1.0*aspectRatio*xmin


#        # Color Version        
#        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=12, xmin=-24.0, xmax=0.0, ymin=0.0,ymax=1.00,nLayersY=5,minStrain=2.0,maxStrain=10.0,mainDir='x')
#        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
#        if (iSim==1):
#            pyMax = np.max(PartY)
#            pyPad = 0.56
#            aspectRatio_list[iStep] = ((pyMax+pyPad)/(-xmin))
#            print('aspectRatio = %0.3f' % ((pyMax+pyPad)/(-xmin)))
#
#        CMAP=([[]])
#        CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,.0,.95,.0],
#                           RGBShift=[[0.0, 0.0, 0.0], 
#                                     [0.0, 0.0, 0.0], 
#                                     [-.2,-.2,0.2], 
#                                     [0.0, 0.0, 0.0]])
#        plt.register_cmap(cmap=CMAP)
#        plt.set_cmap("myColorMap")
        
        
        ## Plot Left: Vel + Disp
        if iStep==0:
            from computeDisplacement import computeDisplacement
            [X_center,Y_center,Mag,BackDisp] = computeDisplacement(dataFolder)
            
            plt.contourf(X_center,Y_center,Mag/BackDisp,np.linspace(.0,2.0,20))
            
            plt.axis([xmin,0.0,0.0,-1.0*aspectRatio*xmin])
            plt.axis('off')


#
#            raw = Output.getData(dataFolder + 'Vx.bin',True)
#            xmin = raw.xmin; ymin = raw.ymin
#            xmax = raw.xmax; ymax = raw.ymax
#            x = np.linspace(raw.xmin,raw.xmax,raw.nx)
#            y = np.linspace(raw.ymin,raw.ymax,raw.ny-1)
#            dx = x[1]-x[0]
#            X, Y = np.meshgrid(x,y)
#            X = X.T
#            Y = Y.T
#            Vx = Output.getData(dataFolder + 'Vx.bin',True).data
#            Vx = 0.5*(Vx[:,0:-1] + Vx[:,1:])
#            Vy = Output.getData(dataFolder + 'Vy.bin',True).data
#            Vy = 0.5*(Vy[0:-1,:] + Vy[1:,:])
#            phase = Output.getData(dataFolder + 'phase.bin',False).data
#            phase = 0.5*(phase[:,0:-1] + phase[:,1:])
#            phase = 0.5*(phase[0:-1,:] + phase[1:,:])
#            
#            Vx = Vx*phase
#            Vy = Vy*phase
#            
#            Vi = np.sqrt(Vx**2+Vy**2)
#            
#            vmin = -1.0
#            vmax = 1.0
#            Vc = Vx.copy()
#            Vc[Vx<vmin] = vmin
#            Vc[Vx>vmax] = vmax
#            
#    
#    
#            plt.streamplot(x,y,Vx.T,Vy.T,color="w", density=2.5,arrowsize=1)       


#        ## Plot Right: Stress orientation
#        elif iSim==2:
#            dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
#            Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
#            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
#            
#            PartX = []
#            PartY = []
#            PartPattern = []
#            
#            ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
#            
#            
#     
#    
#            ymax = 4.0
#            plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
#            plt.axis("off")
#            rx = 1
#            ry = 1
#            
#            phase = Output.getData(dataFolder + 'phase.bin',False).data[::rx,::ry]
#            mask = phase==0
#            
#            Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True).data[::rx,::ry]
#            Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True).data[::rx,::ry]
#            if iSim == 0:
#                raw = Output.getData(dataFolder + 'sigma_xx.bin',True)
#                xmin = raw.xmin; ymin = raw.ymin
#                xmax = raw.xmax; ymax = raw.ymax
#                x = np.linspace(raw.xmin,raw.xmax,raw.nx)
#                y = np.linspace(raw.ymin,raw.ymax,raw.ny)
#                dx = x[1]-x[0]
#                X, Y = np.meshgrid(x,y)
#                X = X.T
#                Y = Y.T
#            Tau = Sxx / Sxy;
#            SII = np.sqrt(Sxx**2+Sxy**2);
#    
#            psi = np.zeros(Sxy.shape)
#            I = Sxy<0.0
#    
#            psi[I] = np.arctan(-Tau[I]+np.sqrt(Tau[I]**2+1.0))
#            psi[~I] = np.arctan(-Tau[~I]-np.sqrt(Tau[~I]**2+1.0))
#            
#    #        psi = np.ma.masked_array(psi, mask_sub)
#    
#            length = 0.2 * phase[::rx,::ry].flatten()
#    
#            xLine = arr([-np.cos(psi.flatten())*length,np.cos(psi.flatten())*length])
#            yLine = arr([-np.sin(psi.flatten())*length,np.sin(psi.flatten())*length])
#            
#            Svec_x = np.cos(psi)
#            Svec_y = np.sin(psi)
#            
#            Svec_x *= phase
#            Svec_y *= phase
#            
#                 
#            
#            
#            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
#            if iSim == 0:
#                CMAP = arr([ [0.4,0.5,0.8,1.0] ])
#                CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP)
#                plt.register_cmap(cmap=CMAP)
#                plt.set_cmap("myColorMap")
#            plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')
#            
#            smin = 0.0*1e6
#            smax = 10.0*1e6
#            SII[SII>smax] = smax
#            SII[SII<smin] = smin
#            n = 100
#            start_points = arr([np.linspace(xmin+dx,xmax-dx,n),
#                                np.linspace(dx,dx,n)]).T
#            
#            colors = arr([[  0,255,255],
#                          [255,255,255],
#                          [255,255,  0]]) / 255.0
#            nTot = 256
#            CMAP = LinearSegmentedColormap.from_list('custom2',colors,N=nTot)       
#            plt.register_cmap(cmap=CMAP)
#            plt.streamplot(x[::rx],y[::ry],Svec_x.T,Svec_y.T,color=SII.T, density=2.0,arrowsize=0.01,cmap='custom2',linewidth=1.0)        
#    