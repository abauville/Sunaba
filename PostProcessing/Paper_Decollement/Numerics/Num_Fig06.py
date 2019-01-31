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
from matplotlib.colors import LinearSegmentedColormap



#chi_list = [ 1, 10, 20, 30, 40, 50, 60, 70, 80]
#chi_list = [ 1, 20, 40,70]
chi_list = [10,15,20,25]
Lambda_list = [60]
nC = len(chi_list)
nL = len(Lambda_list)
nHor = len(Lambda_list)
nVer = len(chi_list)


aspectRatio = 0.3
fig  = Figz_Utils.Figure(106,height=29.7,width=21.0,mode='draft')
#bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.66,leftMarginPad=1.25,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')
#ax = plt.gca()
##plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
#ax.xaxis.tick_top()
#ax.xaxis.set_label_position('top')
#ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
##plt.axis([.0,100.0,100.0,.0])

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
    sampleRate = 1
    pointSize = sampleRate/3.0
else:
    sampleRate = 50
    pointSize = sampleRate/12.0


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
        
        
#        if iSim == 0:
#            raw = Output.getData(dataFolder + 'Vx.bin',True)
#            xmin = raw.xmin; ymin = raw.ymin
#            xmax = raw.xmax; ymax = raw.ymax
#            x = np.linspace(raw.xmin,raw.xmax,raw.nx)
#            y = np.linspace(raw.ymin,raw.ymax,raw.ny-1)
#            dx = x[1]-x[0]
#            X, Y = np.meshgrid(x,y)
#            X = X.T
#            Y = Y.T
#        Vx = Output.getData(dataFolder + 'Vx.bin',True).data
#        Vx = 0.5*(Vx[:,0:-1] + Vx[:,1:])
#        Vy = Output.getData(dataFolder + 'Vy.bin',True).data
#        Vy = 0.5*(Vy[0:-1,:] + Vy[1:,:])
#        phase = Output.getData(dataFolder + 'phase.bin',False).data
#        phase = 0.5*(phase[:,0:-1] + phase[:,1:])
#        phase = 0.5*(phase[0:-1,:] + phase[1:,:])
#        
#        Vx = Vx*phase
#        Vy = Vy*phase
#        
#        Vi = np.sqrt(Vx**2+Vy**2)
#        
#        vmin = -1.0
#        vmax = 1.0
#        Vc = Vx.copy()
#        Vc[Vx<vmin] = vmin
#        Vc[Vx>vmax] = vmax
#        
#
#
#        plt.streamplot(x,y,Vx.T,Vy.T,color="w", density=2.5,arrowsize=1)        

        ymax = 4.0
        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
        plt.axis("off")
        rx = 1
        ry = 1
        
        phase = Output.getData(dataFolder + 'phase.bin',False).data[::rx,::ry]
        mask = phase==0
        
        Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True).data[::rx,::ry]
        Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True).data[::rx,::ry]
        if iSim == 0:
            raw = Output.getData(dataFolder + 'sigma_xx.bin',True)
            xmin = raw.xmin; ymin = raw.ymin
            xmax = raw.xmax; ymax = raw.ymax
            x = np.linspace(raw.xmin,raw.xmax,raw.nx)
            y = np.linspace(raw.ymin,raw.ymax,raw.ny)
            dx = x[1]-x[0]
            X, Y = np.meshgrid(x,y)
            X = X.T
            Y = Y.T
        Tau = Sxx / Sxy;
        SII = np.sqrt(Sxx**2+Sxy**2);

        psi = np.zeros(Sxy.shape)
        I = Sxy<0.0

        psi[I] = np.arctan(-Tau[I]+np.sqrt(Tau[I]**2+1.0))
        psi[~I] = np.arctan(-Tau[~I]-np.sqrt(Tau[~I]**2+1.0))
        
#        psi = np.ma.masked_array(psi, mask_sub)

        length = 0.2 * phase[::rx,::ry].flatten()

        xLine = arr([-np.cos(psi.flatten())*length,np.cos(psi.flatten())*length])
        yLine = arr([-np.sin(psi.flatten())*length,np.sin(psi.flatten())*length])
        
        Svec_x = np.cos(psi)
        Svec_y = np.sin(psi)
        
        Svec_x *= phase
        Svec_y *= phase
        
        
#        xCenter = X[::rx,::ry].flatten()
#        xCenter = arr([xCenter,xCenter])
#        yCenter = Y[::rx,::ry].flatten()
#        yCenter = arr([yCenter,yCenter])
#        
#        plt.plot( xCenter + xLine , yCenter + yLine ,'-w' )
        
#        
#        Xmask = X[::rx,::ry]
#        Xmask = np.ma.masked_array(Xmask, mask)
#        Ymask = Y[::rx,::ry]
#        Ymask = np.ma.masked_array(Ymask, mask)
#        SIImask = SII
        
#        plt.contourf(Xmask,Ymask,SIImask,256,cmap='seismic',vmin=0.0,vmax=10.0*1e6)
        
#        plt.colorbar()  
        
             
        
        
        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
        if iSim == 0:
            CMAP = arr([ [0.4,0.5,0.8,1.0] ])
            CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP)
            plt.register_cmap(cmap=CMAP)
            plt.set_cmap("myColorMap")
        plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')
        
        smin = 0.0*1e6
        smax = 10.0*1e6
        SII[SII>smax] = smax
        SII[SII<smin] = smin
#        SII[phase==0] = np.nan#0.5*(smin+smax) # works when using a colorbar centered on white
#        SIImask = np.ma.masked_array(SII, mask)
        n = 100
#        I = (np.random.rand(n)*PartX.size).astype(np.int)
#        start_points = arr([  PartX[I], PartY[I] ]).T
        start_points = arr([np.linspace(xmin+dx,xmax-dx,n),
                            np.linspace(dx,dx,n)]).T
        
        colors = arr([[  0,255,255],
                      [255,255,255],
                      [255,255,  0]]) / 255.0
        nTot = 256
        CMAP = LinearSegmentedColormap.from_list('custom2',colors,N=nTot)       
        plt.register_cmap(cmap=CMAP)
        plt.streamplot(x[::rx],y[::ry],Svec_x.T,Svec_y.T,color=SII.T, density=2.0,arrowsize=0.01,cmap='custom2',linewidth=1.0)        
#        plt.colorbar()
#        plt.plot(start_points[:,0],start_points[:,1],'or')
##         
        
        
        iSim+=1
