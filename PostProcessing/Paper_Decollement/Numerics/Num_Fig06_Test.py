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
import time


#chi_list = [ 1, 10, 20, 30, 40, 50, 60, 70, 80]
#chi_list = [ 1, 20, 40,70]
chi_list = [20]
Lambda_list = [60]
nC = len(chi_list)
nL = len(Lambda_list)
nHor = len(Lambda_list)
nVer = len(chi_list)


aspectRatio = 0.3
fig  = Figz_Utils.Figure(108,height=29.7,width=21.0,mode='draft')
#bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.66,leftMarginPad=1.25,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')
#ax = plt.gca()
##plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
#ax.xaxis.tick_top()
#ax.xaxis.set_label_position('top')
#ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
##plt.axis([.0,100.0,100.0,.0])

Axes = Figz_Utils.makeAxes(fig,nVer,nHor,aspectRatio=aspectRatio,leftMarginPad=1.5,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')


superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select2/Beta00/"
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
#        outFolder = os.listdir(superRootFolder + superDirList[iSim] + "/Output/")[-1]
        outFolder = 'Out_01989'
#        outFolder = 'Out_01492'
#        outFolder = 'Out_01243'
        dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
        Char = Output.readInput(superRootFolder + superDirList[iSim] + "/Output/" +  'Input/input.json').Char
        timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
        
        PartX = []
        PartY = []
        PartPattern = []
        
        ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
        
      
#        ymax = 4.0
##        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
##        plt.axis("off")
#        
#        rx = 1
#        ry = 1
#        
#        
#        dum = Output.getData(dataFolder + 'phase.bin',False)
#        xmin = dum.xmin
#        xmax = dum.xmax
#        ymin = dum.ymin
#        ymax = dum.ymax
#        nx = dum.nx
#        ny = dum.ny
#        
#        dx = (xmax-xmin)/(nx-1)
#        dy = (ymax-ymin)/(ny-1)
#        
#        XX, YY = np.meshgrid(np.linspace(xmin,xmax,nx),
#                             np.linspace(ymin,ymax,ny))
#        XX = XX.T
#        YY = YY.T
#        XX = XX[::rx,::ry]   
#        YY = YY[::rx,::ry]   
#        
#        
#        phase = Output.getData(dataFolder + 'phase.bin',False).data[::rx,::ry]        
#        strain = Output.getData(dataFolder + 'strain.bin',True).data[::rx,::ry]
#        strain[phase==0] = np.nan
#        vmax = 150.0
#        plt.subplot(211)
#        plt.contourf(XX.T,YY.T,strain.T,np.linspace(.0,vmax,20),vmax=vmax)
#        Ix = 600
#        plt.plot(XX[Ix,:],YY[Ix,:],'r')
#        plt.axis('equal')
#        
##        plt.pcolor(XX,YY,strain.T,vmax=vmax)
#        plt.colorbar()
#        plt.set_cmap('inferno')
#
#        plt.subplot(212)
#        plt.plot(YY[Ix,:],np.cumsum(strain[Ix,:]*dx))
        
        
        
#        PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
        lc = 2e3
        sampleRate = 100
        PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
        PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
        PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sampleRate]/lc
        PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sampleRate]/lc
    
        PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sampleRate]
        
        
        # Compute the displacement
        ## First try
#        strainCriteria = 5.0
#        I = np.argwhere(PartStrain>strainCriteria)
#        
#        PartDisp = np.zeros(PartStrain.shape)
#        
#        for i in I:
#            x = PartX[i]
#            y = PartX[i]
#            plus_empty  = True
#            minus_empty = True
#            
#            I_xini_plus  = np.argwhere(PartXIni>PartXIni[i])
#            if I_xini_plus.shape[0]>0: # if not empty
#                Iclosest_plus  = np.argmin((x-PartX[I_xini_plus ])**2+(y-PartY[I_xini_plus ])**2)
#                plus_empty = False
#            
#            I_xini_minus = np.argwhere(PartXIni<PartXIni[i])
#            if I_xini_minus.shape[0]>0: # if not empty
#                Iclosest_minus = np.argmin((x-PartX[I_xini_minus])**2+(y-PartY[I_xini_minus])**2)
#                minus_empty = False
#                
#            if not(plus_empty) and not(minus_empty):
#                PartDisp[i] = PartXIni[I_xini_plus[Iclosest_plus]] - PartXIni[I_xini_minus[Iclosest_minus]]
#        
#        pointSize = 4.0*pointSize
##        plt.scatter(PartX,PartY,c=PartDisp,s=pointSize,vmin=0.0,edgecolors='None')
#        plt.scatter(PartX,PartY,c=PartXIni,s=pointSize,edgecolors='None')
#        plt.colorbar()
        
        
        ## Second try
        
        ymax = 4.0
        plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
        plt.axis("off")
        
        rx = 5
        ry = 5
        
        
        dum = Output.getData(dataFolder + 'phase.bin',False)
        xmin = dum.xmin
        xmax = dum.xmax
        ymin = dum.ymin
        ymax = dum.ymax
        nx = dum.nx
        ny = dum.ny
        
#        dx = (xmax-xmin)/(nx-1)
#        dy = (ymax-ymin)/(ny-1)
        
        XX, YY = np.meshgrid(np.linspace(xmin,xmax,nx),
                             np.linspace(ymin,ymax,ny))
        XX = XX.T
        YY = YY.T
        XX = XX[::rx,::ry]   
        YY = YY[::rx,::ry]   
        
        dx = XX[1,0]-XX[0,0]
        dy = YY[0,1]-YY[0,0]
        
        strain = Output.getData(dataFolder + 'strain.bin',True).data[::rx,::ry]
        strainCriteria = 5.0
        strainFlat = strain.copy().flatten()
        XFlat = XX.copy().flatten()
        YFlat = YY.copy().flatten()
        I = np.argwhere(strainFlat>strainCriteria)
        
        IJ = np.argwhere(strain>strainCriteria)
        
        notIPart = np.argwhere(PartStrain<=strainCriteria)
        print("A")
        tic = time.time()
        XIni = np.zeros(XX.shape)
        Counter = np.zeros(XX.shape)
        for i in range(PartXIni.shape[0]):
            I = np.int(np.floor((PartX[i]-xmin)/dx))
            J = np.int(np.floor((PartY[i]-ymin)/dy))
            XIni[I,J] += PartXIni[i]
            Counter[I,J]+=1.0
            
#        I = np.floor((PartX-xmin)/dx).astype(np.int)
#        J = np.floor((PartY-ymin)/dy).astype(np.int)
##        for i in range(PartXIni.shape[0]):
#        XIni[I,J] += PartXIni
#        Counter[I,J]+=1.0
#            
        XIni/=Counter
        
        
        # Attribute to points in the shear zone, the value of Xini from the closest point outside the shear zone
#        XIniFlat = np.zeros(XFlat.shape)
#        for i in I:
#            x = XFlat[i]
#            y = YFlat[i]
#            
#            Iclosest = np.argmin((x-PartX[notIPart])**2+(y-PartY[notIPart])**2)
#            XIniFlat[i] = PartXIni[notIPart[Iclosest]]
#            
#        XIni = XIniFlat.reshape(XX.shape)
        print("A: %.2f s" % (time.time()-tic))
        
        
        print("B")
        for i in range(IJ.shape[0]):
            I = IJ[i,0]
            J = IJ[i,1]
            x = XX[I,J]
            y = YY[I,J]
            
            Iclosest = np.argmin((x-PartX[notIPart])**2+(y-PartY[notIPart])**2)
            XIni[I,J] = PartXIni[notIPart[Iclosest]]
            
        
        print("C")
        
        ## Compute the gradient 
        dXIni_dx = XIni[1:,:]-XIni[:-1,:] # On Vx nodes
        dXIni_dy = XIni[:,1:]-XIni[:,:-1] # On Vy nodes
        
        # Interpolate values to cell corners
        dXIni_dx = .5*(dXIni_dx[:,1:]+dXIni_dx[:,:-1])
        dXIni_dy = .5*(dXIni_dy[1:,:]+dXIni_dy[:-1,:])
        
        ## Gradient magnitude
        Mag = np.sqrt(dXIni_dx**2+dXIni_dy**2)
        
        
        dum = .5*(strain[:,1:]+strain[:,:-1])
        strain_center= .5*(dum[1:,:]+dum[:-1,:])
        print("D")
#        IJ = np.argwhere(strain_center>strainCriteria)
        
#        Mag[strain_center<strainCriteria] = 0.0
        
        
#        plt.pcolor(XX,YY,XIni)
#        plt.pcolor(XX,YY,dXIni_dy)
        plt.pcolor(XX,YY,Mag)
        plt.colorbar()
        
            
            
#            I_closest = np.unravel_index(np.argmin(x-XX, axis=None), a.shape)
            
            # find the closest grid point
            
        
        
  
        
        iSim+=1



































