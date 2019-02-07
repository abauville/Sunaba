#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 15:36:53 2019

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
def computeDisplacement(dataFolder):

    Char = Output.readInput(dataFolder + '../../Input/input.json').Char
    timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
    
    PartX = []
    PartY = []
    PartPattern = []
    
#    ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
    
  
    lc = 2e3
    sampleRate = 10
    rx = 1
    ry = 1
    PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
    PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
    PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sampleRate]/lc
    PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sampleRate]/lc

    PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sampleRate]
    
    
    # Compute the displacement
    
#    ymax = 5.0
#    plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
#    plt.axis("off")
    
    
    
    
    dum = Output.getData(dataFolder + 'phase.bin',False)
    xmin = dum.xmin
    xmax = dum.xmax
    ymin = dum.ymin
    ymax = dum.ymax
    nx = dum.nx
    ny = dum.ny
    
    
    XX, YY = np.meshgrid(np.linspace(xmin,xmax,nx),
                         np.linspace(ymin,ymax,ny))
    XX = XX.T
    YY = YY.T
    XX = XX[::rx,::ry]   
    YY = YY[::rx,::ry]   
    
    dx = XX[1,0]-XX[0,0]
    dy = YY[0,1]-YY[0,0]
    
    strain = Output.getData(dataFolder + 'strain.bin',True).data[::rx,::ry]
    strainCriteria = 1.0
    strainFlat = strain.copy().flatten()
    XFlat = XX.copy().flatten()
    YFlat = YY.copy().flatten()
    I = np.argwhere(strainFlat>strainCriteria)
    notI = np.argwhere(strainFlat<=strainCriteria)
    
    IJ = np.argwhere(strain>strainCriteria)
    notIJ = np.argwhere(strain<=strainCriteria)
    
    notIPart = np.argwhere(PartStrain<=strainCriteria)
    
    print("A")
    tic = time.time()
    XIni = np.zeros(XX.shape)
    Counter = np.zeros(XX.shape)

  
    I = np.floor((PartX-xmin)/dx).astype(np.int)
    J = np.floor((PartY-ymin)/dy).astype(np.int)
    ONES = np.ones(XX.shape)
    for i in range(PartXIni.shape[0]):
        XIni[I[i],J[i]] += PartXIni[i]
        Counter[I[i],J[i]]+=1.0
        
    XIni/=Counter
    
    
    XIniFlat = XIni.copy().flatten()
    
    print("A: %.2f s" % (time.time()-tic))
    
    
    print("B")
    tic = time.time()
    for i in range(IJ.shape[0]):
        I = IJ[i,0]
        J = IJ[i,1]
        x = XX[I,J]
        y = YY[I,J]
        
        # version Particles
#            Iclosest = np.argmin((x-PartX[notIPart])**2+(y-PartY[notIPart])**2)
#            XIni[I,J] = PartXIni[notIPart[Iclosest]]
        
#            # version Grid
#            Iclosest = np.argmin((x-XFlat[notI])**2+(y-YFlat[notI])**2)
#            XIni[I,J] = XIniFlat[notI[Iclosest]]
        
        # version Grid2
        w = 15
        Ib = I-w; 
        if Ib<0: 
            Ib = 0
        Ie = I+w; 
        if Ie<0: 
            Ie = 0
        Jb = J-w; 
        if Jb<0: 
            Jb = 0
        Je = J+w; 
        if Je<0: 
            Je = 0
        
        X = XX[Ib:Ie,Jb:Je]
        Y = YY[Ib:Ie,Jb:Je]
        XIni_small = XIni[Ib:Ie,Jb:Je]
        strain_small = strain[Ib:Ie,Jb:Je]
        notI = np.argwhere(strain_small.flatten()<=strainCriteria)
        if (notI.shape[0]>0):
            Iclosest = np.argmin((x-X.flatten()[notI])**2+(y-Y.flatten()[notI])**2)
            XIni[I,J] = XIni_small.flatten()[notI[Iclosest]]
   
    print("B: %.2f s" % (time.time()-tic))
    
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
    
    dum = .5*(XX[:,1:]+XX[:,:-1])
    X_center= .5*(dum[1:,:]+dum[:-1,:])
    
    dum = .5*(YY[:,1:]+YY[:,:-1])
    Y_center= .5*(dum[1:,:]+dum[:-1,:])
    print("D")
    
    ## Because Mag is localized on one cell it is difficult to see, so attribute the closest Mag to cells around
    MagFlat = Mag.copy().flatten()
    I = np.logical_and(MagFlat>0.1,strain_center.flatten()>strainCriteria)
    X_select = X_center.flatten()[I]
    Y_select = Y_center.flatten()[I]
    Mag_select = MagFlat.copy()[I]
    
    
    # 
    I = np.zeros(Mag.shape)
    w = 1 # thicknest of the window
    for i in range(Mag.shape[0]):
        ib = i-w
        ie = i+w
        if ib<0:
            ib=0
        if ie>=Mag.shape[0]:
            ie = Mag.shape[0]
        for j in range(Mag.shape[1]):
            jb = j-w
            je = j+w
            if jb<0:
                jb=0
            if je>=Mag.shape[1]:
                je = Mag.shape[1]
        
            if np.max(Mag[ib:ie,jb:je])>0.5:
                I[i,j] = True
            else:
                I[i,j] = False
                
                
    IJ = np.argwhere(I)
#        IJ = np.argwhere(strain_center>strainCriteria)
    for i in range(IJ.shape[0]):
        I = IJ[i,0]
        J = IJ[i,1]
#            x = XX[I,J]
#            y = YY[I,J]
        x = X_center[I,J]
        y = Y_center[I,J]
    
        Iclosest = np.argmin((x-X_select)**2+(y-Y_select)**2)
        Mag[I,J] = Mag_select[Iclosest]
    
    
    
#        IJ = np.argwhere(strain_center>strainCriteria)
    
#        Mag[strain_center<strainCriteria] = 0.0
    
    
#        plt.pcolor(XX,YY,XIni)
#        plt.pcolor(XX,YY,dXIni_dy)
#        plt.pcolor(XX,YY,Mag)
    Vpush = 10.0 * 0.01 / (365.25*24.0*3600.0)
    H = 2e3
    BackDisp = Vpush*timeSim/H
        
#    plt.contourf(X_center,Y_center,Mag/BackDisp,np.linspace(.0,2.0,20))
##        plt.contourf(X_center,Y_center,strain_center,np.linspace(.0,2.0,20))
##        I = np.argwhere(Mag.flatten()>1)
##        plt.scatter(X_center.flatten()[I],Y_center.flatten()[I],c=Mag.flatten()[I],s=15.0)
#    plt.colorbar()
#    plt.set_cmap('inferno')
    
    
    return [X_center,Y_center,Mag,BackDisp]
        
#            I_closest = np.unravel_index(np.argmin(x-XX, axis=None), a.shape)
        
        # find the closest grid point
        
    
        
  