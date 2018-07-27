#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:09:30 2018

@author: abauville
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import OutputDef as Output





def get_XYandPattern(dataFolder,lc=2.0e3,sampleRate=1, nLayersX = 16, nLayersY=7):
    ## Load and subsample
    # ================================
    PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
    PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
    PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sampleRate]/lc
    PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sampleRate]/lc

    PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sampleRate]

    
    ## Geometrical pattern
    # ================================
#    nLayersY = 7
#    nLayersX = 16
    xmax = 0.0
    xmin = -32.0#np.min(PartXIni)
    YPattern = np.cos(PartYIni*1.0*np.pi*nLayersY+np.pi/2.0)
    maxVal= np.max(YPattern)
    YPattern = np.round((YPattern+maxVal)/(2.0*maxVal))
    
    nColors = nLayersX
    XPattern = PartXIni
    XPattern /= xmax-xmin
    XPattern = -XPattern

    PartPattern= 4.0*np.floor(XPattern * (nLayersX) )
    PartPattern+= 2.0*YPattern
    
                
    ## Strain pattern
    # ================================
    maxStrain = 5.0
    minStrain = .1
    PartStrain -= minStrain
    PartStrain[PartStrain<0.0] = 0.0
    PartStrain /= maxStrain-minStrain
    
    PartStrainLogical = PartStrain>1.0
    PartStrain[PartStrainLogical] = 1.0        
    PartPattern += 1.0*PartStrain



    return (PartX, PartY, PartPattern, nColors)












def getColormap(nColors,name='KellyMod_Layers_Strain',renderer=0,CMAP=np.array([]), shiftHLayerColors = True):
    ## Kenneth Kelly's 22 colour of maximum contrast
#    nColors  = 22
    
    if type(CMAP) is list:
        CMAP=np.array(CMAP)
        
    if (CMAP.size==0):
        if nColors>22:
            raise ValueError("error: nColors should not exceed 22")
        condensedCMAPtemp = np.zeros((nColors,4))
        for i in range (nColors):
            if i ==  0: condensedCMAPtemp[ i,:] = [242, 243, 244, 255]
            if i ==  1: condensedCMAPtemp[ i,:] = [ 34,  34,  34, 255]
            if i ==  2: condensedCMAPtemp[ i,:] = [243, 195,   0, 255]
            if i ==  3: condensedCMAPtemp[ i,:] = [135,  86, 146,255]
            if i ==  4: condensedCMAPtemp[ i,:] = [243, 132,   0,255]
            if i ==  5: condensedCMAPtemp[ i,:] = [161, 202, 241,255]
            if i ==  6: condensedCMAPtemp[ i,:] = [190,   0,  50,255]
    #        if i ==  7: condensedCMAPtemp[ i,:] = [194, 178, 128,255]
            if i ==  7: condensedCMAPtemp[ i,:] = [ 50, 178, 128,255]
    #        if i ==  8: condensedCMAPtemp[ i,:] = [132, 132, 130,255]
            if i ==  8: condensedCMAPtemp[ i,:] = [162,  80, 200,255]
            if i ==  9: condensedCMAPtemp[ i,:] = [  0, 136,  86,255]
            if i == 10: condensedCMAPtemp[ i,:] = [230, 143, 172,255]
            if i == 11: condensedCMAPtemp[ i,:] = [  0, 103, 165,255]
            if i == 12: condensedCMAPtemp[ i,:] = [249, 147, 121,255]
            if i == 13: condensedCMAPtemp[ i,:] = [ 96,  78, 151,255]
            if i == 14: condensedCMAPtemp[ i,:] = [246, 166,   0,255]
            if i == 15: condensedCMAPtemp[ i,:] = [179,  68, 108,255]
            if i == 16: condensedCMAPtemp[ i,:] = [220, 211,   0,255]
            if i == 17: condensedCMAPtemp[ i,:] = [136,  45,  23,255]
            if i == 18: condensedCMAPtemp[ i,:] = [141, 182,   0,255]
            if i == 19: condensedCMAPtemp[ i,:] = [101,  69,  34,255]
            if i == 20: condensedCMAPtemp[ i,:] = [226,  88,  34,255] 
            if i == 21: condensedCMAPtemp[ i,:] = [ 43,  61,  38,255] 
            
        condensedCMAPtemp /= 255.0
        #CMAP = LinearSegmentedColormap.from_list('Kelly',condensedCMAPtemp)
        #plt.register_cmap(cmap=CMAP)
        
        condensedCMAPtemp[:,0:2] *= .5
        condensedCMAPtemp[:,0:2] += .1
        condensedCMAPtemp[:,0  ] += .4
        condensedCMAPtemp[:,1  ] += .35
    
        condensedCMAP = np.zeros((4*nColors,4))
        condensedCMAP[0::4,:] = condensedCMAPtemp #* 0.4 + .3
        condensedCMAP[1::4,:] = condensedCMAPtemp*.25
            
            
    else:
        try:
            if (CMAP.shape[1]!=4):
                raise ValueError("error in getColormap: CMAP must have dimensions nColors*4");
        except IndexError: # Occurs if CMAP contains only one color CMAP.shape = (4,)
            if (CMAP.shape[0]!=4):
                raise ValueError("error in getColormap: CMAP must have dimensions nColors*4");
        nColors = CMAP.shape[0]
        condensedCMAPtemp = CMAP
        
        condensedCMAP = np.zeros((4*nColors,4))
        condensedCMAP[0::4,:] = condensedCMAPtemp #* 0.4 + .3
#        condensedCMAP[1::4,:] = condensedCMAPtemp*.25
        
        
    
        
    if shiftHLayerColors:
        condensedCMAPtemp[:,2] += .1        # Modify the color of the horizontal layers
        condensedCMAPtemp[:,0:1] -= .1      # Modify the color of the horizontal layers
        condensedCMAP[2::4,:] = condensedCMAPtemp*.85   # Darken the color of the horizontal layers
    else: 
        condensedCMAP[2::4,:] = condensedCMAPtemp   # Darken the color of the horizontal layers
    
    condensedCMAP[3::4,:] = condensedCMAPtemp*.1   # Darken the color for strain
    
    condensedCMAP[condensedCMAP>1.0] = 1.0
    condensedCMAP[condensedCMAP<0.0] = 0.0
    condensedCMAP[:,-1] = 1.0
#    nColors = condensedCMAP.shape[0]
    
    if renderer == 0: # renderer is Matplotlib
        CMAP = LinearSegmentedColormap.from_list(name,condensedCMAP)
        return CMAP
        #plt.register_cmap(cmap=CMAP)
    
    if renderer == 1: # renderer is Mayavi
        CMAP = np.zeros((255,4))
        
        n = condensedCMAP.shape[0]
        for i in range(n-1):
            iS = int(round(255 *  i/(n-1) ))
            iE = int(round(255 * (i+1)/(n-1)))
            for j in range(4):
                CMAP[iS:iE,j] = np.linspace(condensedCMAP[i,j],condensedCMAP[i+1,j],int(iE-iS))
        
        CMAP = (CMAP*255.0).astype("uint8")
        return CMAP