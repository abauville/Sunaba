#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:09:30 2018

@author: abauville
"""
#import sys
#sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import OutputDef as Output





def get_XYandPattern(dataFolder,lc=2.0e3,sampleRate=1, nLayersX = 0, nLayersY=7,minStrain=0.1,maxStrain=5.0,xmin = 'auto',xmax = 'auto',ymin='auto',ymax='auto', mainDir='x'):
    ## Load and subsample
    # ================================
    PartX  = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
    PartY  = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
    PartXIni  = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sampleRate]/lc
    PartYIni  = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sampleRate]/lc

    PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data[0::sampleRate]

    if xmin=='auto':
        xmin = np.min(PartX)
    if xmax=='auto':
        xmax = np.max(PartX)
    if ymin=='auto':
        ymin = np.min(PartY)
    if ymax=='auto':
        ymax = np.max(PartY)
        
    ## Geometrical pattern
    # ================================
#    nLayersY = 7
#    nLayersX = 16
#    xmax = 0.0
#    xmin = -32.0#np.min(PartXIni)

    if mainDir == 'x':
        if nLayersX > 0:
            YPattern = PartYIni
            YPattern -= ymin
            YPattern /= ymax-ymin
    #        YPattern = np.cos(PartYIni*1.0*np.pi*nLayersY+np.pi/2.0)
            YPattern = np.cos(YPattern*1.0*np.pi*nLayersY+np.pi/2.0)
            maxVal= np.max(YPattern)
            YPattern = np.round((YPattern+maxVal)/(2.0*maxVal))
            
            nColors = nLayersX
            
            XPattern = PartXIni
#            XPattern -= xmin
            XPattern /= xmax-xmin
        
            
            PartPattern= 4.0*np.floor(-XPattern * (nLayersX) )
    #        
        
#        
            PartPattern+= 2.0*YPattern
            PartPattern = PartPattern%(4*nLayersX)
            
        else:
            PartPattern = np.zeros(PartStrain.shape)
            nColors = 1
    elif mainDir=='y':
        if nLayersY>0:
            PartPattern = np.ones(PartStrain.shape)
            nColors = nLayersY
            
            if nLayersX<=1:
                nLayersX=0
            
            XPattern = PartXIni
            XPattern -= xmin
            XPattern /= xmax-xmin
            XPattern = np.cos(PartXIni*1.0*np.pi*nLayersX+np.pi/2.0)
            maxVal= np.max(XPattern)
            
            if nLayersX>1:
                XPattern = np.round((XPattern+maxVal)/(2.0*maxVal))
            else:
                XPattern = np.round((XPattern)/(2.0*maxVal))
#            XPattern = np.round((XPattern)/(2.0*maxVal))
#            

            
            YPattern = PartYIni
            YPattern /= ymax-ymin
#            YPattern = -YPattern
            
#            YPattern[YPattern>0.95] = 6.0
#            nColors = 6
            
            PartPattern= 4.0*np.floor(YPattern * (nLayersY) )  
            PartPattern+= 2.0*XPattern
            PartPattern = PartPattern%(4*nLayersY)

        else:
            PartPattern = np.zeros(PartStrain.shape)
            nColors = 1
    else:
        raise ValueError("mainDir can take the values 'x' or 'y' only.")
    
                
    ## Strain pattern
    # ================================
    PartStrain -= minStrain
    PartStrain[PartStrain<0.0] = 0.0
    PartStrain /= maxStrain-minStrain
    
    PartStrainLogical = PartStrain>1.0
    PartStrain[PartStrainLogical] = 1.0        
    
    PartPattern += PartStrain



    return (PartX, PartY, PartPattern, nColors)












def getColormap(nColors,name='KellyMod_Layers_Strain',renderer=0,CMAP=np.array([]), 
                darknessFactor=[1.0,0.0,0.5,0.0], 
                RGBShift = [[0.0,0.0,0.0], # Main color
                            [0.0,0.0,0.0], # Main Color strain
                            [0.0,0.0,0.0], # Minor Color
                            [0.0,0.0,0.0]], # Minor Color strain
                ):
    # Darkness Factor is applied before RGB shift
    
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
#        condensedCMAPtemp[:,0  ] += .4
#        condensedCMAPtemp[:,1  ] += .35
        
#        condensedCMAPtemp[:,0:3] *= 0
#        condensedCMAPtemp[:,0:3] += 1.0
#        condensedCMAPtemp = np.ones((nColors,4))
        condensedCMAPtemp[:,0  ] += .4
        condensedCMAPtemp[:,1  ] += .35
    
        
        
#            
            
    else:
        try:
            if (CMAP.shape[1]!=4):
                raise ValueError("error in getColormap: CMAP must have dimensions nColors*4");
        except IndexError: # Occurs if CMAP contains only one color CMAP.shape = (4,)
            if (CMAP.shape[0]!=4):
                raise ValueError("error in getColormap: CMAP must have dimensions nColors*4");
        nColors = CMAP.shape[0]
        condensedCMAPtemp = CMAP
    
    
    
    condensedCMAPtemp[condensedCMAPtemp>1.0] = 1.0
    condensedCMAPtemp[condensedCMAPtemp<0.0] = 0.0
    
    condensedCMAP = np.zeros((4*nColors,4))
    for i in range(4):
        if i%2==0:
            condensedCMAP[i::4,:] = condensedCMAPtemp*darknessFactor[i] #* 0.4 + .3
        else:
            condensedCMAP[i::4,:] = condensedCMAP[i-1::4,:]*darknessFactor[i] #* 0.4 + .3

        condensedCMAP[i::4,0] += RGBShift[i][0]  # Shift Red
        condensedCMAP[i::4,1] += RGBShift[i][1]  # Shift Green
        condensedCMAP[i::4,2] += RGBShift[i][2]  # Shift Blue
    

    
    condensedCMAP[condensedCMAP>1.0] = 1.0
    condensedCMAP[condensedCMAP<0.0] = 0.0
    condensedCMAP[:,-1] = 1.0
#    nColors = condensedCMAP.shape[0]
    
    
    
    if renderer == 0: # renderer is Matplotlib
        res = 12
        CMAP = LinearSegmentedColormap.from_list(name,condensedCMAP,N=res*4*nColors-(res-1))
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