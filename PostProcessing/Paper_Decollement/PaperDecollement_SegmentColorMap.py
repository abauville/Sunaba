#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:09:30 2018

@author: abauville
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
def getColormap(nColors,name='KellyMod_Layers_Strain'):
    ## Kenneth Kelly's 22 colour of maximum contrast
#    nColors  = 22
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
    
    condensedCMAPtemp[:,2] += .1        # Modify the color of the horizontal layers
    condensedCMAPtemp[:,0:1] -= .1      # Modify the color of the horizontal layers
    condensedCMAP[2::4,:] = condensedCMAPtemp*.85   # Darken the color of the horizontal layers
    
    condensedCMAP[3::4,:] = condensedCMAPtemp*.25   # Darken the color for strain
    
    condensedCMAP[condensedCMAP>1.0] = 1.0
    condensedCMAP[condensedCMAP<0.0] = 0.0
    condensedCMAP[:,-1] = 1.0
#    nColors = condensedCMAP.shape[0]
    
    CMAP = LinearSegmentedColormap.from_list(name,condensedCMAP)
    return CMAP
    #plt.register_cmap(cmap=CMAP)