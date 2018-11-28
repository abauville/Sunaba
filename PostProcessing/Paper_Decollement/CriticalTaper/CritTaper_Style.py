#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:27:07 2018

@author: abauville
"""
#import numpy as np
import numpy as np
from numpy import array as arr
from matplotlib.colors import LinearSegmentedColormap

class Style:
    def __init__(self):
        self.Setting = "Paper"
    #    Setting = "Presentation"
        if self.Setting == "Paper":
            self.fontdict = {'family': 'Montserrat',
                        'weight': 'bold',
                        'size': int(11)
                        }
        
        elif self.Setting == "Presentation":
            self.fontdict = {'family': 'Montserrat',
                        'weight': 'bold',
                        'size': int(18)
                        }
        
        self.colormap = "seismic"
        
        
        self.colorRef = arr([.25,.25,.75])
        self.colorBW  = arr([.25,.5,.5])
        self.colorFW  = arr([.75,.25,.25])
    
        self.alphaBW = 0.15
    
        self.colorRef_a = np.concatenate([self.colorRef,[self.alphaBW]])
        self.colorBW_a  = np.concatenate([self.colorBW ,[self.alphaBW]])
        self.colorFW_a  = np.concatenate([self.colorFW ,[self.alphaBW]])
    
    
        
    def getCmap_Type(self):
#        colors = arr([[200, 30, 32],
#                      [248,160, 32],
#                      [248,200,  0],
#                      [  0,200,248],
#                      [ 32,160,248],
#                      [ 32, 30,200]]) / 255.0
        colors = arr([[ 32, 30,200],
                      [ 32,160,248],
                      [  0,200,248],
                      [248,200,  0],
                      [248,160, 32],
                      [200, 30, 32]]) / 255.0
        colors = np.flipud(colors)
        nSeg = colors.shape[0]-1
        
        # algorithm to blend colors
        def blendColorValue(a, b, t):
            return np.sqrt((1 - t) * a**2 + t * b**2)
        
        C = 0
        nSubSegs = np.array([40,1,40,1,40])
        nTot = np.int(np.sum(nSubSegs)) + 1
        blendedColors = np.zeros([nTot,3])
        i = 0
        for iSeg in range(nSeg):
            iSub = 0
            for t in np.linspace(0.0,1.0,nSubSegs[iSeg]+1):
                if (iSeg<nSeg-1 and iSub==nSubSegs[iSeg]):
                    break
                else:
                    blendedColors[i,:] = blendColorValue(colors[iSeg,:],colors[iSeg+1,:],t)
                    i+=1
                    iSub+=1
                
        CMAP = LinearSegmentedColormap.from_list('custom',blendedColors,N=nTot)        
        return (CMAP,blendedColors)

        
        