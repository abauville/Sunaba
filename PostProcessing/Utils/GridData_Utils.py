#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:14:13 2019

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, '../../PostProcessing/Paper_Decollement/CriticalTaper')
import OutputDef as Output
import matplotlib.pyplot as plt
import numpy as np
from numpy import array as arr


def stress_orientation(dataFolder,sampleRateX=1,sampleRateY=1,which='S1'):
    
    
#    phase = Output.getData(dataFolder + 'phase.bin',False).data[::sampleRateX,::sampleRateY]
#    mask = phase==0
    
    Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True).data[::sampleRateX,::sampleRateY]
    Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True).data[::sampleRateX,::sampleRateY]
    Tau = Sxx / Sxy;
    I = Sxy<0.0
    
        
    SII = np.sqrt(Sxx**2+Sxy**2);
    psi = np.zeros(Sxy.shape)

    psi[I] = np.arctan(-Tau[I]+np.sqrt(Tau[I]**2+1.0))
    psi[~I] = np.arctan(-Tau[~I]-np.sqrt(Tau[~I]**2+1.0))
    
    if which == 'S1':
        daijoubu = 1
    elif which == 'S3':
        psi+=np.pi/2.0
    else:
        raise ValueError("which must be 'S1' or 'S3'.")        
    
#        psi = np.ma.masked_array(psi, mask_sub)

#    length = 0.2 * phase[::sampleRateX,::sampleRateY].flatten()

#    xLine = arr([-np.cos(psi.flatten())*length,np.cos(psi.flatten())*length])
#    yLine = arr([-np.sin(psi.flatten())*length,np.sin(psi.flatten())*length])
    
    Svec_x = np.cos(psi)
    Svec_y = np.sin(psi)
    
#    Svec_x *= phase
#    Svec_y *= phase
        
    return [Svec_x,Svec_y,SII,psi]
    