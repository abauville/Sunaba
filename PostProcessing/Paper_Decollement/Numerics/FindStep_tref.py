#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:47 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '..//CriticalTaper')
sys.path.insert(0, '../')
import OutputDef as Output
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern
from CritTaper_utils import Taper

# Misc
# =========================================
Style = CritTaper_Style.Style()

#   Define chi_list
# =========================================
#thisFile_chi_list = [1, 20, 60, 10, 40, 80]
chi_list = arr([1, 10, 20, 30, 40, 50, 60, 70, 80])
Lambda_list = arr([40,60,80])


#chi_list = arr([1, 10, 20, 40, 60, 80])
#Lambda_list = arr([60])



nC = len(chi_list)
nSim = nC





## Figure xFault
# ============================================

superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select/Beta00/"
Setup = Output.readInput(superRootFolder + "Weak01/Lambda60"  +  '/Output/Input/input.json')

yr = 365.25*24.0*3600.0
kyr = 1000.0 * yr
beta = 0.0
Lambda = 0.6
DataFolder = "/Users/abauville/Output/Paper_Decollement/Figz/Data/"


#loadedData_a = np.load(DataFolder + "locSlopes_Beta%02d_Lambda%02d_WeakDense_winSize64.npy" % (beta*10.0*180.0/np.pi,Lambda*100.0)).item(0)    
loadedData_b = np.load(DataFolder + "locSlopes_Beta%02d_all.npy" % (beta*10.0*180.0/np.pi)).item(0)



plot=0
folder = []
allSteps = []

for chi in chi_list/100.0:
    for Lambda in Lambda_list/100.0:

        thisData = loadedData_b["Lambda%02d_chi%02d" % (Lambda*100,chi*100)]
        locSlopes = thisData["locSlopes"]
        lenPrisms = thisData["lenPrisms"]
        tSteps = thisData["tSteps"]
        xFront = thisData["xFront"]
        
    
        list_front = [900]
        steps = np.zeros(len(list_front))
        for i in range(len(list_front)):
            Istep = np.argmin(np.abs(Setup.Grid.nxC-xFront - list_front[i]))
            steps[i] = tSteps[Istep]
        print("chi=%02d, Lambda=%02d, step=%i" % (chi*100, Lambda*100, steps[0]))
        
        folder.append("wWater/Beta00/Weak%02d/Lambda%02d/Output/Out_%05d" % (chi*100,Lambda*100,steps[0]))
        
#        subSteps = np.linspace(0,steps[0],9,dtype=np.int)
#        for subStep in subSteps[1:-1]:   
#            folder.append("wWater/Beta00/Weak%02d/Lambda%02d/Output/Out_%05d" % (chi*100,Lambda*100,subStep))
            
#        print(list(subSteps[1::]))
        allSteps.append(int(steps[0]))
    
    

    
# end iSim
print(list(allSteps))
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    