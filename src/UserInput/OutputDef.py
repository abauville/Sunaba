#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 15:41:07 2017

@author: abauville
"""

from InputDef import Frozen
import numpy as np
import json
import os


def readJson(Filename):
    with open(Filename) as data_file:    
        state = json.load(data_file)
        return state


class dataSet(Frozen):
     _Frozen__List = ["nx","ny","xmin","xmax","ymin","ymax","charUnit","data"]
     def __init__(self,isDimensional=False):
         self.nx        = 0
         self.ny        = 0
         self.xmin      = 0
         self.xmax      = 0
         self.ymin      = 0
         self.ymax      = 0
         self.charUnit  = 0
         self.data      = 0

    
def getData(FileName):
    
    
    myDataSet = dataSet()
    
    #with open(FileName, mode='rb') as file: # b is important -> binary
    #    fileContent = file.read()
    
    
    # note: making a custom ndtype would be more elegant
    f = open(FileName, "rb")
    nx = np.fromfile(f, dtype=np.int32, count=1, sep='')
    
    f.seek(4, os.SEEK_SET)
    ny = np.fromfile(f, dtype=np.int32, count=1, sep='')

    f.seek(8, os.SEEK_SET)
    xmin = np.fromfile(f, dtype=np.double, count=1, sep='')

    f.seek(16, os.SEEK_SET)
    xmax = np.fromfile(f, dtype=np.double, count=1, sep='')
    
    f.seek(24, os.SEEK_SET)
    ymin = np.fromfile(f, dtype=np.double, count=1, sep='')

    f.seek(32, os.SEEK_SET)
    ymax = np.fromfile(f, dtype=np.double, count=1, sep='')

    f.seek(40, os.SEEK_SET)
    charUnit = np.fromfile(f, dtype=np.double, count=1, sep='')
    
    f.seek(48, os.SEEK_SET)    
    data = np.fromfile(f, dtype=np.double, count=-1, sep='')
    data = np.reshape(data, (ny,nx));
    f.close()
    
    
    myDataSet.nx = nx[0]
    myDataSet.ny = ny[0]
    myDataSet.xmin = xmin[0]
    myDataSet.xmax = xmax[0]
    myDataSet.ymin = ymin[0]    
    myDataSet.ymax = ymax[0]   
    myDataSet.charUnit = charUnit[0]
    myDataSet.data = data
    
    
    return myDataSet













