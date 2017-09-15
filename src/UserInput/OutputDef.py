#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 15:41:07 2017

@author: abauville
"""

from InputDef import Frozen
import InputDef as Input
import MaterialsDef as Materials

import numpy as np
import json
import os
from collections import namedtuple


def readJson(Filename):
    with open(Filename) as data_file:    
        state = json.load(data_file)
        return state
    
def readState(Filename):
    with open(Filename) as data_file:    
        myFileData = json.load(data_file,object_hook=lambda d: namedtuple('state', d.keys())(*d.values()))
        return myFileData
    
def readInput(Filename):
    with open(Filename) as data_file:    
        myFile = json.load(data_file)

        Setup = Input.Setup()
        
        Setup.Char.setFromDict(myFile['Char'])
        
        
        Setup.Grid.setFromDict(myFile['Grid'])
        Setup.IC.setFromDict(myFile['IC'])
        
        Setup.Numerics.setFromDict(myFile['Numerics'])
        Setup.Output.setFromDict(myFile['Output'])
        Setup.Particles.setFromDict(myFile['Particles'])
        Setup.Physics.setFromDict(myFile['Physics'])
        Setup.Visu.setFromDict(myFile['Visu'])
        
        Setup.Description = myFile['Description']
        
        Setup.BC.Stokes.setFromDict(myFile['BC']['Stokes'])
        Setup.BC.Thermal.setFromDict(myFile['BC']['Thermal'])
        
        
        for key in myFile['MatProps']:
            Setup.MatProps[key] = Input.Material()
            Setup.MatProps[key].setFromDict(myFile['MatProps'][key])
            Setup.MatProps[key].vDisl = Materials.DislocationCreep   ("Off")
            Setup.MatProps[key].vDiff = Materials.DiffusionCreep     ("Off")
            Setup.MatProps[key].vPei  = Materials.PeierlsCreep       ("Off")
            Setup.MatProps[key].vDisl.setFromDict(myFile['MatProps'][key]['vDisl'])
            Setup.MatProps[key].vDiff.setFromDict(myFile['MatProps'][key]['vDiff'])
            Setup.MatProps[key].vPei .setFromDict(myFile['MatProps'][key]['vPei' ])
            
        #Note: Geometry is left empty for the moment (because it involves different kinds of classes)
            
        return Setup


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

    def interpFromCellsToCells(self):
        # assign new data matrix
        CellData = np.zeros( ( self.nx+1, self.ny+1 ) )
        
        # interpolate from nodes to inner cell centers
        tempData                = (self.data[:-1,:  ] + self.data[1: ,:  ]) / 2.0
        CellData[1:-1,1:-1]     = (tempData[:  ,:-1] + tempData[:  ,1: ]) / 2.0
        
        # copy inner data to sides
        CellData[ 0,: ] = CellData[ 1,: ]     # bottom
        CellData[-1,: ] = CellData[-2,: ]     # top
        CellData[ :, 0] = CellData[ :, 1]     # left
        CellData[ :,-1] = CellData[ :,-2]     # right    
        
        dx = (self.xmax - self.xmin)/self.nx
        dy = (self.ymax - self.ymin)/self.ny   
         
        # update dimensions
        self.nx += 1
        self.ny += 1
        self.xmin -= dx
        self.xmax += dx
        self.ymin -= dy
        self.ymax += dy
        
        self.data = CellData
        
        
        
def getData(FileName):
    myDataSet = dataSet()
    
    #with open(FileName, mode='rb') as file: # b is important -> binary
    #    fileContent = file.read()
    
    
    # note: making a custom ndtype would be more elegant
    f = open(FileName, "rb")
    nx = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]
    
    f.seek(4, os.SEEK_SET)
    ny = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]

    f.seek(8, os.SEEK_SET)
    xmin = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(16, os.SEEK_SET)
    xmax = np.fromfile(f, dtype=np.double, count=1, sep='')[0]
    
    f.seek(24, os.SEEK_SET)
    ymin = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(32, os.SEEK_SET)
    ymax = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(40, os.SEEK_SET)
    charUnit = np.fromfile(f, dtype=np.double, count=1, sep='')[0]
    
    f.seek(48, os.SEEK_SET)    
    data = np.fromfile(f, dtype=np.double, count=-1, sep='')
    f.close()
    data = np.reshape(data, (nx,ny),order='F')
#    data = np.reshape(data, (ny,nx))
#    data = np.transpose(data)
    
    myDataSet.nx = nx
    myDataSet.ny = ny
    myDataSet.xmin = xmin
    myDataSet.xmax = xmax
    myDataSet.ymin = ymin    
    myDataSet.ymax = ymax   
    myDataSet.charUnit = charUnit
    myDataSet.data = data
    
    return myDataSet

def getDataMatrix(FileName):    
    f = open(FileName, "rb")
    nx = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]
    
    f.seek(4, os.SEEK_SET)
    ny = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]

    f.seek(48, os.SEEK_SET)    
    data = np.fromfile(f, dtype=np.double, count=-1, sep='')
    f.close()
    data = np.reshape(data, (nx,ny),order='F')

    return data    
    
def getDataInfo(FileName):
    myDataSet = dataSet()
    
    # note: making a custom ndtype would be more elegant
    f = open(FileName, "rb")
    nx = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]
    
    f.seek(4, os.SEEK_SET)
    ny = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]

    f.seek(8, os.SEEK_SET)
    xmin = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(16, os.SEEK_SET)
    xmax = np.fromfile(f, dtype=np.double, count=1, sep='')[0]
    
    f.seek(24, os.SEEK_SET)
    ymin = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(32, os.SEEK_SET)
    ymax = np.fromfile(f, dtype=np.double, count=1, sep='')[0]

    f.seek(40, os.SEEK_SET)
    charUnit = np.fromfile(f, dtype=np.double, count=1, sep='')[0]
    
    myDataSet.nx = nx
    myDataSet.ny = ny
    myDataSet.xmin = xmin
    myDataSet.xmax = xmax
    myDataSet.ymin = ymin    
    myDataSet.ymax = ymax   
    myDataSet.charUnit = charUnit
    
    return myDataSet

def getData_OneCell(FileName,ix,iy):
    f = open(FileName, "rb")
    nx = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]
    
    f.seek(4, os.SEEK_SET)
    ny = np.fromfile(f, dtype=np.int32, count=1, sep='')[0]

    if (ix>nx-1):
        raise ValueError( "ix=%i is larger than nx=%i" % (ix,nx) )
    if (iy>ny-1):
        raise ValueError( "iy=%i is larger than ny=%i" % (iy,ny) )


    Position = 48 + (ix+iy*nx)*8
    f.seek(Position, os.SEEK_SET)    
    data = np.fromfile(f, dtype=np.double, count=1, sep='')[0]
    f.close()

    return data



def getNumberOfOutFolders(rootFolder):
    return len(os.listdir(rootFolder)) - 2
    









