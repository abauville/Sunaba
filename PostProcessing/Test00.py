# Output reading test

import sys
import os
sys.path.insert(0, '../src/UserInput')
#from InputDef import *

import struct

import matplotlib.pyplot as plt
import numpy as np



import json
from pprint import pprint

modelStateFile = '../../Output_StokesFD/Test00/timeStep_00001/modelState.json'

with open(modelStateFile) as data_file:    
    data = json.load(data_file)

#pprint(data)

dataFile = '../../Output_StokesFD/Test00/timeStep_00200/Viscosity.bin'
with open(dataFile, mode='rb') as file: # b is important -> binary
    fileContent = file.read()
    

# note: making a custom ndtype would be more elegant
f = open(dataFile, "rb")
nx = np.fromfile(f, dtype=np.int32, count=1, sep='')
nx = nx[0]
f.seek(4, os.SEEK_SET)
ny = np.fromfile(f, dtype=np.int32, count=1, sep='')
ny = ny[0]
f.seek(8, os.SEEK_SET)

xmin = np.fromfile(f, dtype=np.double, count=1, sep='')
xmin = xmin[0]
f.seek(16, os.SEEK_SET)
xmax = np.fromfile(f, dtype=np.double, count=1, sep='')
xmax = xmax[0]
f.seek(24, os.SEEK_SET)
ymin = np.fromfile(f, dtype=np.double, count=1, sep='')
ymin = ymin[0]
f.seek(32, os.SEEK_SET)
ymax = np.fromfile(f, dtype=np.double, count=1, sep='')
ymax = ymax[0]
f.seek(40, os.SEEK_SET)


Char_quantity = np.fromfile(f, dtype=np.double, count=1, sep='')
Char_quantity = Char_quantity[0]
f.seek(48, os.SEEK_SET)

data = np.fromfile(f, dtype=np.double, count=-1, sep='')
data = np.reshape(data, (ny,nx));
f.close()

x = np.linspace(xmin,xmax,nx)


y = np.linspace(ymin,ymax,ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)
data = np.transpose(data)

plt.pcolor(xv,yv,np.log10(data))
plt.axis('equal')
plt.show()




