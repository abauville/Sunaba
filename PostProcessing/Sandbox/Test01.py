# Output reading test

import sys
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
#from pprint import pprint

#import scipy

import OutputDef as Output

# Set file
# =====================
rootFolder = "/Users/abauville/GoogleDrive/StokesFD_Output/OutputTest/"
outFolder = "Out_00001/"
dataFolder = rootFolder + outFolder




# Read data
# =====================
state = Output.getModelState(dataFolder + '/modelState.json')
khi = Output.getData(dataFolder + 'Khi.bin')




# Define grid
# =====================
thisData = khi
x = np.linspace(thisData.xmin,thisData.xmax,thisData.nx)
y = np.linspace(thisData.ymin,thisData.ymax,thisData.ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)




# Plot
# =====================
data = np.transpose(thisData.data)
plt.pcolor(xv,yv,np.log10(data))
plt.axis('equal')
plt.show()




