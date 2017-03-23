# Output reading test

import sys
sys.path.insert(0, '../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
#from pprint import pprint

#import scipy
import json

import OutputDef as Output

import InputDef as Input

from matplotlib.colors import LinearSegmentedColormap

from math import pi, tan



# Units
# =====================
degree = pi/180
MPa = 1e6
km = 1e3


# Colormap
# =====================
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.75 , 1.0, 1.0),
                   (1.0 , 1.0, 1.0)),

         'green': ((0.0 , 1.0, 1.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.0 , 1.0, 1.0)),

         'blue':  ((0.0 , 1.0, 1.0),
                   (0.25 , 1.0, 1.0),
                   (0.5 , 0.0, 0.0),
                   (1.0 , 0.0, 0.0))
        }


CMAP = LinearSegmentedColormap('StokesFD', cdict1)
plt.register_cmap(cmap=CMAP)
plt.set_cmap('StokesFD')
#plt.set_cmap("gray")



# Set file
# =====================
rootFolder = "/Users/abauville/StokesFD_Output/TestDarcy/"
simFolder  = ""
inFolder  = "Input/"
outFolder = "Out_01300/"
#outFolder = "Out_00000/"


# Read parameters of this simulation
# =====================
inputFile = Output.readJson(rootFolder + simFolder + inFolder + '/input.json')
s = inputFile["Description"]

# Read strain rate data
# =====================

#dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'Porosity.bin')
#vmin = .325
#vmax = .375

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'Pf.bin')
dataSet2     = Output.getData(rootFolder + simFolder + outFolder + 'Pc.bin')

#dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'Pf.bin')

vmin = -0.1
vmax = 0.1
#data = dataSet.data
data = dataSet.data + dataSet2.data

xmin    = dataSet.xmin
xmax    = dataSet.xmax
ymin    = dataSet.ymin
ymax    = dataSet.ymax
W       = xmax-xmin
H       = ymax-ymin

# Define grid
# =====================
x = np.linspace(xmin,xmax,dataSet.nx)
y = np.linspace(ymin,ymax,dataSet.ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)



# Shifting
# =====================
xShift = 0#WShift*ix + xShiftIni
yShift = 0#HShift*iy + yShiftIni

#try:


# =====================        
#plt.pcolor(xv+xShift,yv+yShift,np.log10(data), vmin=0, vmax=2.0)
#                plt.pcolor(xv+xShift,yv+yShift,(Sxx),vmin=-10000,vmax=10000)
plt.clf()
plt.pcolor(xv+xShift,yv+yShift,data,vmin=vmin, vmax=vmax)
#plt.pcolor(xv+xShift,yv+yShift,dataSet.data)
plt.axis('equal')
plt.colorbar()
#plt.plot([xmin, xmax, xmax, xmin, xmin]+xShift , [ymin, ymin, Hsed_nondim, Hsed_nondim+tan(+thisSurfaceAngle*degree)*xmin, ymin]+yShift, "k", linewidth=0.5)
#                except:
#                    simFolder = "HFac" + str(round(thisHFac)) + "_SA" + str(round(thisSurfaceAngle)) + "_C" + str(round(thisCohesion)) + "_FA" + str(round(thisFrictionAngle)) + "/"
#                    print(simFolder + "could not be loaded")


 

#        plt.yticks(HShift*np.arange(iy+1) + yShiftIni + H/5.0, Syst_surfaceAngleAngle)
#        plt.xticks(WShift*np.arange(ix+1) + xShiftIni - W/2.0, Syst_frictionAngle)
#        
#        plt.ylabel("surface angle [°]")
#        plt.xlabel("friction angle [°]")
#        
#        plt.title( "Hsed = " + str(thisHSed) + " m" + ", Cohesion = " + str(thisCohesion/1e6) + " MPa")
#        
#        plt.axis([(W-WShift), (ix)*WShift+xShiftIni+(WShift-W), (H-HShift), (iy+1)*HShift+yShiftIni+(HShift-H)])
#        plt.colorbar()
#        plt.show()
#        
#        
#        



        #plt.savefig("HFac" + str(round(thisHFac)) + "_C" + str(round(thisCohesion)) + ".png", dpi=500)
        
        













