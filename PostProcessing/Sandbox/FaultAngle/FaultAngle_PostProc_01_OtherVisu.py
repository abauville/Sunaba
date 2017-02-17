# Output reading test

import sys
sys.path.insert(0, '../../../src/UserInput')
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
rootFolder = "/Users/abauville/StokesFD_Output/FaultAngle_Systematics_Batch03/"
simFolder  = "HSed1_SA0_C1000_FA4/"
inFolder  = "Input/"
outFolder = "Out_00000/"


# Read parameters of this simulation
# =====================
inputFile = Output.readJson(rootFolder + simFolder + inFolder + '/input.json')
s = inputFile["Description"]
json_acceptable_string = s.replace("'", "\"")
Params = json.loads(json_acceptable_string)

Hsed            = Params["Hsed"]
cohesion        = Params["cohesion"]
frictionAngle   = Params["friction_angle"]
surfaceAngle    = Params["surface_angle"]

HFac            = Hsed/(1.0*km)

Hsed_nondim     = Hsed / inputFile["Char"]["length"]


# Read strain rate data
# =====================

Eii     = Output.getData(rootFolder + simFolder + outFolder + 'strainRate.bin')
xmin    = Eii.xmin
xmax    = Eii.xmax
ymin    = Eii.ymin
ymax    = Eii.ymax
W       = xmax-xmin
H       = Hsed_nondim

# Define grid
# =====================
thisData = Eii
x = np.linspace(xmin,xmax,thisData.nx)
y = np.linspace(ymin,ymax,thisData.ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)
















#Syst_HFac               = [0.1, 1.0, 10.0]
#Syst_surfaceAngleAngle  = [0, 2, 4, 6, 8, 10, 12] # in degrees
#Syst_frictionAngle      = [2, 4, 8, 16, 32]    # in degrees
#Syst_cohesion           = [0.1*MPa, 1.0*MPa, 10.0*MPa, 100.0*MPa]



Syst_HSed               = [1]
Syst_surfaceAngleAngle  = [2] # in degrees
Syst_frictionAngle      = [4]    # in degrees
Syst_cohesion           = [1.0*MPa]


no_simulation = len(Syst_HSed) * len(Syst_surfaceAngleAngle) * len(Syst_frictionAngle) * len(Syst_cohesion)
                
#print("number of simulations: " + str(no_simulation))

# counters



#thisHFac        = Syst_HFac[1]
#thisCohesion    = Syst_cohesion[2]



WShift = 1.1*W
HShift = 1.1*H

xShiftIni = -xmin
yShiftIni = -ymin

simCounter = 0

for thisHSed in Syst_HSed:
    
    for thisCohesion in Syst_cohesion:
        plt.clf()
        
        plt.axis('equal')
        print("\nHSed: " + str(thisHSed) + ", C: " + str(thisCohesion) )
        
        ix = -1
        iy = -1
        
        for thisSurfaceAngle in Syst_surfaceAngleAngle:
            iy+=1
            ix = -1
            for thisFrictionAngle in Syst_frictionAngle:
                simCounter += 1
                ix+=1
                print("   sim #: " + str(simCounter) + "/" + str(no_simulation) + ", SA: " + str(thisSurfaceAngle) + ", FA: " + str(thisFrictionAngle))
                simFolder = "HSed" + str(round(thisHSed)) + "_SA" + str(round(thisSurfaceAngle)) + "_C" + str(round(thisCohesion)) + "_FA" + str(round(thisFrictionAngle)) + "/"
                # Shifting
                # =====================
                xShift = WShift*ix + xShiftIni
                yShift = HShift*iy + yShiftIni
                
                #try:



                # Get deviatoric stresses and Pressure
                # =====================
                thisData = Output.getData(rootFolder + simFolder + outFolder + 'sigma_xx.bin')
                Tau_xx = np.ma.masked_array(thisData.data, yv > (Hsed_nondim+tan(+thisSurfaceAngle*degree)*xv)) # Hide the air
                                         
                                         
                thisData = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
                P = np.ma.masked_array(thisData.data, yv > (Hsed_nondim+tan(+thisSurfaceAngle*degree)*xv)) # Hide the air
                
                thisData = Output.getData(rootFolder + simFolder + outFolder + 'sigma_xy.bin')
                thisData.interpFromNodesToCells()
                Sxy = np.ma.masked_array(thisData.data, yv > (Hsed_nondim+tan(+thisSurfaceAngle*degree)*xv)) # Hide the air
                                         
                Sxx =  Tau_xx - P
                Syy = -Tau_xx - P
                Syx = Sxy                         
                                        
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                # Plot
                # =====================        
                #plt.pcolor(xv+xShift,yv+yShift,np.log10(data), vmin=0, vmax=2.0)
#                plt.pcolor(xv+xShift,yv+yShift,(Sxx),vmin=-10000,vmax=10000)
                plt.pcolor(xv+xShift,yv+yShift,np.log10(np.abs(Tau_xx/Sxy)),vmin=-1, vmax=2)
                plt.plot([xmin, xmax, xmax, xmin, xmin]+xShift , [ymin, ymin, Hsed_nondim, Hsed_nondim+tan(+thisSurfaceAngle*degree)*xmin, ymin]+yShift, "k", linewidth=0.5)
#                except:
#                    simFolder = "HFac" + str(round(thisHFac)) + "_SA" + str(round(thisSurfaceAngle)) + "_C" + str(round(thisCohesion)) + "_FA" + str(round(thisFrictionAngle)) + "/"
#                    print(simFolder + "could not be loaded")
                
        
         
        
        plt.yticks(HShift*np.arange(iy+1) + yShiftIni + H/5.0, Syst_surfaceAngleAngle)
        plt.xticks(WShift*np.arange(ix+1) + xShiftIni - W/2.0, Syst_frictionAngle)
        
        plt.ylabel("surface angle [°]")
        plt.xlabel("friction angle [°]")
        
        plt.title( "Hsed = " + str(thisHSed) + " m" + ", Cohesion = " + str(thisCohesion/1e6) + " MPa")
        
        plt.axis([(W-WShift), (ix)*WShift+xShiftIni+(WShift-W), (H-HShift), (iy+1)*HShift+yShiftIni+(HShift-H)])
        plt.colorbar()
        plt.show()
        
        
        
        
        
        
        #plt.savefig("HFac" + str(round(thisHFac)) + "_C" + str(round(thisCohesion)) + ".png", dpi=500)
        
        













