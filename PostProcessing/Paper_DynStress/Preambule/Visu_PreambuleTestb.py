# Output reading test

import sys
import os
sys.path.insert(0, './../../../src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
#from pprint import pprint

#import scipy
import json

import OutputDef as Output

import InputDef as Input

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

from math import pi
from numpy import sin, tan, cos

import matplotlib as mpl


##             Units
## =====================================
m       = 1.0
s       = 1.0
K       = 1.0
kg      = 1.0

cm      = 0.01      * m
km      = 1000.0    * m

mn      = 60        * s
hour    = 60        * mn
day     = 24        * hour
yr      = 365       * day
Kyr     = 1e3       * yr
Myr     = 1e6       * yr

Pa      = kg/m/s/s
MPa     = 1e6 * Pa
GPa     = 1e9 * Pa

degree = pi/180

screenAdjust =  1.1 # somehow the figure is created a bit too small on my screen
cm2inch = 0.393701 * screenAdjust

dpi = 200 # same as default
cm2pt = cm2inch * dpi



# Colormap
# =====================
# Colormap Strain rate
# c-b-k-r-y
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5 , 0.0, 0.0),
                   (0.75, 1.0, 1.0),
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

cdict2 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5 , 0.0, 0.0),
                   (0.75, 1.0, 1.0),
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

## Colormap Pressure
## w-k-b-r
#cdict1 = {'red':  ((0.0 , 1.0, 1.0),
#                   (0.33, 0.1, 0.1),
#                   (0.66, 0.0, 0.0),
#                   (1.0 , 1.0, 0.0)),
#
#         'green': ((0.0 , 1.0, 1.0),
#                   (0.33, 0.1, 0.1),
#                   (0.66, .25, .25),
#                   (1.0 , .00, 0.0)),
#
#         'blue':  ((0.0 , 1.0, 1.0),
#                   (0.33, 0.0, 0.0),
#                   (0.66, 1.0, 1.0),
#                   (1.0 , 0.0, 0.0))
#        }


## darkB-b-c-w-y-r-darkR
#cdict2 = {'red':  ((0.0     , 0.0 , 0.0),
#                   (1.0/6.0 , 0.0 , 0.0),
#                   (2.0/6.0 , 0.2 , 0.2),
#                   (3.0/6.0 , 1.0 , 1.0),
#                   (4.0/6.0 , 1.0 , 1.0),
#                   (5.0/6.0 , 1.0 , 1.0),
#                   (1.0 , 0.33, 1.0)),
#
#         'green': ((0.0     , 0.0 , 0.0),
#                   (1.0/6.0 , 0.0 , 0.0),
#                   (2.0/6.0 , 0.8 , 0.8),
#                   (3.0/6.0 , 1.0 , 1.0),
#                   (4.0/6.0 , 0.8 , 0.8),
#                   (5.0/6.0 , 0.0 , 0.0),
#                   (1.0     , 0.0 , 0.0)),
#
#         'blue':  ((0.0     , 0.33, 0.33),
#                   (1.0/6.0 , 1.0 , 1.0),
#                   (2.0/6.0 , 1.0 , 1.0),
#                   (3.0/6.0 , 1.0 , 1.0),
#                   (4.0/6.0 , 0.2 , 0.2),
#                   (5.0/6.0 , 0.0 , 0.0),
#                   (1.0     , 0.0 , 0.0))
#        }



## Cmap of strain for the thin-thick paper
#cdict2 = {'red':  ((0.0     , 1.0 , 1.0),
#                   (1.0/4.0 , 0.25, 0.25),
#                   (2.0/4.0 , 0.2 , 0.2),
#                   (3.0/4.0 ,  .9 , .9),
#                   (1.0     , 0.8 , .8)),
#
#         'green': ((0.0     , 1.0 , 1.0),
#                   (1.0/4.0 , 0.25, 0.25),
#                   (2.0/4.0 , 0.5 , 0.5),
#                   (3.0/4.0 ,  .8 , .8),
#                   (1.0     , 0.1 , 0.1)),
#
#         'blue':  ((0.0     ,1.0 , 1.0),
#                   (1.0/4.0 , .7 ,  .7),
#                   (2.0/4.0 , .5 , .5),
#                   (3.0/4.0 , .0 , .0),
#                   (1.0     , .1 , .1))
#        }


## Cmap of strain for the thin-thick paper -- more contrast
#cdict2 = {'red':  ((0.0     , 0.0 , 0.0),
#                   (1.0/4.0 , 0.25, 0.25),
#                   (2.0/4.0 , 0.25 , 0.25),
#                   (3.0/4.0 , .95 , .95),
#                   (1.0     , 1.0 , 1.0)),
#
#         'green': ((0.0     , 0.0 , 0.0),
#                   (1.0/4.0 , 0.25, 0.25),
#                   (2.0/4.0 , 0.5 , 0.5),
#                   (3.0/4.0 ,  .85, .85),
#                   (1.0     , 0.0 , 0.0)),
#
#         'blue':  ((0.0     ,0.0 , 0.0),
#                   (1.0/4.0 , .75 ,  .75),
#                   (2.0/4.0 , .5 , .5),
#                   (3.0/4.0 , .0 , .0),
#                   (1.0     , .0 , .0))
#        }

# w-blueish-y-oramge
#cdict2 = {'red':  ((0.0     , 1.0 , 1.0),
#                   (1.0/3.0 , 0.25, 0.25),
#                   (2.0/3.0 ,  .9 , .9),
#                   (1.0     , 0.8 , .8)),
#
#         'green': ((0.0     , 1.0 , 1.0),
#                   (1.0/3.0 , 0.25, 0.25),
#                   (2.0/3.0 ,  .8 , .8),
#                   (1.0     , 0.1 , 0.1)),
#
#         'blue':  ((0.0     ,1.0 , 1.0),
#                   (1.0/3.0 , .7 ,  .7),
#                   (2.0/3.0 , .0 , .0),
#                   (1.0     , .1 , .1))
#        }


#cdict2 = {'red':  ((0.0     , 1.0 , 1.0),
#                   (1.0/3.0 ,  .9 , .9),
#                   (2.0/3.0 , 0.25, 0.25),
#                   (1.0     , 1.0 , 1.0)),
#
#         'green': ((0.0     , 1.0 , 1.0),
#                   (1.0/3.0 ,  .8 , .8),
#                   (2.0/3.0 , 0.25, 0.25),
#                   (1.0     , 0.1 , 0.1)),
#
#         'blue':  ((0.0     ,1.0 , 1.0),
#                   (1.0/3.0 , .0 , .0),
#                   (2.0/3.0 , .7 ,  .7),
#                   (1.0     , .1 , .1))
#        }




#
#
#cdict1 = {'red':  ((0.0 , 0.0, 0.0),
#                   (0.50, 0.0, 0.0),
#                   (1.0 , 0.0, 0.0)),
#
#         'green': ((0.0 , 0.0, 0.0),
#                   (0.50, 1.0, 1.0),
#                   (1.0 , 0.0, 0.0)),
#
#         'blue':  ((0.0 , 0.0, 0.0),
#                   (0.5 , 0.0, 0.0),
#                   (1.0 , 1.0, 1.0))
#        }

CMAP = LinearSegmentedColormap('StrainRate', cdict1)
plt.register_cmap(cmap=CMAP)
CMAP = LinearSegmentedColormap('Pressure', cdict2)
plt.register_cmap(cmap=CMAP)
plt.set_cmap('StrainRate')




# Set file
# =====================
RFac = 1
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_Test/"
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test/dt_stressFac_1.0e-05/"
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_Stronger_Seed_10timesWeaker/dt_stressFac_1.0e-03/"
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_Stronger_Seed_Save/dt_stressFac_1.0e-03/"
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/dtDependence/Test_CleanSave/dt_stressFac_1.0e-03/"
simFolder  = ""
inFolder  = "Input/"
nSteps = Output.getNumberOfOutFolders(rootFolder) -1;
#nSteps = 150
outFolder = "Out_%05d/" % (nSteps-1)


# Read parameters of this simulation
# =====================
#inputFile = Output.readJson(rootFolder + simFolder + inFolder + '/input.json')
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char
CharExtra = Input.CharExtra(Char)


# Read strain rate data
# =====================
#dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'khi.bin')
#khi = dataSet.data

#dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
#P = dataSet.data

dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
SII = dataSet.data

C = Setup.MatProps['1'].cohesion
phi = Setup.MatProps['1'].frictionAngle

backStrainRate = abs(Setup.BC.Stokes.backStrainRate)

xmin    = dataSet.xmin * Setup.Char.length
xmax    = dataSet.xmax * Setup.Char.length
ymin    = dataSet.ymin * Setup.Char.length
ymax    = dataSet.ymax * Setup.Char.length
nx      = dataSet.nx
ny      = dataSet.ny
Wbox       = xmax-xmin 
Hbox       = ymax-ymin
dx      = Wbox/nx
dy      = Hbox/ny

# Define grid
# =====================
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)

xv, yv = np.meshgrid(x,y)
xv = np.transpose(xv)
yv = np.transpose(yv)

Hmatrix = 1000.0;




# Set figure
# =====================
pageW       = 21.0      * cm2inch
pageH       = 25.0      * cm2inch #29.7*cm2inch
pageLPad    = 2.0       * cm2inch
pageRPad    = 3.0       * cm2inch
pagePad     = pageLPad+pageRPad

nFields      = 3 # Strain rate Pressure, Stress

nSim        = 3
simWPad     = 0.5       * cm2inch
simHPad     = .75       * cm2inch
simW        = ((pageW-pagePad) - (nSim-1)*simWPad) / nSim                      # width
simH        = simW * Hbox/Wbox                                                  # height
simB        = pageH-pagePad-simH                                                # bottom



# Vertical colorbar
colorBarH   =  simH
colorBarHPad=  .0       * cm2inch
colorBarWPad= .25       * cm2inch
colorBarW   = .3        * cm2inch
colorBarB = simB



graphH      = 5.0       * cm2inch
graphLPad   = .8       * cm2inch
graphW      = pageW - pagePad - graphLPad + colorBarWPad + colorBarW - simW - graphLPad - .25*cm2inch
graphHPad   = simHPad#1.0       * cm2inch
graphB      = simB-2*simH-2*simHPad - graphHPad - graphH

zoomGraphH      = 5.0       * cm2inch
#zoomGraphLPad   = 1.0       * cm2inch
zoomGraphW      = simW
zoomGraphWPad   = simWPad + colorBarWPad + colorBarW + -.2 * cm2inch
zoomGraphHPad   = simHPad#1.0       * cm2inch
zoomGraphB      = simB-2*simH-2*simHPad - graphHPad - graphH




thisFig = plt.figure(3)#,figsize = (pageW,pageH))
thisFig.set_size_inches(pageW,pageH)



#thisFig = plt.figure(1,figsize = (0.1,1))
plt.clf()
axSimEII = dict()
axSimP   = dict()
axSimSII = dict()
#letters = "abcdefghijklmno"
letters = "ABCDEFGHIJKLMNO"
for iSim in range(0,nSim):
    axSimEII["%i" % iSim] = plt.axes([(pageLPad+iSim*simW+iSim*simWPad)/pageW,simB/pageH,simW/pageW,simH/pageH])
    axSimEII["%i" % iSim].set_xticks([])
    axSimEII["%i" % iSim].set_yticks([])
#    plt.text(xmin+Wbox/30.0,ymin+Hbox/30.0,letters[iSim])
    
for iSim in range(0,nSim):
    axSimP  ["%i" % iSim] = plt.axes([(pageLPad+iSim*simW+iSim*simWPad)/pageW,(simB-simHPad-simH)/pageH,simW/pageW,simH/pageH])
    axSimP  ["%i" % iSim].set_xticks([])
    axSimP  ["%i" % iSim].set_yticks([])
    
for iSim in range(0,nSim):
    axSimSII["%i" % iSim] = plt.axes([(pageLPad+iSim*simW+iSim*simWPad)/pageW,(simB-2*simHPad-2*simH)/pageH,simW/pageW,simH/pageH])
    axSimSII["%i" % iSim].set_xticks([])
    axSimSII["%i" % iSim].set_yticks([])
    
#axColorbar = plt.axes([(pagePad+colorBarWPad)/pageW,colorBarB/pageH,colorBarW/pageW,colorBarH/pageH])

axGraph = plt.axes([(pageLPad+graphLPad)/pageW, graphB/pageH, graphW/pageW, graphH/pageH])
axZoomGraph = plt.axes([(pageLPad+graphLPad+graphW+zoomGraphWPad)/pageW, zoomGraphB/pageH, zoomGraphW/pageW, zoomGraphH/pageH])
axGraph.spines['right'].set_visible(False)
axGraph.spines['top'].set_visible(False)
axZoomGraph.spines['right'].set_visible(False)
axZoomGraph.spines['top'].set_visible(False)


axGraphTitle = plt.axes([(pageLPad)/pageW, graphB/pageH, graphW/pageW, graphH/pageH],frameon=False,faceColor='None')
plt.xticks([])
plt.yticks([])

#dictOption = dict(loc='left',fontName='Times New Roman')
plt.sca(axSimEII["0"])
plt.title("1) Strain rate",loc='left',fontName='Times New Roman',VerticalAlignment='center') 
plt.sca(axSimP  ["0"])
plt.title("2) Pressure",loc='left',fontName='Times New Roman',VerticalAlignment='center') 
plt.sca(axSimSII["0"])
plt.title("3) Deviatoric stress",loc='left',fontName='Times New Roman',VerticalAlignment='center') 
plt.sca(axGraphTitle)
plt.title("4) Stress evolution",loc='left',fontName='Times New Roman',VerticalAlignment='center') 











# Extracting data
# =====================
computePostProc = True
if computePostProc:
    # Choose a cell to monitor
    # =====================
    halfSpan = 76
    ixCellMin = 0#76 - halfSpan
    ixCellMax = nx#ixCellMin + 2*halfSpan
    iyCell = np.argmin( np.abs(y-Hmatrix/2.0) )
    ixCell = int(round((ixCellMin+ixCellMax)/2.0))
    
    # Declare time series arrays
    # =====================
    SII_t = np.zeros(nSteps) # Second invariant of the deviatoric stress tensor
    P_t = np.zeros(nSteps)
    time_t = np.zeros(nSteps)
    dt_t = np.zeros(nSteps)
    Sy_t = np.zeros(nSteps) # yield stress
    I_t = np.zeros(nSteps,dtype=int)
    
    # Loop over time
    # =====================
    it = -1
    for iFolder in range(0,nSteps,1) :
        it += 1
        if (np.mod(iFolder,100)==0):
                print("iStep = %i/%i" % (iFolder, nSteps))
        outFolder   = "Out_%05d/" % (iFolder)
        State       = Output.readState(rootFolder + simFolder + outFolder + "modelState.json")
    
        dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
        SII = dataSet.data * CharExtra.stress
        subset_SII  = SII[ixCellMin:ixCellMax+1,iyCell]
        
        dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
        P = dataSet.data * CharExtra.stress
        subset_P  = P[ixCellMin:ixCellMax+1,iyCell]
        
#        I = np.argmin(subset_P) # find the minimum khi
        I = 68
        SII_t[it] = subset_SII[I] # get the stress corresponding to the minimum value
        
        P_t[it] = subset_P[I] # get the stress corresponding to the minimum value
        Sy_t[it] = C*np.cos(phi) + P_t[it]*np.sin(phi)
        
        time_t[it] = State.time * Setup.Char.time + State.dt
        dt_t[it] = State.dt #* Char.time
        
        I_t[it]     = I
    #endfor
    
#    np.savez("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/PreambuleTest",SII_t=SII_t,P_t=P_t,Sy_t=Sy_t,time_t=time_t,dt_t=dt_t,I_t=I_t,ixCell=ixCell,iyCell=iyCell)
else:
    loadedData = np.load("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/PreambuleTest.npz");
#    loadedData = np.load("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Save/Preambule.npz");
    SII_t   = loadedData["SII_t"]
    P_t     = loadedData["P_t"]
    Sy_t    = loadedData["Sy_t"]
    time_t  = loadedData["time_t"]
    dt_t    = loadedData["dt_t"]
    I_t     = loadedData["I_t"]
    ixCell  = loadedData["ixCell"]
    iyCell  = loadedData["iyCell"]
    
#endif





# Analytical yield stress
# =====================  
#P = Setup.Physics.Pback
#Sy_back = C*np.cos(phi) + P*np.sin(phi)
S3 = Setup.Physics.Pback
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
Sy_back = (S1-S3)/2.0
P_Lim = (S1+S3)/2.0

# Find interesting values
# =====================
I_sigmaMax = np.argmax(SII_t)
#I_sigmaMin = I_sigmaMax + np.argmin(SII_t[I_sigmaMax:])
I_sigmaMin = I_sigmaMax + np.argmin(SII_t[I_sigmaMax:I_sigmaMax+200])

Delta_timeSoft = time_t[I_sigmaMin] - time_t[I_sigmaMax]


# Units and characteristic values
# =====================
EII = backStrainRate
eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
G = Setup.MatProps['1'].G
#t = 0
t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
charTime = t
#charTime = Sy_back / (2*G*EII * np.exp(-G/eta*t));
timeMaxwell = eta/G
timePlotUnit = charTime#timeMaxwell#Kyr

Pback = Setup.Physics.Pback
sigmaPlotUnit = Pback#MPa











# Plotting 2D
# =====================
StepList = [int(nSteps/8.0),I_sigmaMax,I_sigmaMin]
for iSim in range(0,nSim):
    outFolder = "Out_%05d/" % (StepList[iSim]-1)
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'strainRate.bin')
    EII = dataSet.data * 1.0/Char.time
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
    P = dataSet.data * CharExtra.stress
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
    SII = dataSet.data * CharExtra.stress
    cAx_EIIMin = -2.0
    cAx_EIIMax = +2.0
    plt.sca(axSimEII["%i" % iSim])

    plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(EII / backStrainRate),vmin=cAx_EIIMin, vmax=cAx_EIIMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
    plt.axis('tight')
    plt.set_cmap("StrainRate")
    
    cAx_PMin = 1.0
    cAx_PMax = 3.0
    plt.sca(axSimP  ["%i" % iSim])
    plt.pcolor(xv - dx/2.0,yv - dy/2.0,P/sigmaPlotUnit,vmin=cAx_PMin,vmax=cAx_PMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
    plt.axis('tight')
    plt.set_cmap("Pressure")

    cAx_SIIMin = 1.0
    cAx_SIIMax = 1.75
    plt.sca(axSimSII["%i" % iSim])
    plt.pcolor(xv - dx/2.0,yv - dy/2.0,SII/sigmaPlotUnit,vmin=cAx_SIIMin,vmax=cAx_SIIMax) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
    plt.set_cmap("Pressure")
    plt.axis('equal')

    
# Add index of subfigures
# =====================
sizeFacBox = 6.0
dictOption = dict(color='w',horizontalAlignment='left',verticalAlignment='baseline',fontName='PT Mono',fontWeight='bold')
xText = xmin+Hbox/sizeFacBox/6.0
yText = ymax-Hbox/sizeFacBox/1.25
for iSim in range(0,nSim):
     plt.sca(axSimEII["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox*2.0,xmin+Hbox/sizeFacBox*2.0]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+0*nSim] + ". $t_%i$" % iSim,dictOption)
     
     plt.sca(axSimP  ["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox,xmin+Hbox/sizeFacBox]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+1*nSim],dictOption)
     
     plt.sca(axSimSII["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox,xmin+Hbox/sizeFacBox]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+2*nSim],dictOption)



# Colorbar
# =====================
# Colorbar EII
# ===============     

axColorbarEII = plt.axes([(pageW-pageRPad+colorBarWPad)/pageW,colorBarB/pageH,colorBarW/pageW,colorBarH/pageH])
axColorbarEII.yaxis.label.set_fontname("Courier")

List_cAx = np.arange(cAx_EIIMin,cAx_EIIMax+1)
List_cAx_10pow = []
for i in range(0,List_cAx.size):
    List_cAx_10pow.append("$10^{%i}$" % (List_cAx[i]))

axColorbarEII.xaxis.set_label_position('top')
axColorbarEII.set_xlabel("$\dot{\\epsilon}_{II}/\dot{\\epsilon}_{back}$",verticalAlignment='center',labelpad=.12*cm2pt)

plt.sca(axSimEII["0"])
Cbar = plt.colorbar(cax=axColorbarEII, orientation='vertical')
Cbar.set_ticks(List_cAx)

axColorbarEII.tick_params(direction='in')     
axColorbarEII.set_yticklabels(List_cAx_10pow)

# Colorbar P
# ===============    
axColorbarP   = plt.axes([(pageW-pageRPad+colorBarWPad)/pageW,(colorBarB-simH-simHPad)/pageH,colorBarW/pageW,colorBarH/pageH])
axColorbarP.yaxis.label.set_fontname("Courier")

axColorbarP.xaxis.set_label_position('top')
axColorbarP.set_xlabel("$P/P_{back}$",verticalAlignment='center',labelpad=.12*cm2pt)
plt.sca(axSimP  ["0"])
Cbar = plt.colorbar(cax=axColorbarP, orientation='vertical')
axColorbarP.tick_params(direction='in')     
List_cAx = np.arange(cAx_PMin,cAx_PMax+1)
Cbar.set_ticks(List_cAx)

# Colorbar SII
# ===============    
axColorbarSII   = plt.axes([(pageW-pageRPad+colorBarWPad)/pageW,(colorBarB-2*simH-2*simHPad)/pageH,colorBarW/pageW,colorBarH/pageH])
axColorbarSII.yaxis.label.set_fontname("Courier")

axColorbarSII.xaxis.set_label_position('top')
axColorbarSII.set_xlabel("$\\tau_{II}/P_{back}$",verticalAlignment='center',labelpad=.12*cm2pt)
axColorbarSII.tick_params(direction='in')     
plt.sca(axSimSII  ["0"])
Cbar = plt.colorbar(cax=axColorbarSII, orientation='vertical')
List_cAx = np.linspace(cAx_SIIMin,cAx_SIIMax,4)
ListLabel_cAx = ["%.1f" % number for number in List_cAx]
ListLabel_cAx[1] = ""
ListLabel_cAx[3] = ""
Cbar.set_ticks(List_cAx)
plt.sca(axColorbarSII)
plt.yticks(List_cAx,ListLabel_cAx)












# Plotting Zoom Graph
# =====================
plt.sca(axZoomGraph)
zoomCoord = np.array([1.6, 1.7001, 1.0, np.max(P_t/sigmaPlotUnit)*1.05 ])
plt.fill(zoomCoord[[0,0,1,1]],zoomCoord[[2,3,3,2]],color=[.9,.9,.9])
plt.plot(time_t/timePlotUnit,(P_t)/sigmaPlotUnit,'-b')
plt.plot(time_t/timePlotUnit,Sy_t/sigmaPlotUnit,'--k',linewidth=4)
plt.plot(time_t/timePlotUnit,SII_t/sigmaPlotUnit,'-r')

plt.plot((0,time_t[-1]/timePlotUnit), np.array((Sy_back,Sy_back))/sigmaPlotUnit,'--k')


#plt.legend(["$\\tau_{y}$","$\\tau_{II}$","$P$","$\\tau_{y}$ at $P_{back}$"])
plt.xlabel("$time/time_c$")

#plt.ylabel("Stress/$P_{back}$")


plt.axis(zoomCoord)

#plt.tick_params(axis='y',which='both',pad=.5)
yTickList = np.arange(zoomCoord[2],zoomCoord[3],.25)
yTickLabelList  = ["%.1f" % number for number in yTickList]
#for i in range(len(yTickLabelList)):
#    if round(yTickList[i]) != yTickList[i]:
#        yTickLabelList[i] = ""
yTickLabelList[1::2] = [""]*len(yTickLabelList[1::2])
plt.yticks(yTickList,yTickLabelList)

xTickList = np.arange(zoomCoord[0],zoomCoord[1],.05)
xTickLabelList = ["%.1f" % number for number in xTickList]
xTickList = np.append(xTickList, time_t[StepList[1:]]/timePlotUnit)
xTickLabelList.extend(("", ""))
plt.text(time_t[StepList[1]]/timePlotUnit,zoomCoord[2]+(zoomCoord[3]-zoomCoord[2])/20.0,"$t_1$",horizontalAlignment='center')
plt.text(time_t[StepList[2]]/timePlotUnit,zoomCoord[2]+(zoomCoord[3]-zoomCoord[2])/20.0,"$t_2$",horizontalAlignment='center')
#plt.xaxis.set_minor_locator(minorLocator)
plt.xticks(xTickList,xTickLabelList)



plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='minor',      # both major and minor ticks are affected
    direction='in') # labels along the bottom edge are of

print("Delta_timeSoft = %.f Kyr" % (Delta_timeSoft/Kyr))
print("sigmaMax = %.f MPa, sigmaMin = %.f MPa" % (SII_t[I_sigmaMax]/MPa, SII_t[I_sigmaMin]/MPa))
print("sigmaMaxTime = %.f Kyr, sigmaMinTime = %.f Kyr" % (time_t[I_sigmaMax]/Kyr, time_t[I_sigmaMin]/Kyr))




# Plotting Graph
# =====================
plt.sca(axGraph)

plt.fill(zoomCoord[[0,0,1,1]],zoomCoord[[2,3,3,2]],color=[.9,.9,.9])

plt.plot(time_t/timePlotUnit,(P_t)/sigmaPlotUnit,'.b')
plt.plot(time_t/timePlotUnit,Sy_t/sigmaPlotUnit,'--k',linewidth=4)
plt.plot(time_t/timePlotUnit,SII_t/sigmaPlotUnit,'-r')

plt.plot((0,time_t[-1]/timePlotUnit), np.array((Sy_back,Sy_back))/sigmaPlotUnit,'--k')


plt.legend(["$P$","$\\tau_{y}$","$\\tau_{II}$","$\\tau_{yback}$"],ncol=2,frameon=False,loc='upper left',bbox_to_anchor=(.05,1.03))
plt.xlabel("$time/time_c$")
plt.ylabel("Stress/$P_{back}$")

graphCoord = [0.0, time_t[-1]/timePlotUnit, 0.0, np.max(P_t/sigmaPlotUnit)*1.05 ]

plt.axis(graphCoord)
yTickList = np.arange(graphCoord[2],graphCoord[3],.5)
#yTickLabelList = np.array2string(yTickList)
yTickLabelList  = ["%.0f" % number for number in yTickList]
yTickLabelList[1::2] = [""]*len(yTickLabelList[1::2])
plt.yticks(yTickList,yTickLabelList)
xTickList = np.arange(graphCoord[0],graphCoord[1],.5)
xTickLabelList = ["%.1f" % number for number in xTickList]
xTickList = np.append(xTickList, time_t[StepList]/timePlotUnit)
xTickLabelList.extend(("", "", ""))
plt.text(time_t[StepList[0]]/timePlotUnit,graphCoord[2]+(graphCoord[3]-graphCoord[2])/20.0,"$t_0$",horizontalAlignment='center')
plt.text(time_t[StepList[1]]/timePlotUnit,graphCoord[2]+(graphCoord[3]-graphCoord[2])/20.0,"$t_1$",horizontalAlignment='right')
plt.text(time_t[StepList[2]]/timePlotUnit,graphCoord[2]+(graphCoord[3]-graphCoord[2])/20.0,"$t_2$",horizontalAlignment='left')
plt.xticks(xTickList,xTickLabelList)
plt.tick_params(axis='both',direction='in')

# Add index of subfigures
# =====================
plt.sca(axGraph)
dictOption = dict(color='w',horizontalAlignment='center',verticalAlignment='center',fontName='PT Mono',fontWeight='bold')
sizeFacBox = 10.0

#HWRatio = (graphCoord[3]-graphCoord[2])/(graphCoord[1]-graphCoord[0])
HWRatio = graphH/graphW
x0 = graphCoord[0]
x1 = graphCoord[0]+(graphCoord[1]-graphCoord[0])/sizeFacBox*HWRatio
y0 = graphCoord[3]
y1 = graphCoord[3]-(graphCoord[3]-graphCoord[2])/sizeFacBox
xText = (x0+x1)/2
yText = (y0+y1)/2
plt.fill([x0,x0,x1,x1],[y0,y1,y1,y0],'k')
plt.text(xText,yText,letters[iSim+2*nSim+1],dictOption)

plt.sca(axGraph)
dictOption = dict(color='w',horizontalAlignment='center',verticalAlignment='center',fontName='PT Mono',fontWeight='bold')
sizeFacBox = 10.0

#HWRatio = (graphCoord[3]-graphCoord[2])/(graphCoord[1]-graphCoord[0])
plt.sca(axZoomGraph)
HWRatio = zoomGraphH/zoomGraphW
x0 = zoomCoord[0]
x1 = zoomCoord[0]+(zoomCoord[1]-zoomCoord[0])/sizeFacBox*HWRatio
y0 = zoomCoord[3]
y1 = zoomCoord[3]-(zoomCoord[3]-zoomCoord[2])/sizeFacBox
xText = (x0+x1)/2
yText = (y0+y1)/2
plt.fill([x0,x0,x1,x1],[y0,y1,y1,y0],'k')
plt.text(xText,yText,letters[iSim+2*nSim+2],dictOption)

plt.tick_params(axis='both',direction='in')

print("Delta_timeSoft = %.f Kyr" % (Delta_timeSoft/Kyr))
print("sigmaMax = %.f MPa, sigmaMin = %.f MPa" % (SII_t[I_sigmaMax]/MPa, SII_t[I_sigmaMin]/MPa))
print("sigmaMaxTime = %.f Kyr, sigmaMinTime = %.f Kyr" % (time_t[I_sigmaMax]/Kyr, time_t[I_sigmaMin]/Kyr))

#plt.savefig("/Users/abauville/Dropbox/01_Papers/DynStressPaper/Figures/Preambule.png")

