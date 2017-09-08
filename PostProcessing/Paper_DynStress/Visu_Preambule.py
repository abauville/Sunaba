# Output reading test

import sys
import os
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
RFac = 1
#rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_Test/"
rootFolder = "/Users/abauville/Work/Paper_DynStress/Output/Preambule_TestSave/"
simFolder  = ""
inFolder  = "Input/"
nSteps = Output.getNumberOfOutFolders(rootFolder);
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
pageW        = 21.0      * cm2inch
pageH        = 25.0      * cm2inch #29.7*cm2inch
pagePad     = 2.0       * cm2inch


nSim        = 3
simWPad     = 0.5       * cm2inch
simHPad     = 1.5       * cm2inch
simW        = ((pageW-2*pagePad) - (nSim-1)*simWPad) / nSim                      # width
simH        = simW * Hbox/Wbox                                                  # height
simB        = pageH-pagePad-simH                                                # bottom


colorBarH   = .4        * cm2inch
colorBarHPad= .25       * cm2inch
colorBarWPad= 2         * cm2inch
colorBarW   = (pageW-2*pagePad) - 2*colorBarWPad
colorBarB   = (simB-2*simHPad-2*simH-colorBarHPad-colorBarH)


graphH      = 5.0       * cm2inch
graphW      = pageW - 2*pagePad
graphHPad   = 1.5       * cm2inch
graphB      = colorBarB - graphHPad - graphH


thisFig = plt.figure(1)#,figsize = (pageW,pageH))
thisFig.set_size_inches(pageW,pageH)
#thisFig = plt.figure(1,figsize = (0.1,1))
plt.clf()
axSimEII = dict()
axSimP   = dict()
axSimSII = dict()
#letters = "abcdefghijklmno"
letters = "ABCDEFGHIJKLMNO"
for iSim in range(0,nSim):
    axSimEII["%i" % iSim] = plt.axes([(pagePad+iSim*simW+iSim*simWPad)/pageW,simB/pageH,simW/pageW,simH/pageH])
    axSimEII["%i" % iSim].set_xticks([])
    axSimEII["%i" % iSim].set_yticks([])
#    plt.text(xmin+Wbox/30.0,ymin+Hbox/30.0,letters[iSim])
    
for iSim in range(0,nSim):
    axSimP  ["%i" % iSim] = plt.axes([(pagePad+iSim*simW+iSim*simWPad)/pageW,(simB-simHPad-simH)/pageH,simW/pageW,simH/pageH])
    axSimP  ["%i" % iSim].set_xticks([])
    axSimP  ["%i" % iSim].set_yticks([])
    
for iSim in range(0,nSim):
    axSimSII["%i" % iSim] = plt.axes([(pagePad+iSim*simW+iSim*simWPad)/pageW,(simB-2*simHPad-2*simH)/pageH,simW/pageW,simH/pageH])
    axSimSII["%i" % iSim].set_xticks([])
    axSimSII["%i" % iSim].set_yticks([])
    
#axColorbar = plt.axes([(pagePad+colorBarWPad)/pageW,colorBarB/pageH,colorBarW/pageW,colorBarH/pageH])

axGraph = plt.axes([pagePad/pageW, graphB/pageH, graphW/pageW, graphH/pageH])


plt.sca(axSimEII["0"])
#dictOption = dict(loc='left',fontName='Times New Roman')
plt.title("1) Strain rate",loc='left',fontName='Times New Roman') 

    
#    
## Horiztonal colorbar
#colorBarH   =  .3       * cm2inch
#colorBarHPad=  .2       * cm2inch
#colorBarWPad= 4         * cm2inch
#colorBarW   = (pageW-2*pagePad) - 2*colorBarWPad
#colorBarB = simB + simH + colorBarHPad
##colorBarB   = (simB-2*simHPad-2*simH-colorBarHPad-colorBarH)
#axColorbar = plt.axes([(pagePad+colorBarWPad)/pageW,colorBarB/pageH,colorBarW/pageW,colorBarH/pageH])


# Vertical colorbar
colorBarH   =  simH
colorBarHPad=  .0       * cm2inch
colorBarWPad= .25       * cm2inch
colorBarW   = .3        * cm2inch
colorBarB = simB
#colorBarB   = (simB-2*simHPad-2*simH-colorBarHPad-colorBarH)
axColorbar = plt.axes([(pageW-pagePad+colorBarWPad)/pageW,colorBarB/pageH,colorBarW/pageW,colorBarH/pageH])

   












# Extracting data
# =====================
computePostProc = False
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
    for it in range(0,nSteps) :
        outFolder   = "Out_%05d/" % (it)
        State       = Output.readState(rootFolder + simFolder + outFolder + "modelState.json")
    
        dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
        SII = dataSet.data * CharExtra.stress
        subset_SII  = SII[ixCellMin:ixCellMax+1,iyCell]
        
        dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
        P = dataSet.data * CharExtra.stress
        subset_P  = P[ixCellMin:ixCellMax+1,iyCell]
        
        I = np.argmin(subset_P) # find the minimum khi
        SII_t[it] = subset_SII[I] # get the stress corresponding to the minimum value
        
        P_t[it] = subset_P[I] # get the stress corresponding to the minimum value
        Sy_t[it] = C*np.cos(phi) + P_t[it]*np.sin(phi)
        
        time_t[it] = State.time * Setup.Char.time + State.dt
        dt_t[it] = State.dt #* Char.time
        
        I_t[it]     = I
    #endfor
    
    np.savez("./Save/Preambule",SII_t=SII_t,P_t=P_t,Sy_t=Sy_t,time_t=time_t,dt_t=dt_t,I_t=I_t,ixCell=ixCell,iyCell=iyCell)
else:
    loadedData = np.load("./Save/Preambule.npz");
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
P = Setup.Physics.Pback
Sy_back = C*np.cos(phi) + P*np.sin(phi)

# Find interesting values
# =====================
I_sigmaMax = np.argmax(SII_t)
I_sigmaMin = I_sigmaMax + np.argmin(SII_t[I_sigmaMax:])

Delta_timeSoft = time_t[I_sigmaMin] - time_t[I_sigmaMax]


# Units and characteristic values
# =====================
EII = backStrainRate
eta = Setup.MatProps['1'].getRefVisc(0.0,1.0,EII)
G = Setup.MatProps['1'].G
#t = 0
t = eta/G * np.log(2*eta*EII / (2*eta*EII - Sy_back ));
charTime = Sy_back / (2*G*EII * np.exp(-G/eta*t));
timeMaxwell = eta/G
timePlotUnit = charTime#timeMaxwell#Kyr

Pback = Setup.Physics.Pback
sigmaPlotUnit = Pback#MPa



# Plotting 2D
# =====================
StepList = [int(nSteps/3),I_sigmaMax,I_sigmaMin]
for iSim in range(0,nSim):
    outFolder = "Out_%05d/" % (StepList[iSim]-1)
    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'strainRate.bin')
    EII = dataSet.data * 1.0/Char.time
#    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'P.bin')
#    P = dataSet.data * CharExtra.stress
#    dataSet     = Output.getData(rootFolder + simFolder + outFolder + 'sigma_II.bin')
#    SII = dataSet.data * CharExtra.stress
    
    plt.sca(axSimEII["%i" % iSim])
    plt.pcolor(xv - dx/2.0,yv - dy/2.0,np.log10(EII / backStrainRate),vmin=-2, vmax=+2) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
    plt.axis('tight')
    
#    plt.sca(axSimP  ["%i" % iSim])
#    plt.pcolor(xv - dx/2.0,yv - dy/2.0,P/sigmaPlotUnit,vmin=1.5,vmax=3.0) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
#    plt.axis('equal')
#    
#    plt.sca(axSimSII["%i" % iSim])
#    plt.pcolor(xv - dx/2.0,yv - dy/2.0,SII/sigmaPlotUnit,vmin=1.0,vmax=2.0) # -dx/2.0 because pcolor takes the coordinate given as the lower left corner instead of the center
#    plt.plot(xv[I_t[StepList[iSim]],iyCell], yv[I_t[StepList[iSim]],iyCell],'ow',markerFaceColor='none')
#    plt.axis('equal')

    
# Plotting 2D
# =====================
sizeFacBox = 6.0
dictOption = dict(color='w',horizontalAlignment='left',verticalAlignment='baseline',fontName='PT Mono',fontWeight='bold')
xText = xmin+Hbox/sizeFacBox/6.0
yText = ymax-Hbox/sizeFacBox/1.25
for iSim in range(0,nSim):
     plt.sca(axSimEII["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox,xmin+Hbox/sizeFacBox]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+0*nSim],dictOption)
     
     plt.sca(axSimP  ["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox,xmin+Hbox/sizeFacBox]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+1*nSim],dictOption)
     
     plt.sca(axSimSII["%i" % iSim])
     plt.fill([xmin,xmin,xmin+Hbox/sizeFacBox,xmin+Hbox/sizeFacBox]-dx/2,[ymax,ymax-Hbox/sizeFacBox,ymax-Hbox/sizeFacBox,ymax]-dy/2,'k')
     plt.text(xText,yText,letters[iSim+2*nSim],dictOption)

    
    
#    plt.colorbar()

# Colorbar
# =====================
     
## Horizontal colorbar
#plt.sca(axColorbar)
#Cbar = plt.colorbar(cax=axColorbar, orientation='horizontal')
#axColorbar.yaxis.set_label_position('left')
#axColorbar.xaxis.set_ticks_position('top')
#Cbar.set_ticks([-2,-1,0,1,2])
#CbarLims = Cbar.get_clim()
#plt.text(1.02,0.5,"$\dot{\\epsilon}_{II}/\dot{\\epsilon}_{back}$",verticalAlignment='center')
#axColorbar.tick_params(direction='in')
     
# Vertical colorbar
#plt.sca(axColorbar)
Cbar = plt.colorbar(cax=axColorbar, orientation='vertical')
#plt.sca(axSimEII["%i" % iSim])
#Cbar = plt.colorbar()
#axColorbar.yaxis.set_label_position('left')
#axColorbar.xaxis.set_ticks_position('top')
Cbar.set_ticks([-2,-1,0,1,2])
#CbarLims = Cbar.get_clim()
#plt.text(1.02,0.5,"$\dot{\\epsilon}_{II}/\dot{\\epsilon}_{back}$",verticalAlignment='center')
axColorbar.tick_params(direction='in')     





# Plotting Graph
# =====================
#plt.figure(1)
#plt.subplot(2,1,2)
#plt.clf()
plt.sca(axGraph)





plt.plot(time_t/timePlotUnit,(P_t)/sigmaPlotUnit,'-b')
plt.plot(time_t/timePlotUnit,Sy_t/sigmaPlotUnit,'--k',linewidth=4)
plt.plot(time_t/timePlotUnit,SII_t/sigmaPlotUnit,'-r')


plt.plot((0,time_t[-1]/timePlotUnit), np.array((Sy_back,Sy_back))/sigmaPlotUnit,'--k')

#for iSim in range(0,nSim):
#    plt.plot(np.array((1,1))*time_t[StepList[iSim]]/timePlotUnit, np.array((0,2.5)),'--k')

#plt.plot(time_t[I_sigmaMax]/timePlotUnit,SII_t[I_sigmaMax]/sigmaPlotUnit,'or',markerSize=8,markerFaceColor='none')
#plt.plot(time_t[I_sigmaMin]/timePlotUnit,SII_t[I_sigmaMin]/sigmaPlotUnit,'ob',markerSize=8,markerFaceColor='none')


plt.legend(["$\\tau_{y}$","$\\tau_{II}$","$P$","$\\tau_{y}$ at $P_{back}$"])
#plt.title("$P_{back}$ = %.f MPa, C = %.f MPa, G = %.f GPa" % (Setup.Physics.Pback/MPa, Setup.MatProps['1'].cohesion/MPa, Setup.MatProps['1'].G/GPa))
plt.xlabel("$time/time_c$")
plt.ylabel("Stress/$P_{back}$")


print("Delta_timeSoft = %.f Kyr" % (Delta_timeSoft/Kyr))
print("sigmaMax = %.f MPa, sigmaMin = %.f MPa" % (SII_t[I_sigmaMax]/MPa, SII_t[I_sigmaMin]/MPa))
print("sigmaMaxTime = %.f Kyr, sigmaMinTime = %.f Kyr" % (time_t[I_sigmaMax]/Kyr, time_t[I_sigmaMin]/Kyr))
