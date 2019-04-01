# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, '/work/G10501/abauville/Software/StokesFD/src/UserInput')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan, arctan
from PaperDecollement_Utils import getColormap, get_XYandPattern
#from pprint import pprint

#import scipy

import OutputDef as Output

import InputDef as Input
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.image as mpimg
from matplotlib.colors import Normalize
#from scipy import interpolate
#from scipy import ndimage

from Units import *
import time




## Create the folder tree
# ================================
Weak = 60
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/NoTopo/Beta0/Weak%i/" % Weak
# superRootFolder = "/work/G10501/abauville/Paper_Decollement/NoTopo/Beta0/Weak%i/" % Weak
superRootFolder = "/Users/abauville/Output/ListricDecollement/OutputTemp/"

#superDirList = os.listdir(superRootFolder)
#try:
#    superDirList.remove('.DS_Store')
#except ValueError:
#    print("dummy print: no .DS_Store")

#superDirList = ["Lambda83_Hc050_PfW05_GFac005_Beta00"]
superDirList = ["Lambda83_Hc025_PfW05_GFac005_Beta00"]


sampleRate = 5
pointSize = sampleRate/5.0

rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
#if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]
#rootFolder += subFolder + "/"
DirList = os.listdir(rootFolder)
try:
    DirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
try:
    DirList.remove('Input')
except ValueError:
    print("dummy print: no Input")


rootFolders = [''] * len(superDirList) # initialize to the right size
nStepsList = np.zeros(len(superDirList),dtype='int16')
for i in range(len(superDirList)):
    rootFolders[i] = superRootFolder + superDirList[i] + "/Output/"
#    subFoldersList = os.listdir(rootFolders[i])
#    try:
#        subFoldersList.remove('.DS_Store')
#    except ValueError:
#        Dummy=0
#    rootFolders[i] += subFoldersList[0] + "/"
#    
    stepFolderList = os.listdir(rootFolders[i])
    try:
        stepFolderList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
#Setup = Output.readInput(rootFolder +  'Input/input.json')
#s = Setup.Description
#Char = Setup.Char

pushVel = 10.0*cm/yr

Compute = True
if Compute:    
    it=-1
    goldenRatio = (1.0+np.sqrt(5)/2.0)

    # set figure
    
    cm2inch = 0.393701
    figW = 2.0*18.0 * cm2inch -0.01
    figH = figW/goldenRatio -0.005
    plt.close()
    fig = plt.figure(1,figsize=(figW,figH))
    mngr = plt.get_current_fig_manager()
#     to put it into the upper left corner for example:
    FigCoord = mngr.window.geometry().getRect()
    mngr.window.setGeometry(FigCoord[0]-2500,FigCoord[1]-200,FigCoord[2],FigCoord[3])
    Ax1 = plt.axes([0,0,1,1])
        
        
        
                
    # define the limits of the axis
    padAx = 0.1
    ypMin = -padAx
    xpMin = -16.0-padAx
    xpMax = padAx  
    ypMax = (xpMax-xpMin)/goldenRatio-padAx
              
    
    
    # set a font dict
    font = {'family': 'Montserrat',
    'weight': 'bold',
    'size': 20
    }
    
    
    
    nSim = 1#len(superDirList)
#    nSim = 11
    iSim0 = nSim-1
    Hsed = 2.0 * km
    nSteps = nStepsList[iSim0]
    i0 = 0#nSteps-1
    #jump-2
#    iStep = i0
    jump = 1
#    frame = int(i0/jump)
    frame = i0
    for iStep in range(i0,nSteps,jump):

        plt.cla()
        iSub = 0
        for iSim in range(iSim0,nSim):
            tic = time.time()
            iSub+=1
            it+=1
            
            ## Set Folder
            # ================================
            
#            outFolder = os.listdir(rootFolder)[-1]
            outFolder = 'Out_%05d' % (iStep)
            dataFolder = rootFolders[iSim] + outFolder + "/"
            print(rootFolders[iSim])
            
               
            ## Get Data and Pattern
            # ================================
            Setup = Output.readInput(rootFolders[iSim] +  'Input/input.json')
            Char = Setup.Char
            
            

            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=12, xmin=-24.0, xmax=0.0, ymin=0.0,ymax=1.00,nLayersY=5,minStrain=0.0,maxStrain=1.0,mainDir='x')
            
    
            ## Create the colormap many random colors
            # ================================
            if iSub==1:
                # CMAP = getColormap(nColors,"myColorMap",renderer)
                CMAP=([[]])
                CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,.0,.95,.0],
                           RGBShift=[[0.0, 0.0, 0.0], 
                                     [0.0, 0.0, 0.0], 
                                     [-.2,-.2,0.2], 
                                     [0.0, 0.0, 0.0]])
                plt.register_cmap(cmap=CMAP)
                plt.set_cmap("myColorMap")
    
    
            ## Plot
            # ===============================
            outFolder = "../Movies/wWater/Beta00/Weak%02d/" % (Weak) + superDirList[iSim] + "/"
            try:
                os.makedirs(outFolder)       
            except FileExistsError:
                daijoubu = 1
                
            plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
            
#            plt.axis("equal")

            plt.axis([xpMin,xpMax,ypMin,ypMax]) 
            plt.axis("off")
            plt.draw()
            fig.canvas.flush_events()  
            
            # plt.text(xpMin+padAx,1.45,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),fontdict=font)
#            plt.title(superDirList[iSim])
#                plt.pause(0.0001)
#            plt.savefig(outFolder + "Frame_%05d" % frame,dpi=200)
            print("iStep = %i/%i, toc = %.2f" % (iStep, nSteps-1, time.time()-tic))
            frame += 1
                #
#                module_manager = scene.children[0].children[0]
#                module_manager.scalar_lut_manager.show_legend = True
        
