# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
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
Weak = 50
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/NoTopo/Beta0/Weak%i/" % Weak
#superDirList = os.listdir(superRootFolder)
#try:
#    superDirList.remove('.DS_Store')
#except ValueError:
#    print("dummy print: no .DS_Store")

superDirList = ["Hc1.000_Lambda90"]

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
nStepsList = np.zeros(len(superDirList),dtype='int16');
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
    renderMatplotlib = 0
    renderMayavi = 1
    renderer = 0 # 0 Matplotlib; 1 Mayavi
    if renderer == renderMatplotlib:
        # set figure
        
        cm2inch = 0.393701
        figW = 2.0*18.0 * cm2inch
        figH = figW/goldenRatio
        plt.close()
        plt.figure(1,figsize=(figW,figH))
        mngr = plt.get_current_fig_manager()
    #     to put it into the upper left corner for example:
        FigCoord = mngr.window.geometry().getRect()
        mngr.window.setGeometry(FigCoord[0]-2500,FigCoord[1]-200,FigCoord[2],FigCoord[3])
        Ax1 = plt.axes([0,0,1,1])
        
    elif renderer == renderMayavi:
        import mayavi.mlab as mlab
        try:
            mlab.close()
        except:
            Dummy = 0    
        scene = mlab.figure(1,size=(1980,1080+35)) # +35 to take into account the GUI...
#        scene = mlab.figure(1)
        mlab.clf()
        scene.scene.x_plus_view()
        scene.scene.parallel_projection = True
        scene.scene.background = (1.0,1.0,1.0)
        scene.scene.disable_render = False
        
        
                
    # define the limits of the axis
    padAx = 0.1
    ypMin = -padAx
    xpMin = -16.0
    xpMax = padAx  
    ypMax = (xpMax-xpMin)/goldenRatio-padAx
              
    
    
    # set a font dict
    font = {'family': 'Montserrat',
    'weight': 'bold',
    'size': 16
    }
    
    
    
    nSim = 1#len(superDirList)
#    nSim = 11
    iSim0 = nSim-1
    Hsed = 2.0 * km
    nSteps = nStepsList[iSim0]
    i0 = 141#nSteps-1
    #jump-2
#    iStep = i0
    jump = 1
#    frame = int(i0/jump)
    frame = 0
    for iStep in range(i0,nSteps,jump):
        if renderer == renderMatplotlib:
            plt.cla()
        elif renderer == renderMayavi:
            mlab.clf()
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
            print("iStep = %i/%i" % (iStep, nSteps-1))
               
            ## Get Data and Pattern
            # ================================
            Setup = Output.readInput(rootFolders[iSim] +  'Input/input.json')
            Char = Setup.Char
            
            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=1)
            
    
            ## Create the colormap many random colors
            # ================================
            if iSub==1:
                CMAP = getColormap(nColors,"myColorMap",renderer)
                if renderer == renderMatplotlib:
                    plt.register_cmap(cmap=CMAP)
                    plt.set_cmap("myColorMap")
    
    
            ## Plot
            # ===============================
            outFolder = "/Users/abauville/Output/Paper_Decollement/Movies/NoTopo/Beta0/Weak%i/" % (Weak) + superDirList[iSim] + "/"
            try:
                os.makedirs(outFolder)       
            except FileExistsError:
                daijoubu = 1
            if renderer == renderMatplotlib:
                plt.scatter(PartX,PartY,c=PartPattern,s=0.5,vmin=0.0,vmax=4*nColors-1)      
                
    #            plt.axis("equal")
                
                plt.axis([xpMin,xpMax,ypMin,ypMax]) 
                
                plt.text(xpMin+padAx,1.40,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),fontdict=font)
    #            plt.title(superDirList[iSim])
#                plt.pause(0.0001)
                plt.savefig(outFolder + "Frame_%05d" % frame,dpi=200)
#                
            elif renderer == renderMayavi:
                
                x = np.zeros(PartX.size)
                glyph0 = mlab.points3d(x, PartX, PartY, PartPattern, mask_points=1,vmin=0.0,vmax=4.0*nColors-1.0,scale_mode="none",mode="cone",resolution=4,scale_factor=0.1)
                glyph0.actor.property.lighting = False
                glyph0.glyph.glyph_source.glyph_source.capping = False
                glyph0.module_manager.scalar_lut_manager.lut.table = CMAP
#                
                text = mlab.text(0.025,0.35,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),color=(0.0,0.0,0.0),width=0.08)
#                text = engine.scenes[0].children[0].children[0].children[1]
#                text.property.shadow_offset = array([ 1, -1])
                text.property.bold = True

                scene.scene.camera.position = [39.5,-6.6, 3.5]
                scene.scene.camera.focal_point = [0.0, scene.scene.camera.position[1], scene.scene.camera.position[2]]
                scene.scene.camera.view_angle = 30.0
                scene.scene.camera.view_up = [0.0, 0.0, 1.0]
                scene.scene.camera.clipping_range = [38.09508309975636, 41.38114149447358]
                scene.scene.camera.compute_view_plane_normal()
                scene.scene.camera.zoom(2.5)
                
#                scene.scene.render()
#                mlab.draw()
#                mlab.show()
#                scene.scene.save(u'/Users/abauville/Output/Paper_Decollement/Movies/Static/Test00/Frame_%05d.png' % frame)
                scene.scene.save(outFolder + 'Frame_%05d.png' % frame)

                toc = time.time()
                toc -= tic
                print("step: %.1f s" % toc)
            frame += 1
                #
#                module_manager = scene.children[0].children[0]
#                module_manager.scalar_lut_manager.show_legend = True
        
