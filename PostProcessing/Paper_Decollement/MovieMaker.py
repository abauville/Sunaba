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





## Create the folder tree
# ================================
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Beta0/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[1]
rootFolder += subFolder + "/"
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
    subFoldersList = os.listdir(rootFolders[i])
    try:
        subFoldersList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    rootFolders[i] += subFoldersList[0] + "/"
    
    stepFolderList = os.listdir(rootFolders[i])
    try:
        stepFolderList.remove('.DS_Store')
    except ValueError:
        Dummy=0
    nStepsList[i] = len(stepFolderList)-1

# Read parameters of this simulation
# =====================
Setup = Output.readInput(rootFolder +  'Input/input.json')
s = Setup.Description
Char = Setup.Char

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
        scene = mlab.figure(1,size=(1980,1080))
#        scene = mlab.figure(1)
        mlab.clf()
        scene.scene.x_plus_view()
        scene.scene.parallel_projection = True
        scene.scene.background = (1.0,1.0,1.0)
        scene.scene.disable_render = False
        
        
                
    # define the limits of the axis
    padAx = 0.1
    ypMin = -padAx
    ypMax = 5.0
    xpMin = -(ypMax-ypMin)*goldenRatio+padAx
    xpMax = padAx            
    
    
    # set a font dict
    font = {'family': 'Montserrat',
    'weight': 'bold',
    'size': 16
    }
    
    
    
#    nSim = len(superDirList)
    nSim = 10
    iSim0 = nSim-1
    Hsed = 2.0 * km
    nSteps = nStepsList[-1]
    i0 = 746#nSteps-1#jump-2
#    iStep = i0
    jump = 1000000
    frame = 0
    for iStep in range(i0,nSteps,jump):
        if renderer == renderMatplotlib:
            plt.cla()
        elif renderer == renderMayavi:
            mlab.clf()
        iSub = 0
        for iSim in range(iSim0,nSim):
            iSub+=1
            it+=1
            
            ## Set Folder
            # ================================
            
#            outFolder = os.listdir(rootFolder)[-1]
            outFolder = 'Out_%05d' % (iStep)
            dataFolder = rootFolders[iSim] + outFolder + "/"
            print(rootFolder)
            print("iStep = %i/%i" % (iStep, nSteps-1))
               
            ## Get Data and Pattern
            # ================================
            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=100)
            
    
            ## Create the colormap many random colors
            # ================================
            if iSub==1:
                CMAP = getColormap(nColors,"myColorMap",renderer)
                if renderer == renderMatplotlib:
                    plt.register_cmap(cmap=CMAP)
                    plt.set_cmap("myColorMap")
    
    
            ## Plot
            # ================================
            
            
            if renderer == renderMatplotlib:
                plt.scatter(PartX,PartY,c=PartPattern,s=15.0,vmin=0.0,vmax=4*nColors-1)      
                
    #            plt.axis("equal")
                
                plt.axis([xpMin,xpMax,ypMin,ypMax]) 
                
                plt.text(xpMin+padAx,1.25,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),fontdict=font)
    #            plt.title(superDirList[iSim])
                plt.pause(0.0001)
#                plt.savefig("/Users/abauville/Output/Paper_Decollement/Movies/Test/Frame_%05d" % frame,dpi=200)
#                
            elif renderer == renderMayavi:
                
                x = np.zeros(PartX.size)
                glyph0 = mlab.points3d(x, PartX, PartY, PartPattern, mask_points=1,vmin=0.0,vmax=4.0*nColors-1.0,scale_mode="none",mode="cone",resolution=8,scale_factor=0.1)
                glyph0.actor.property.lighting = False
                glyph0.glyph.glyph_source.glyph_source.capping = False
                glyph0.module_manager.scalar_lut_manager.lut.table = CMAP
#                
                text = mlab.text(0.025,0.25,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),color=(0.0,0.0,0.0),width=0.15)
#                text = engine.scenes[0].children[0].children[0].children[1]
#                text.property.shadow_offset = array([ 1, -1])
                text.property.bold = True
#                
                scene.scene.camera.position = [39.5,-5.4, 2.8]
                scene.scene.camera.focal_point = [0.0, scene.scene.camera.position[1], scene.scene.camera.position[2]]
                scene.scene.camera.view_angle = 30.0
                scene.scene.camera.view_up = [0.0, 0.0, 1.0]
                scene.scene.camera.clipping_range = [38.09508309975636, 41.38114149447358]
                scene.scene.camera.compute_view_plane_normal()
                scene.scene.camera.zoom(3.5)
#                scene.scene.render()
                mlab.draw()
                mlab.show()
#                scene.scene.save(u'/Users/abauville/Output/Paper_Decollement/Movies/Test/Frame_%05d.png' % frame)


            frame += 1
                #
#                module_manager = scene.children[0].children[0]
#                module_manager.scalar_lut_manager.show_legend = True
        
