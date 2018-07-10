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
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Static2/Beta0/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    
## Replace superDirFolder with a custom subset
# ================================
#superDirList = ['C1.0_Weak10_Lambda40',
#                'C1.0_Weak20_Lambda40',
#                'C5.0_Weak10_Lambda40',
#                'C5.0_Weak20_Lambda40',
#                'C1.0_Weak20_Lambda40_G20',
#                'C1.0_Weak20_Lambda40_G20_swIni0.5_swEnd_1.5']

#superDirList = ['C1.0_Weak10_Lambda80',
#                'C1.0_Weak20_Lambda80',
#                'C5.0_Weak10_Lambda80',
#                'C5.0_Weak20_Lambda80',
#                'C5.0_Weak10_Lambda80_G20_swIni0.1_swEnd_1.0',
#                'C5.0_Weak10_Lambda80_G20_swIni0.5_swEnd_1.5',
#                'C5.0_Weak10_Lambda80_G5_swIni0.5_swEnd_1.5',
#                'C5.0_Weak20_Lambda80_G5_swIni0.5_swEnd_1.5',
#                'C5.0_Weak20_Lambda80_timeFac3']
#    
#    
    
#superDirList = ['Hc0.062_Weak20_Lambda0',
#                'Hc0.062_Weak20_Lambda60',
#                'Hc0.062_Weak20_Lambda90',
#                'Hc0.500_Weak20_Lambda0',
#                'Hc0.500_Weak20_Lambda60',
#                'Hc0.500_Weak20_Lambda90', 
#                'Hc2.000_Weak20_Lambda0',                
#                'Hc2.000_Weak20_Lambda60',                                               
#                'Hc2.000_Weak20_Lambda90']

superDirList = ['Hc0.062_Weak10_Lambda0',
                'Hc0.062_Weak10_Lambda60',
                'Hc0.062_Weak10_Lambda90',
                'Hc0.500_Weak10_Lambda0',
                'Hc0.500_Weak10_Lambda60',
                'Hc0.500_Weak10_Lambda90', 
                'Hc2.000_Weak10_Lambda0',                
                'Hc2.000_Weak10_Lambda60',                                               
                'Hc2.000_Weak10_Lambda90']

    
rootFolder = superRootFolder + superDirList[0] + "/Output/"
subFolder = os.listdir(rootFolder)[0]
if subFolder == ".DS_Store": subFolder = os.listdir(rootFolder)[1]

#
#rootFolder += subFolder + "/"
#DirList = os.listdir(rootFolder)
#try:
#    DirList.remove('.DS_Store')
#except ValueError:
#    print("dummy print: no .DS_Store")
#try:
#    DirList.remove('Input')
#except ValueError:
#    print("dummy print: no Input")


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
    goldenRatio = (1.0+np.sqrt(5))/2.0
    renderMatplotlib = 0
    renderMayavi = 1
    renderer = 0 # 0 Matplotlib; 1 Mayavi
    if renderer == renderMatplotlib:
        nrows= 3
        ncols = 3
        # set figure
        cm2inch = 0.393701
#        figW = 2.0*18.0 * cm2inch
#        figH = 2.0*figW/goldenRatio
        
        figW = 2.0*25.0 * cm2inch
        figH = figW/goldenRatio
        
        plt.close("all")
        plt.figure(1,figsize=(figW,figH))
#        plt.subplots(nrows=nrows, ncols=ncols, figsize=(figW,figH))
        mngr = plt.get_current_fig_manager()
    #     to put it into the upper left corner for example:
        FigCoord = mngr.window.geometry().getRect()
        mngr.window.setGeometry(FigCoord[0]-2500,FigCoord[1]-200,FigCoord[2],FigCoord[3])
#        Ax1 = plt.axes([0.1,0.1,0.8,0.8])
        
#        plt.tight_layout()
        
#        axPos = plt.gca().get_position().extents
#        axPos = plt.gca().get_tightbbox(plt.gcf().canvas.renderer).extents
#        aspectRatio = (axPos[2]-axPos[0])/(axPos[3]-axPos[1])
        aspectRatio = 3.5

    
#        axPos = plt.subplot(nrows,2,1).get_position().extents
#        aspectRatio = (axPos[2]-axPos[0])/(axPos[3]-axPos[1])
        # define the limits of the axis
        padAx = 0.1
        ypMin = -padAx
        ypMax = 4.5
    #    xpMin = -(ypMax-ypMin)*goldenRatio+padAx
        xpMin = -(ypMax-ypMin)*aspectRatio+padAx
        xpMax = padAx    
        
        
        Ax = [];#dict()
        LeftPad = 0.06
        RightPad = 0.04
        TopPad = 0.05
        BottomPad = 0.05
        xPad        = 0.0125
        yPad        = 0.0125
        TotalWidth  = 1.0-LeftPad-RightPad
        TotalHeight = 1.0-TopPad-BottomPad
        Width       = (TotalWidth  - (ncols-1)*xPad)/ncols
        Height      = Width * (ypMax-ypMin)/(xpMax-xpMin)
        for i in range(nrows):
            for j in range(ncols):
                Ax.append(plt.axes([LeftPad+xPad*(j)+Width*(j),1.0-(TopPad+yPad*(i)+Height*(i+1))*figW/figH,Width,Height*figW/figH]))
        
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
        
        
                
            
    

    
    # set a font dict
    font = {'family': 'Montserrat',
    'weight': 'bold',
    'size': 16
    }
    
    
    
    nSim = len(superDirList)
#    nSim = 7
    iSim0 = 0#nSim-1
    Hsed = 2.0 * km
    nSteps = 747#nStepsList[-1]
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
            
            outFolder = os.listdir(rootFolders[iSim])[-1]
#            if (iSim==6):
#                outFolder = os.listdir(rootFolders[iSim-1])[-1]
#            outFolder = 'Out_%05d' % (iStep)
            dataFolder = rootFolders[iSim] + outFolder + "/"
            print(rootFolders[iSim])
            #print("iStep = %i/%i" % (iStep, nSteps-1))
            print(outFolder)  
            ## Get Data and Pattern
            # ================================
            Char = Output.readInput(rootFolders[iSim] +  'Input/input.json').Char
            timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=25, nLayersX=1, nLayersY=0.00)
            
    
            ## Create the colormap many random colors
            # ================================
            if iSub==1:
                CMAP = getColormap(nColors,"myColorMap",renderer)
                if renderer == renderMatplotlib:
                    plt.register_cmap(cmap=CMAP)
                    plt.set_cmap("myColorMap")
#            else:
#                break
    
            ## Plot
            # ================================
            
            
            if renderer == renderMatplotlib:
#                ax = plt.subplot(nrows,ncols,iSub,anchor="SE")
                ax = plt.sca(Ax[iSub-1])
                
                plt.scatter(PartX,PartY,c=PartPattern,s=0.5,vmin=0.0,vmax=4*nColors-1)      
                
                
#                plt.axis("equal")
                
#                plt.axis("scaled",anchor="SE")
                plt.axis([xpMin,xpMax,ypMin,ypMax]) 
#                plt.ylim(ypMin,ypMax)
#                plt.axis([xpMin,xpMax,ypMin,ypMax]) 
                
                plt.text(xpMin+padAx,2.25,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),fontdict=font)
                
                if ((iSub-1)%ncols==0 and (iSub-1)<ncols): #upper left corner
                    plt.text(xpMin-4.0*padAx,ypMax-5.0*padAx,"Hc",fontdict=font,horizontalAlignment='right',verticalAlignment='top')
                    plt.text(xpMin+0.0*padAx,ypMax-5.0*padAx,"$\mathbf{\\lambda}$",fontdict=font,horizontalAlignment='left',verticalAlignment='baseline')
                    
                if ((iSub-1)%ncols==0):
#                    plt.text(xpMin+padAx,2.25,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),fontdict=font)
                    Hc = float(superDirList[iSim][2:7])
                    if (Hc<1.0):
                        txtString = "1/%i" % (round(1.0/Hc))
                    else:
                        txtString = "%i" % (round(Hc))
                    plt.text(xpMin-4.0*padAx,0.0,txtString,fontdict=font,horizontalAlignment='right',verticalAlignment='baseline')
                if ((iSub-1)<ncols):
                    if superDirList[iSim][-2] == "a":
                        lambdaFac = float(superDirList[iSim][-1:]) / 100.0
                    else:
                        lambdaFac = float(superDirList[iSim][-2:]) / 100.0
                    txtString = "%.1f" % lambdaFac
                    plt.text((xpMin+xpMax)/2.0,ypMax-5.0*padAx,txtString,fontdict=font,horizontalAlignment='center',verticalAlignment='baseline')
                
                plt.axis("off")
                plt.pause(0.0001)
#                plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/Systematics_Beta0_Weak20",dpi=600)
#                
            elif renderer == renderMayavi:
                
                x = np.zeros(PartX.size)
                glyph0 = mlab.points3d(x, PartX, PartY, PartPattern, mask_points=1,vmin=0.0,vmax=4.0*nColors-1.0,scale_mode="none",mode="cone",resolution=8,scale_factor=0.1)
                glyph0.actor.property.lighting = False
                glyph0.glyph.glyph_source.glyph_source.capping = False
                glyph0.module_manager.scalar_lut_manager.lut.table = CMAP
#                
                text = mlab.text(0.025,0.25,"SHORT. = %02.1f" % (timeSim*pushVel/Hsed),color=(0.0,0.0,0.0),width=0.25)
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
        
