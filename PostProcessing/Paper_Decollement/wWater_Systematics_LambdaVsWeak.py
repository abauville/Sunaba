#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Output reading test

import sys
import os
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
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
import CritTaper_Style
#from scipy import interpolate
#from scipy import ndimage
import CritTaper_dataMaker
from Units import *
from numpy import array as arr

#weakList = arr([1, 5, 10, 20, 30, 40, 50, 60, 70, 80])
weakList = arr([1, 10, 20, 30, 40, 50, 60, 70, 80])
#weakList = arr([40, 40 ,40])
nW = len(weakList)
LambdaList = arr([00,40,60,80])
#LambdaList = arr([60])
nL = len(LambdaList)

## Get Data from CritTaper theory analalysis
Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
#CMAP_ref = plt.get_cmap(Style.colormap)._lut
CMAP_ref = arr([[0.0,0.0,0.0]])

#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_Ref_all - alphas_WB_up_all 

IL = np.zeros(nL,dtype=np.int)
IW = np.zeros(nW,dtype=np.int)

i = 0
for Lambda in LambdaList/100.0:
    IL[i] = np.argmin(abs(LambdaRef_list-Lambda))
    i+=1

i = 0
for weak in weakList/100.0:
    IW[i] = np.argmin(abs(chi_list-weak))
    i+=1
    
    
beta = 0.0
alphas_diff_bar  = np.zeros([nL,nW])
IRefColorMap_Big = np.zeros([nL,nW],dtype=np.int)
IRefColorMap = np.zeros(nL*nW,dtype=np.int)
#taper_angles = np.zeros([nL,nW])


# find colors in the colormap
cmin = -1.0
cmax = 1.0

nC = CMAP_ref.shape[0]
refCmapValues = np.linspace(cmin,cmax,nC)

iiL = 0
for iL in IL:
    iiW = 0
    for iW in IW:
        iB = np.argmin(np.abs(betas_all[iL,iW,:]-beta))
#        alphas_diff[iL,iW] = alphas_diff_all[iL,iW,iB]
#        alphas_width[iL,iW] = alphas_WB_up_all[iL,iW,iB] - alphas_WB_low_all[iL,iW,iB]
#        alphas_WB_up[iL,iW] = alphas_WB_up_all[iL,iW,iB]
#        taper_angles[iL,iW] = betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB]
        alphas_diff_bar[iiL,iiW] = alphas_diff_all[iL,iW,iB]/(betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB])

#        if (Lambdas[iL,iW]+chis[iL,iW]<1.0):
#            alphas_diff_bar[iiL,iiW] = (Lambdas_Ref_all[iL,iW,iB]-1.0+chis_all[iL,iW,iB])
#        else:
#            alphas_diff_bar[iiL,iiW] = Lambdas_Ref_all[iL,iW,iB]*chis_all[iL,iW,iB]
        IRefColorMap_Big[iiL,iiW] = np.argmin(np.abs(refCmapValues-alphas_diff_bar[iiL,iiW]))
        iiW+=1
    iiL+=1




## Create the folder tree
# ================================
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Static2/Beta0/"
Weak = 50
wWater = ""#"_w_Water"
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/NoTopo/Beta0/Weak%i%s/" % (Weak, wWater)
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta00/"
superDirList = os.listdir(superRootFolder)
try:
    superDirList.remove('.DS_Store')
except ValueError:
    print("dummy print: no .DS_Store")
    
    
ProductionMode = True


if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 10
    pointSize = sampleRate/100.0
else:
    sampleRate = 2000
    pointSize = 1.0

superDirList = []
i = 0
iW = 0
for weak in weakList:
    iL = 0
    for Lambda in LambdaList:
        superDirList.append("Weak%02d/Lambda%02d" % (weak, Lambda))
        IRefColorMap[i] = IRefColorMap_Big[iL,iW]
        i+=1
        iL+=1
    iW+=1




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
        nrows= nW
        ncols = nL
        # set figure
        cm2inch = 0.393701
#        figW = 2.0*18.0 * cm2inch
#        figH = 2.0*figW/goldenRatio
        
#        figW = 1.0*21.0 * cm2inch
#        figH = figW/goldenRatio
#        figW = 29.7 * cm2inch
#        figH = 21.0*cm2inch#figW/goldenRatio
        figW = 29.7 * cm2inch 
        figH = 21.0 * cm2inch#21.0*cm2inch#figW/goldenRatio
        
#        plt.close("all")
        fig = plt.figure(1,figsize=(figW,figH))
        plt.clf()
        fig.set_dpi(220)
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
        ypMax = 5.0
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
    'size': 22
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
            
            PartX = []
            PartY = []
            PartPattern = []
            
            ax = plt.sca(Ax[iSub-1])
            
            PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
            plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
            
            
#            
#            dataSet  = Output.getData(dataFolder + 'strain.bin',True)
#            X = np.linspace(dataSet.xmin,dataSet.xmax,dataSet.nx)
#            Y = np.linspace(dataSet.ymin,dataSet.ymax,dataSet.ny)
#            X, Y = np.meshgrid(X,Y)
#            strain = dataSet.data
#            rate = 1
#            plt.pcolor(X.T[::rate,::rate],Y.T[::rate,::rate],strain[::rate,::rate],vmin=0.1,vmax=5.0)

            
#            Index = [0, 0, 0, 0,
#                     0, 0, 0, 0,
#                     0, 0, 0, 1,
#                     0, 0, 1, 1,
#                     0, 1, 1, 2,
#                     1, 1, 1, 2,
#                     1, 1, 2, 2]
#
##            CMAP=np.array([CMAP_ref[IRefColorMap[iSim],:]])
##            CMAP_all = arr([[219, 59, 38,255],
##                            [239,189, 64,255],
##                            [ 70,159,248,255]])/255.0
#            CMAP = arr([CMAP_all[Index[iSim],:]])
#            CMAP = arr([.95,.85,.5,1.0])
            CMAP_all = arr([[.8,1.,.85,1.0],
                            [1.,.8,.85,1.0],
                            [.82,.8,.95,1.0]])
            CMAP = arr([CMAP_all[2]])
#            CMAP = arr([CMAP_all[Index[iSim],:]])
            CMAP = getColormap(nColors,"myColorMap",renderer,CMAP=CMAP,shiftHLayerColors=False)
            plt.register_cmap(cmap=CMAP)
            plt.set_cmap("myColorMap")
            
            
            
#                CMAP = getColormap(nColors,"myColorMap",renderer,CMAP=Color)
#                if renderer == renderMatplotlib:
#                    plt.register_cmap(cmap=CMAP)
#                    plt.set_cmap("myColorMap")
#            else:
#                break
    
            ## Plot
            # ================================
            
            
#                ax = plt.subplot(nrows,ncols,iSub,anchor="SE")
            
            
            #plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1)      
            
#                break
#                plt.axis("equal")
            
#                plt.axis("scaled",anchor="SE")
            plt.axis([xpMin,xpMax,ypMin,ypMax]) 
#                plt.ylim(ypMin,ypMax)
#                plt.axis([xpMin,xpMax,ypMin,ypMax]) 
            
            
            Sediment = Output.readInput(rootFolders[iSim] +  'Input/input.json').MatProps['1']
            coh = Sediment.cohesion
            phi = Sediment.frictionAngle
            g = 9.81
            rho = Sediment.rho0
            
            Hc2 = coh/(rho*g*tan(phi))      / Hsed
            Hc3 = coh/(rho*g)               /Hsed

#                plt.text(xpMin+padAx,2.25,"SHORT. = %02.1f, Hc2 = %.3f, c = %.1fMPa" % (timeSim*pushVel/Hsed, Hc2, coh/1e6),fontdict=font)
            
            if ((iSub-1)%ncols==0 and (iSub-1)<ncols): #upper left corner
                plt.text(xpMin-4.0*padAx,ypMax+25.0*padAx,"$\mathbf{\\chi}$ $[\%]$",fontdict=font,horizontalAlignment='right',verticalAlignment='top',rotation=90)
                plt.text(xpMin+0.0*padAx,ypMax-5.0*padAx,"$\mathbf{\\lambda} \; [\%]$",fontdict=font,horizontalAlignment='left',verticalAlignment='baseline')
                
            if ((iSub-1)%ncols==0):
                weak = float(superDirList[iSim][4:6])
#                    if (Hc<1.0):
#                        txtString = "1/%i" % (round(1.0/Hc))
#                    else:
#                        txtString = "%i" % (round(Hc))
                txtString = "%i" % (weak)
                plt.text(xpMin-4.0*padAx,0.0,txtString,fontdict=font,horizontalAlignment='right',verticalAlignment='baseline')
            if ((iSub-1)<ncols):
                if superDirList[iSim][-2] == "a":
                    lambdaFac = float(superDirList[iSim][-1:]) 
                else:
                    lambdaFac = float(superDirList[iSim][-2:]) 
                txtString = "%i" % lambdaFac
                plt.text((xpMin+xpMax)/2.0,ypMax-5.0*padAx,txtString,fontdict=font,horizontalAlignment='center',verticalAlignment='baseline')
            
            plt.axis("off")
#               
            frame += 1
                #
#                module_manager = scene.children[0].children[0]
#                module_manager.scalar_lut_manager.show_legend = True
        
plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/Systematics_LambdaVsWeak3",dpi=500)