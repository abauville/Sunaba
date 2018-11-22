#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:47 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '..//CriticalTaper')
sys.path.insert(0, '../')
import OutputDef as Output
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern
from CritTaper_utils import Taper

# Misc
# =========================================
Style = CritTaper_Style.Style()

#   Define chi_list
# =========================================
#thisFile_chi_list = [1, 20, 60, 10, 40, 80]
thisFile_chi_list = [1, 20, 60]

nC = len(thisFile_chi_list)
nSim = nC


#tSteps_list_Fig02 = arr([[214, 428, 641, 855, 1069, 1282, 1496, 1710],
#                         [248, 495, 743, 990, 1238, 1486, 1733, 1981],
#                         [134, 269, 403, 538,  672,  806,  941, 1075],
##                         [  0, 250, 500, 750, 1000, 1250, 1500, 2000], # chi 05
#                         [241, 482, 722, 963, 1204, 1444, 1685, 1926], # chi 10
#                         [176, 352, 528, 704,  880, 1056, 1232, 2014],
#                         [176, 352, 528, 704,  880, 1056, 1232, 696]]) 

tSteps_list_Fig02 = arr([[216, 432, 648, 864, 1080, 1296, 1512, 1728],
                     [249, 497, 746, 994, 1243, 1492, 1740, 1989],
                     [134, 269, 403, 538,  672,  806,  941, 1075]]) 


#  Get colors of Type
# =========================================
beta= 0.0 * np.pi/180.0
Lambda = 0.6
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0))
Lambdas = loadedData["Lambdas"][()]
chis    = loadedData["chis"][()]
floatType = loadedData["floatType"][()]
alphas_diff = loadedData["alphas_diff"][()]
taper_angles = loadedData["taper_angles"][()]

CMAP, colorList_Type = Style.getCmap_Type()
Type_list = np.linspace(-1.0,2.0,colorList_Type.shape[0])

ITs = np.zeros(nSim,np.int)
Types = np.zeros(nSim)
for iSim in range(nSim):
    Lambda = 0.6#Lambdas[iSim]
    chi = thisFile_chi_list[iSim]/100.0
    IL = np.argmin(np.abs(Lambdas-Lambda))
    IC = np.argmin(np.abs(chis-chi))
    
    Type = floatType[IL,IC]
    Types[iSim] = Type
    ITs[iSim] = np.argmin(np.abs(Type_list-Type))
#        print(IC)
    
    







#  Figure
# =========================================
#aspectRatio = 1.0/5.0
fig             = Figz_Utils.Figure(103,height=29.7,width=21.0,mode='production')
#fig             = Figz_Utils.Figure(303,height=29.7,width=21.0,mode='draft')
#fig             = Figz_Utils.Figure(103,height=29.7)

bottomMarginPad = 2.0


aspectRatioDrawing = .4
aspectRatioAlpha = 1.0
aspectRatioFault = 1.5

yShiftTop = 0.5
#yShift0 = AxesAlpha['info']['plotsHeight']+AxesAlpha['info']['yPad']
yShift = yShiftTop
xPad = 0.25
yPad = 10.0



# =====

yShift0 = 0.0
yShift += yShift0
nCol = 3
nRow = int(np.ceil(nSim/nCol))

nRow = 1

AxesDrawing      = Figz_Utils.makeAxes(fig,nRow,nCol,aspectRatio=aspectRatioDrawing,
                                      leftMarginPad=0.0,rightMarginPad=0.0,
                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                      xPad=xPad,yPad=yPad,
                                      setAspectRatioBasedOn='x')


yShift0 = AxesDrawing['info']['plotsHeight']+0.05#AxesAlpha['info']['yPad']
yShift += yShift0
AxesxFault      = Figz_Utils.makeAxes(fig,nRow,nCol,aspectRatio=aspectRatioFault,
                                      leftMarginPad=0.0,rightMarginPad=0.0,
                                      topMarginPad=yShift,bottomMarginPad = bottomMarginPad,
                                      xPad=xPad,yPad=3.78,
                                      setAspectRatioBasedOn='x')



#   File system
# =========================================
beta = 0.0
#superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater/Beta%02d/" % round(beta*180.0/np.pi*10.0)
superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output/wWater_Select2/Beta00/"
superDirList = []
i = 0
iW = 0
for iSim in range(nSim):
    superDirList.append("Weak%02d/Lambda%02d" % (thisFile_chi_list[iSim], 60))

#  Production mode
# =========================================
ProductionMode = True
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/30.0
else:
    sampleRate = 10
    pointSize = sampleRate/100.0


    

#  Colors xFaults
# =========================================
colors_xFault=arr([[0.2,0.2,0.2,1.0],
                   [.9,0.2,.9,1.0],
                   [.6,0.3,.2,1.0]])

colorWedge = arr([.85,.85,.9,1.0])
    
#  Plot wedge drawings
# =========================================

for iSim in range(nSim):

    dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + "Out_%05d/" % tSteps_list_Fig02[iSim,-1]
    Char = Output.readInput(superRootFolder + superDirList[iSim] + '/Input/input.json').Char
    timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
    
    PartX = []
    PartY = []
    PartPattern = []
    
    iRow = int(np.ceil((iSim+1)/nCol))
    iCol = iSim%nCol + 1

    ax = plt.sca(AxesDrawing['%i%i' % (iRow,iCol)])
    
    PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,minStrain=1.0,maxStrain=5.0)
#    PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=1, xmin=-4.0, xmax=0.0, ymin=0.0,ymax=1.15,nLayersY=6,minStrain=1.0,maxStrain=5.0,mainDir='y')
    plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')      
    
    
    xmin = -15.25
    xmax = 0.0
    
    
   
    
    
    plt.axis([xmin,0.0,0.0,-1.0*aspectRatioDrawing*xmin])
#    plt.axis("off")

    CMAP = colorWedge
#    CMAP[:-1] = colorList_Type[ITs[iSim],:]
    
    

    RGBShift = [[0.0,0.0,0.0],
                [0.0,0.0,0.0],
                [0.0,0.0,0.0],
                [0.0,0.0,0.0]]
    
    CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,darknessFactor=[1.0,0.0,1.0,0.0],RGBShift=RGBShift)
   

#    CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False)    
    
    plt.register_cmap(cmap=CMAP)
    plt.set_cmap("myColorMap")

    plt.axis('off')
    





## Figure xFault
# ============================================

Setup = Output.readInput(superRootFolder + superDirList[0] +  '/Input/input.json')

yr = 365.25*24.0*3600.0
kyr = 1000.0 * yr
beta = 0.0
Lambda = 0.6
DataFolder = "/Users/abauville/Output/Paper_Decollement/Figz/Data/"


loadedData_a = np.load(DataFolder + "locSlopes_Beta%02d_Lambda%02d_WeakDense_winSize64.npy" % (beta*10.0*180.0/np.pi,Lambda*100.0)).item(0)    
#loadedData_b = np.load(DataFolder + "locSlopes_Beta%02d_all.npy" % (beta*10.0*180.0/np.pi)).item(0)



alphas_Ref = np.zeros(nSim)
alphas_WF = np.zeros(nSim)
alphas_WB_up = np.zeros(nSim)
alphas_WB_low = np.zeros(nSim)
iSim0 = 0

transparency = 0.2
Color  = arr([[.85,.15,.25],[.25,.5,.5],[.85,.15,.25]])
Color_w_transparency = arr([[.25,.5,.5,transparency],
                            [.85,.15,.25,transparency],
                            [.85,.15,.25,transparency]])


    
    
    

plot=0
for iSim in range(iSim0,nSim):  
    print("iSim = %i/%i" % (iSim,nSim))

    Lambda = 0.6#Lambdas[iSim]
    chi = thisFile_chi_list[iSim]/100.0
    iRow = int(np.ceil((iSim+1)/nCol))
    iCol = iSim%nCol + 1
    
    

    
    
    
    
#    if iSim == 3:
    thisData = loadedData_a["Lambda%02d_chi%02d" % (Lambda*100,chi*100)]
#    else:
#        thisData = loadedData_b["Lambda%02d_chi%02d" % (Lambda*100,chi*100)]
        
    locSlopes = thisData["locSlopes"]
    lenPrisms = thisData["lenPrisms"]
    time_list = thisData["time_list"]
    tSteps = thisData["tSteps"]
    xFront = thisData["xFront"]
    xBase = thisData["xBase"]
    xMid = thisData["xMid"]
    
    x0 = abs(xmin)
    x1 = 0.0
#    x0 = time_list[0]/kyr
#    x1 = time_list[-1]/kyr
    lastStep = int(tSteps_list_Fig02[iSim,-1]+1)
    y0 = 0
    y1 = time_list[lastStep-1]/kyr
    
    nSteps = len(tSteps)
    


    plt.sca(AxesxFault['%i%i' % (iRow,iCol)])    
    ax = plt.gca()
    ax.patch.set_facecolor([0.0,0.0,0.0,0.0])
    if iSim == 0:
        plt.text(x0+(x1-x0)*0.5,y0+(y1-y0)*0.015,'x position')
        plt.text(x0+(x1-x0)*0.97,y0+(y1-y0)*0.95,'time',horizontalAlignment='center',rotation=90)
    Letters='ABCDEF'
    if Types[iSim]>-0.05:
#        plt.fill(x0+(x1-x0)*arr([.0,.43,.43,.0]),y0+(y1-y0)*arr([.0,.0,.06,.06]),color=colorList_Type[ITs[iSim],:])
        plt.text(x0+(x1-x0)*0.01,y0+(y1-y0)*0.015,Letters[iSim] + '. Type=%i' % abs(Types[iSim]),weight='bold')
    else:
        plt.fill(x0+(x1-x0)*arr([.0,.445,.445,.0]),y0+(y1-y0)*arr([.0,.0,.06,.06]),color=colorList_Type[ITs[iSim],:])
        plt.text(x0+(x1-x0)*0.01,y0+(y1-y0)*0.015,Letters[iSim] + '. Type=%.1f' % Types[iSim],weight='bold')
    plt.xticks([])
    plt.yticks([])
    
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    lc = 2.0e3
    dx = (Setup.Grid.xmax-Setup.Grid.xmin)/(Setup.Grid.nxC-1)/lc
    

#    color = colorList_Type[ITs[iSim],:]
#    color = [.85,.85,.88]
    if iSim == 4:
        colorWedge2 = colorWedge*1.12
        colorWedge2[colorWedge2>1.0] = 1.0
        plt.fill((Setup.Grid.nxC-np.concatenate([xFront[:lastStep],[Setup.Grid.nxC,Setup.Grid.nxC]]))*dx,np.concatenate([time_list[:lastStep],[time_list[lastStep-1],0]])/kyr,color=colorWedge2,linewidth=1.0)
        plt.fill((Setup.Grid.nxC-np.concatenate([xFront[:1409],[Setup.Grid.nxC,Setup.Grid.nxC]]))*dx,np.concatenate([time_list[:1409],[time_list[1408],0]])/kyr,color=colorWedge,linewidth=1.0)
    else:     
        plt.fill((Setup.Grid.nxC-np.concatenate([xFront[:lastStep],[Setup.Grid.nxC,Setup.Grid.nxC]]))*dx,np.concatenate([time_list[:lastStep],[time_list[lastStep-1],0]])/kyr,color=colorWedge,linewidth=1.0)
        

    plt.plot((Setup.Grid.nxC-xMid  [:lastStep])*dx,time_list[:lastStep]/kyr,'-',color=colors_xFault[1,:-1],linewidth=1.5,markersize=.5)
    plt.plot((Setup.Grid.nxC-xBase [:lastStep])*dx,time_list[:lastStep]/kyr,'-',color=colors_xFault[2,:-1],linewidth=1.5,markersize=.5)


    if iSim < 3:
        for iStep in range(len(tSteps_list_Fig02[0,:])):
            color = [1.0,1.0,.8]
            color = colorWedge.copy()#[.8,.8,.7]
            color[:-1]*=.6
            I = tSteps_list_Fig02[iSim,iStep]
            dI = tSteps_list_Fig02[iSim,1]-tSteps_list_Fig02[iSim,0]
#            plt.plot([x1,x1-(x1-x0)*.025],time_list[I]/kyr*arr([1.0,1.0]),'-k',linewidth=1.0)
#            plt.plot([x0,x1],time_list[I]/kyr*arr([1.0,1.0]),':k',linewidth=.5)
            plt.plot([x0,(Setup.Grid.nxC-xFront[I])*dx],time_list[I]/kyr*arr([1.0,1.0]),'-',color=color,linewidth=.5)
            if iSim==2:
                if iStep<len(tSteps_list_Fig02[0,:])-1:
#                    plt.text((Setup.Grid.nxC-xFront[I])*dx*0.975,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'$t_{%i}$'%(iStep+1),color=color,size=12)
                    plt.text((Setup.Grid.nxC-xFront[I+int(.25*dI)])*dx+1.5,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'$t_{%i}$'%(iStep+1),color=color,size=12)
                else:
#                    plt.text((Setup.Grid.nxC-xFront[I])*dx*0.975,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'$t_{ref}$',color=color,size=12)
                    plt.text(-xmin-.2,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'$t_{ref}$',color=color,size=12)


            # write legend
            if iSim==2:
                if iStep==6:
                    plt.text((Setup.Grid.nxC-xFront[I])*dx-1.6,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'front position',color='k',size=10,rotation=293,horizontalAlignment='center')
                
            if iSim==1:
                if iStep==6:
                    plt.text((Setup.Grid.nxC-xMid[I])*dx+.8,time_list[I]/kyr+time_list[lastStep-1]/kyr*.08,'fault @ y=0.8',color=colors_xFault[1,:],size=10,rotation=42,horizontalAlignment='center')
                    plt.text((Setup.Grid.nxC-xBase[I])*dx+.3,time_list[I]/kyr+time_list[lastStep-1]/kyr*.025,'fault @ y=0.0',color=colors_xFault[2,:],size=10,rotation=39,horizontalAlignment='center')
    
    
    elif iSim == 3 or iSim == 5:
        I = tSteps_list_Fig02[iSim,len(tSteps_list_Fig02[0,:])-1]
        plt.text(-xmin-.2,time_list[I]/kyr+time_list[lastStep-1]/kyr*.015,'$t_{ref}$',color=color,size=12)
    elif iSim == 4:
        I = 1408
        plt.plot([x0,(Setup.Grid.nxC-xFront[I])*dx],time_list[I]/kyr*arr([1.0,1.0]),'-',color=color,linewidth=.5)
        plt.text(-xmin-1.5,time_list[I]/kyr+time_list[I]/kyr*.015,'$t_{ref}$',color=color,size=12)
        
#        plt.ylim([x1,x0])
#        plt.xlim([abs(xmin),0])

    plt.ylim([y0,y1])
    plt.xlim([x0,x1])

    list_front = [200,300,400,500,600,700,800,900]
    steps = np.zeros(len(list_front))
    print("iSim = %i" % iSim)
    for i in range(len(list_front)):
        Istep = np.argmin(np.abs(Setup.Grid.nxC-xFront - list_front[i]))
        steps[i] = tSteps[Istep]
    print(steps)

#    print(np.round(np.linspace(0,steps[-1],9)))

    
    ax = plt.sca(AxesDrawing['%i%i' % (iRow,iCol)])
    xIntersectBase = -(Setup.Grid.nxC-xBase [lastStep-1])*dx
    xIntersectMid  = -(Setup.Grid.nxC-xMid  [lastStep-1])*dx
        
#    if iSim==0:
#        plt.plot([xmin, xmin*0.95],[0.87,0.87],'-',color='w',linewidth=2.0)
#        plt.plot([xmin, xmin*0.95],[0.05,0.05],'-',color='w',linewidth=2.0)
#        plt.plot([xmin, xmin*0.95],[0.87,0.87],'-',color=colors_xFault[1,:],linewidth=1.0)
#        plt.plot([xmin, xmin*0.95],[0.05,0.05],'-',color=colors_xFault[2,:],linewidth=1.0)
        
    plt.plot(xIntersectMid +arr([-.5,.5]),[0.87,0.87],'-',color='w',linewidth=2.0)
    plt.plot(xIntersectBase+arr([-.5,.5]),[0.05,0.05],'-',color='w',linewidth=2.0)
    plt.plot(-(Setup.Grid.nxC-xBase [lastStep-1])*arr([dx,dx]),[0.0,0.05+0.4],'-',color='w',linewidth=2.0)
    plt.plot(-(Setup.Grid.nxC-xMid  [lastStep-1])*arr([dx,dx]),[0.0,0.87+0.4],'-',color='w',linewidth=2.0)
        
    plt.plot(xIntersectMid +arr([-.3,.3]),[0.87,0.87],'-',color=colors_xFault[1,:],linewidth=1.0)
    plt.plot(xIntersectBase+arr([-.3,.3]),[0.05,0.05],'-',color=colors_xFault[2,:],linewidth=1.0)    
    plt.plot(-(Setup.Grid.nxC-xBase [lastStep-1])*arr([dx,dx]),[0.0,0.05+0.4],'-',color=colors_xFault[2,:],linewidth=1.0)
    plt.plot(-(Setup.Grid.nxC-xMid  [lastStep-1])*arr([dx,dx]),[0.0,0.87+0.4],'-',color=colors_xFault[1,:],linewidth=1.0)
        
    

    
# end iSim
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    