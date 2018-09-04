#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""


import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr

Style = CritTaper_Style.Style()
#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

deg = 180.0/pi

fig    = Figz_Utils.Figure(6,height=13.0,mode='draft')
#fig    = Figz_Utils.Figure(6,height=13.0,mode='draft')
Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.33,leftMarginPad=1.5)
#Axes['12'].axis('off')
plt.sca(Axes['11'])




for iTaper in range(1):
    
    betas = betas_all[iTaper,:,:]
    alphas_diff = alphas_diff_all[iTaper,:,:]
    chis = chis_all[iTaper,:,:]
    
    beta_outline = np.concatenate((betas[0,:],betas[1:-2,-1],betas[-1,-1::-1],betas[-2::-1,0]))
    chi_outline = np.concatenate((chis[0,:],chis[1:-2,-1],chis[-1,-1::-1],chis[-2::-1,0]))
    
    CS = plt.contour(betas*deg, chis, alphas_diff*deg, [0.0,1e10])
    plt.plot([0.0,0.0],[0.0,1.0],"--k",linewidth=1)


    #plt.sca(ax32)
    # 1. Extract the beta values from the contour plot on ax21
    beta_contour = CS.allsegs[0][0][:,0]/180.0*pi
    chi_contour  = CS.allsegs[0][0][:,1]    
    
    alpha_contour = np.zeros(beta_contour.shape)
    for iB in range(len(beta_contour)):
        alpha_contour[iB]   = Taper_Ref[iTaper].findAlpha(beta_contour[iB],"average")




#
#plt.figure(6)
plt.cla()
n = np.int(nLambda/2.0)
x0 = -45
x1 = 105
y0 = 0
y1 = 100
plt.axis([x0,x1,y0,y1])
for iPlot in range(2):
    if iPlot==0:
        I = np.arange(0,n+1)
    else:
        I = np.arange(n,nLambda)
        
    
    betas = betas_all[I,0,:]
    alphas_Ref = alphas_Ref_all[I,0,:]
    taper_angles = betas+alphas_Ref
    Lambdas = Lambdas_Ref_all[I,0,:]
        
    
#    plt.contour(betas*180.0/pi,Lambdas,taper_angles*180.0/pi,colors=[[.8,.8,.8,1.0]],levels=np.arange(0.0,60.0,2.0))
    CS = plt.contour(betas*180.0/pi,Lambdas*100.0,taper_angles*180.0/pi,colors="k",levels=np.arange(0.0,60.0,5.0))
    plt.plot([0.0,0.0],[0.0,100.0],'--',linewidth=0.5)
    fmt = {}
    if iPlot==0:
#        ax.clabel(CS, CS.levels, inline=True)
        Ilevels = np.arange(0,CS.levels.size)
        for angle in CS.levels:
            if plt.rcParams["text.usetex"]:
                fmt[angle] = r'%.0f\°' % angle
            else:
                fmt[angle] = '%.0f°' % angle
    else:
        tapers_contour = (beta_contour+alpha_contour)*180.0/pi
        chiLevels = np.ones(CS.levels.shape)
        Ilevels = []
        i = 0
        
        for angle in CS.levels:
            chiLevels[i] = chi_contour[np.argmin(abs(tapers_contour-angle))]
            if chiLevels[i]<chi_contour[-1]:
                Ilevels.append(i)
                if plt.rcParams["text.usetex"]:
                    fmt[angle] = r'%.0f\%%' % (chiLevels[i] * 100.0)
                else:
                    fmt[angle] = '%.0f%%' % (chiLevels[i] * 100.0)
            i+=1
    
    textColor = ["r","b"]    
    cLabel = plt.clabel(CS, CS.levels[Ilevels], inline=True, fmt=fmt, fontsize=11)
    for thisText in cLabel:
        thisText._color = textColor[iPlot]
        thisText._fontproperties._family = 'Montserrat'


## Add legend
plt.text(89,90,'$\\alpha + \\beta$',verticalAlignment='baseline')
plt.text(89,82,'$\\chi_{transition}$',verticalAlignment='baseline')
plt.text(82,90,'5°',family='Montserrat',color=textColor[0],fontsize=11,verticalAlignment='baseline')
plt.text(82,82,'5%',family='Montserrat',color=textColor[1],fontsize=11,verticalAlignment='baseline')

# set a font dict
font = {'family': 'Montserrat',
'weight': 'bold',
'size': np.int(16)
}
xlabel = plt.xlabel("$\\beta$ [°]",fontdict=font,fontsize=12)
#xlabel.set_verticalalignment('center')

xlabel.set_position(xlabel.get_position()+arr([0.0,+1.0]))
ylabel = plt.ylabel("$\\lambda$ [%]",fontdict=font,fontsize=12)
ylabel.set_verticalalignment('center')

#plt.text(x0-(x1-x0)*0.08,y1-(y1-y0)*0.015,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
#plt.text(x1-(x1-x0)*0.06,y0-(y1-y0)*0.06,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)

xTicks = np.arange(-45.0,105.0,15.0)
xTickLabels = []
for i in range(xTicks.size):
    xTickLabels.append("%.0f" % xTicks[i])
plt.xticks(xTicks, xTickLabels)

#yTicks = np.linspace(0,100,11)
#yLabels = []
#for y in yTicks:
#    yLabels.append('%i' % y)
#yLabels[-1] = ''
#plt.yticks(yTicks,yLabels)