#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
from numpy import pi,sin,cos
import matplotlib.pyplot as plt
import CritTaper_dataMaker
#from CritTaper_utils import Taper
import CritTaper_Style
import Figz_Utils
from numpy import array as arr
#from CritTaper_WedgeVisu import plotWedge


## Create window, Style, etc...
Style = CritTaper_Style.Style()

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_WB_up_all - alphas_Ref_all
deg = 180.0/pi


#fig = Figz_Utils.Figure(5,mode="draft",height=11.0)
fig         = Figz_Utils.Figure(5,height=11.0)
Axes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.0,leftMarginPad=1.5)
#Axes['11'].axis('off')
AxesW = Axes['info']['plotsWidth']
AxesH = Axes['info']['plotsHeight']
AxesxPad = Axes['info']['xPad']
AxeslPad = Axes['info']['leftMarginPad']
AxesrPad = Axes['info']['rightMarginPad']
cBaryPad    = -1.1
cBarlPad = 1.5
cBarrPad = 0.25
cBarAxes   = Figz_Utils.makeAxes(fig,1,aspectRatio=0.1,leftMarginPad=AxeslPad+cBarlPad,rightMarginPad=AxesrPad+2*AxesxPad+2*AxesW+cBarrPad,topMarginPad=AxesH+cBaryPad)

edgeColor = ["r","r","r"]
edgeColorWeak = ["g","m"]
faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
linestyle = ["-","-"]
i = 0


#plt.subplot(212)
iCount = 0
#colors = ["r","g","b","y","m"]
colors = np.random.rand(nLambda,4)
colors[:,-1] = 1.0


chiList = [0.25, 0.99]
#axList = [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33]
axList = [Axes['11'],Axes['12'],Axes['13']]

for iTaper in range(nLambda):
    
    betas = betas_all[iTaper,:,:]
    alphas_diff = alphas_diff_all[iTaper,:,:]
    alphas_Ref = alphas_Ref_all[iTaper,:,:]
    alphas_WB_up = alphas_WB_up_all[iTaper,:,:]
    alphas_WB_low = alphas_WB_low_all[iTaper,:,:]
    chis = chis_all[iTaper,:,:]
    taper_angles = betas+alphas_Ref
    beta_outline = np.concatenate((betas[0,:],betas[1:-2,-1],betas[-1,-1::-1],betas[-2::-1,0]))
    chi_outline = np.concatenate((chis[0,:],chis[1:-2,-1],chis[-1,-1::-1],chis[-2::-1,0]))

    if iTaper in [nLambda-6]:

        plt.sca(Axes['11'])
#        plt.pcolor(betas*deg, chis, alphas_diff/taper_angles,vmin=-1.0,vmax=1.0)
        plt.contourf(betas*deg, chis, alphas_diff/taper_angles,np.linspace(-1.0,1.0,1000),vmin=-1.00001,vmax=1.00001)
        plt.set_cmap("seismic")


    if iTaper in range(0,nLambda,5):
        plt.sca(Axes['12'])
#        CS = plt.contour(betas*deg, chis, alphas_diff*deg, [0.0,1e10],colors=[colors[iTaper,:]])
        CS = plt.contour(betas*deg, chis, alphas_diff, [0.0,1e10],colors='k')
#        plt.plot([0.0,0.0],[0.0,1.0],"--k",linewidth=1)
        Lambda = Lambdas_Ref_all[iTaper,0,0]        
        fmt = {}
        for level in CS.levels:
            if iTaper == 0:
                if plt.rcParams["text.usetex"]:
                    fmt[level] = r'$\\lambda=$%.0f \%%' % (Lambda * 100.0)
                else:
                    fmt[level] = '$\\lambda=$%.0f %%' % (Lambda * 100.0)
            else:
                if plt.rcParams["text.usetex"]:
                    fmt[level] = r'%.0f \%%' % (Lambda * 100.0)
                else:
                    fmt[level] = '%.0f %%' % (Lambda * 100.0)
        
        cLabel = plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)
            
    if iTaper==0:            
        plt.sca(Axes['13'])
        # 1. Extract the beta values from the contour plot on ax21
        beta_contour = CS.allsegs[0][0][:,0]/180.0*pi
        chi_contour  = CS.allsegs[0][0][:,1]    
        # 2. tpr.findAlpha for those beta values
        alpha_contour = np.zeros(beta_contour.shape)
        for iB in range(len(beta_contour)):
            alpha_contour[iB]   = Taper_Ref[iTaper].findAlpha(beta_contour[iB],"average")
            
        plt.plot((beta_contour+alpha_contour)*180.0/pi,chi_contour,"-k")
        plt.text(10,0.525,'for all $\\lambda$',rotation=55)
        

    

plt.sca(Axes['11'])
plt.axis([-10.0,70.0,0.0,1.0])
cbar = plt.colorbar(orientation='horizontal',cax=cBarAxes['11'], ticks=[-1, 0, 1])
plt.sca(cBarAxes['11'])
#plt.text(0.5,-1.2,"$\\frac{\\Delta \\alpha}{(\\alpha_{ref}+\\beta)}$",horizontalAlignment='center',size=12)
plt.text(0.5,+1.2,"$\\bar{\\Delta \\alpha}$",horizontalAlignment='center',size=11)
#plt.text(0.25,0.5,"ext.",horizontalAlignment='center',verticalAlignment='center',size=10,color='w')
#plt.text(0.75,0.5,"compr.",horizontalAlignment='center',verticalAlignment='center',size=10,color='w')
#cbar.ax.set_xticks([-1.0,0.0,1.0])
#cbar.ax.set_xticks([-1.0,0.0,1.0])
cbar.ax.set_xticklabels([-1,0,1])


#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.sca(Axes['11'])
plt.text(-3,0.8,"extensional",rotation = 80,color = 'w')
plt.text(20.0,0.6,"compressional",rotation = 0)
plt.xlabel("$\\beta$ [°]")
plt.ylabel("$\\chi$")
#    
plt.sca(Axes['12'])
plt.axis([-45,25,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("$\\beta$  [°]")


plt.sca(Axes['13'])
plt.axis([0.0,25.0,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("$\\alpha+\\beta$  [°]")
#iCount+=1

for ax in axList:
        
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


Axes['12'].axes.get_yaxis().set_ticklabels([])
Axes['13'].axes.get_yaxis().set_ticklabels([])
