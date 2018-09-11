#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
## Lambda vs chi @ beta=0
fig    = Figz_Utils.Figure(7,height=13.0,mode='production')
#fig    = Figz_Utils.Figure(7,height=13.0,mode='draft')
#AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
#AxesDum['12'].axis('off')
Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=10.5,topMarginPad = 1.0)
#Axes   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=1.5,topMarginPad = 1.0,xPad = 2.0z)
AxesW = Axes['info']['plotsWidth']
AxesH = Axes['info']['plotsHeight']
AxesxPad = Axes['info']['xPad']
AxeslPad = Axes['info']['leftMarginPad']
AxesrPad = Axes['info']['rightMarginPad']
AxestPad = Axes['info']['topMarginPad']
AxesbPad = Axes['info']['bottomMarginPad']
cBaryPad    = 0.0
cBartPad = .75
cBarbPad = 0.0
cBarlPad = 0.4
cBarrPad = 0.0
cBarW = 0.4
#cBarAxes   = Figz_Utils.makeAxes(fig,1,aspectRatio=0.15,leftMarginPad=AxeslPad+cBarlPad,rightMarginPad=AxesrPad+1*AxesxPad+1*AxesW+cBarrPad,topMarginPad=AxesH+cBaryPad)
cBarAxes   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+AxesW+cBarlPad,
                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+AxesW+cBarlPad)-cBarW),
                                       topMarginPad=AxestPad+cBartPad,
                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)
#Axes['12'].axis('off')
plt.sca(Axes['11'])






Lambdas = Lambdas_Ref_all[:,:,0]
chis = chis_all[:,:,0]
#alphas_diff = alphas_diff_all[iTaper,:,:]
#Lambdas = np.zeros((nLambda,nChi))
alphas_diff = np.zeros((nLambda,nChi))
alphas_width = np.zeros((nLambda,nChi))
#chis = np.zeros((nLambda,nChi))
taper_angles = np.zeros((nLambda,nChi))
alphas_WB_up = np.zeros((nLambda,nChi))
beta = 0.0
for iL in range(nLambda):
    for iW in range(nChi):
        iB = np.argmin(abs(betas_all[iL,iW,:]-beta))
        alphas_diff[iL,iW] = alphas_diff_all[iL,iW,iB]
        alphas_width[iL,iW] = alphas_WB_up_all[iL,iW,iB] - alphas_WB_low_all[iL,iW,iB]
        alphas_WB_up[iL,iW] = alphas_WB_up_all[iL,iW,iB]
        taper_angles[iL,iW] = betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB]
        


CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(-1.0001,1.0001,20),vmin=-1.00,vmax=1.00)

#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff,np.linspace(-1.0001,1.0001,1001),vmin=-1.00,vmax=1.00)
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,shading='Gouraud',vmin=-1.0,vmax=1.0)

# = np.zeros((nLambda,nChi))

#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff-taper_angles, levels = [-1000.0, 0.0, 1000.0])
#ax.clabel(CS, CS.levels, inline=True, fontsize=16)

plt.text(75,90,"Extensional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)
plt.text(35,20,"Compressional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)

plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
plt.ylabel("$\\mathbf{\\chi}$ [%]",weight='bold',verticalAlignment='center')

#plt.plot(alphas_width*180.0/pi)
#


ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])






#cbar = plt.colorbar()
#cbar.set_ticks([-1.0,0.0,1.0])
plt.sca(Axes['11'])
#plt.axis([-10.0,70.0,0.0,1.0])
cbar = plt.colorbar(cax=cBarAxes['11'], ticks=[-1, 0, 1])

plt.sca(cBarAxes['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)



cbar.ax.set_xticklabels([-1.0,"",1.0])

