#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import sys
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
## Lambda vs chi @ beta=0
fig    = Figz_Utils.Figure(77,height=13.0,mode='production')
#fig    = Figz_Utils.Figure(7,height=13.0,mode='draft')
#AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
#AxesDum['12'].axis('off')
#Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=10.5,topMarginPad = 1.0)
Axes   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=1.5,topMarginPad = 1.0,xPad = 3.0)
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



deg = 180.00/np.pi


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
        
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,1000)

CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(-1.0001,1.0001,20),vmin=-1.00,vmax=1.00)
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff*deg,20,vmin=-20.0,vmax=20.0)
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
#cbar = plt.colorbar(cax=cBarAxes['11'])

plt.sca(cBarAxes['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)



#cbar.ax.set_xticklabels([-1.0,"",1.0])



## Add indication of max delta alpha
plt.sca(Axes['11'])
chis_alpha_diff_min = chi_list[np.argmin(alphas_diff/taper_angles,axis=1)]
chis_alpha_diff_max = chi_list[np.argmax(alphas_diff/taper_angles,axis=1)]
#Lambdas_alpha_diff_0 = LambdaRef_list[np.argmax(np.abs(alphas_diff[0:-1,:]),axis=0)]
#
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_max*100.0,'--k')
#plt.plot(Lambdas_alpha_diff_0*100.0,chi_list*100.0)
#plt.cla()
#plt.contour(alphas_diff)

#plt.sca(Axes['12'])
bDalpha = alphas_diff/taper_angles
dum = (bDalpha[1:,:]-bDalpha[0:-1,:])/((chi_list[1]-chi_list[0])*100.0)
dAlpha_dChi = (dum[:,1:] + dum[:,0:-1])/2.0
dum = (bDalpha[:,1:]-bDalpha[:,0:-1])/((LambdaRef_list[1]-LambdaRef_list[0])*100.0)
dAlpha_dLambda = (dum[1:,:] + dum[0:-1,:])/2.0

#plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, dAlpha_dChi,vmin=-.1,vmax=.1)
#plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, dAlpha_dLambda,1000,vmin=-.1,vmax=.1)
#CS = plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, -np.sqrt(dAlpha_dLambda**2+dAlpha_dChi**2),1000)
vGrad = np.sqrt(dAlpha_dLambda**2+dAlpha_dChi**2)



#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(.00,0,2),vmin=-1.00,vmax=1.00)
#plt.colorbar()
r = 5
#plt.set_cmap('plasma')
dum = (Lambdas[0:-1,:] + Lambdas[1:,:])/2.0
Lambdas_centered = (dum[:,0:-1] + dum[:,1:])/2.0

dum = (chis[0:-1,:] + chis[1:,:])/2.0
chis_centered = (dum[:,0:-1] + dum[:,1:])/2.0

#CS = plt.contour(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, vGrad,100)
#CS = plt.contour(Lambdas_centered*100.0, chis_centered*100.0, vGrad,100)


#plt.quiver(Lambdas[0:-1:r,0:-1:r]*100.0, chis[0:-1:r,0:-1:r]*100.0, dAlpha_dLambda[::r,::r],dAlpha_dChi[::r,::r],scale=.5)
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_max*100.0,'--w')
#plt.plot(Lambdas_alpha_diff_0[:-1]*100.0,chi_list[:-1]*100.0,'--w')


chi_list_centered = (chi_list[0:-1] + chi_list[1:])/2.0
LambdaRef_list_centered = (LambdaRef_list[0:-1] + LambdaRef_list[1:])/2.0
chis_vGrad_min = chi_list_centered [np.argmin(vGrad,axis=1)]



plt.sca(Axes['11'])
plt.plot(LambdaRef_list_centered*100.0,chis_vGrad_min*100.0,'--k')





#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')
#
#Lambdas_alpha_diff_max = LambdaRef_list[np.argmax(alphas_diff,axis=0)]
#Lambdas_alpha_diff_min = LambdaRef_list[np.argmin(alphas_diff,axis=0)]
#plt.plot(Lambdas_alpha_diff_max*100.0,chi_list*100.0,'--k')
#plt.plot(Lambdas_alpha_diff_min*100.0,chi_list*100.0,'--k')
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')





Type = np.zeros([nLambda-1,nChi-1])
for iL in range(nLambda-1):
    chi_boundary = chis_vGrad_min[iL]
    for iC in range(nChi-1):
        chi = chi_list_centered[iC]
        if chi<chi_boundary:
            Type[iL][iC] = 0
        else:
            if alphas_diff[iL][iC]<0.0:
                Type[iL][iC]=2
            else:
                Type[iL][iC]=1
    #end iC
#end iL
    
from matplotlib.colors import LinearSegmentedColormap
plt.sca(Axes['12'])
plt.cla()
plt.contourf(Lambdas_centered*100.0,chis_centered*100.0,Type,vmin=0,vmax=7)
CMAP = [[.99,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,1.0]]
#CMAP = (CMAP*255.0).astype("uint8")
#CMAP = getColormap(1,"myColorMap",0,CMAP=CMAP)
#plt.register_cmap(cmap=CMAP)
nThresholds = 10
#LinearSegmentedColormap.from_list(name='myColorMap2', colors = CMAP, N = nThresholds)
nThresholds = 10
colors=[(0.75, 0.15, 0.15), (1, 0.75, 0.15), (0.15, 0.75, 0.15)]
#cmap = LinearSegmentedColormap.from_list(name='custom', colors = colors, N=nThresholds)

#plt.register_cmap(name='custom')
plt.set_cmap("Set2")
#plt.colorbar()
#plt.axis([.0,100.0,100.0,.0])
#plt.contour(Lambdas[0:-1,0:-1],chis[0:-1,0:-1],Type)
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])
                


ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])







