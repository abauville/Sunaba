#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

## Lambda vs chi @ beta=0
plt.figure(7)
plt.clf()
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

CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,levels = np.linspace(-1.0,1.0,1000))
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,shading='Gouraud',vmin=-1.0,vmax=1.0)

# = np.zeros((nLambda,nChi))

#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff-taper_angles, levels = [-1000.0, 0.0, 1000.0])
#ax.clabel(CS, CS.levels, inline=True, fontsize=16)

plt.text(75,90,"Extensional",fontdict=font,horizontalAlignment="center",verticalAlignment="center",color="w")
plt.text(35,20,"Compressional",fontdict=font,horizontalAlignment="center",verticalAlignment="center",color="w")

plt.xlabel("$\\lambda$")
plt.ylabel("$\\chi$")

#plt.plot(alphas_width*180.0/pi)

#plt.ylim(100.0,0.0)

cbar = plt.colorbar()

ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])