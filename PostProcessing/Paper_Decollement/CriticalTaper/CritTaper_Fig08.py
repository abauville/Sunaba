#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)


alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

## Lambda vs chi @ beta=0
plt.figure(8)
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
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff*180.0/pi, levels = [-1000.0, 0.0, 1000.0])
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_diff*180.0/pi,shading='flat',vmin=-20.0,vmax=20.0)
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_width*180.0/pi,shading='flat',vmin=-20.0,vmax=20.0)
plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,shading='flat',vmin=-1.0,vmax=1.0)
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_WB_up*180.0/pi,shading='flat',vmin=-40.0,vmax=40.0)
# = np.zeros((nLambda,nChi))
plt.colorbar()
#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff-taper_angles, levels = [-1000.0, 0.0, 1000.0])
#ax.clabel(CS, CS.levels, inline=True, fontsize=16)

#plt.text(75,75,"Type 2",fontdict=font,horizontalAlignment="center",verticalAlignment="center",color="w")
#plt.text(25,25,"Type 1",fontdict=font,horizontalAlignment="center",verticalAlignment="center",color="w")

plt.xlabel("$\\lambda$")
plt.ylabel("$\\chi$")

#plt.plot(alphas_width*180.0/pi)

plt.ylim(100.0,0.0)