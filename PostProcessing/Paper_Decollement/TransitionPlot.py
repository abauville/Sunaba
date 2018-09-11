#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:56:03 2018

@author: abauville
"""
import sys
sys.path.insert(0, './CriticalTaper')
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import CritTaper_Style
import Figz_Utils
from numpy import array as arr


Style = CritTaper_Style.Style()

#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

deg = 180.0/pi

#fig = Figz_Utils.Figure(4,mode="draft",height=9.25)
fig = Figz_Utils.Figure(11,height=20.25)
plt.clf()
myAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00)

#
#chiShortList = [0.1, 0.3, .6, .999]
#plt.clf()
plt.sca(myAxes['11'])
LambdaShortList = np.linspace(0.0,1.0,25)
nLambdaShort = len(LambdaShortList)
#LambdaShortList = arr([0.6])
ILambdas = []
for Lambda in LambdaShortList:
    ILambdas.append(np.argmin(np.abs(LambdaRef_list-Lambda)))
    
beta = 0
alphas_WB_up = np.zeros(nChi)
alphas_WB_low = np.zeros(nChi)
alphas_Ref = np.zeros(nChi)
alphas_WF = np.zeros(nChi)
alphas_WB_up_max = np.zeros(nLambdaShort)
chis_WB_up_max = np.zeros(nLambdaShort)
alphas_WB_maxWidth = np.zeros(nLambdaShort)
chis_WB_maxWidth = np.zeros(nLambdaShort)

alphas_RefWB_up  = np.zeros(nLambdaShort)
chis_RefWB_up = np.zeros(nLambdaShort)

il = -1
for iL in ILambdas:
    il += 1
    for iC in range(nChi):
        iB = np.argmin(abs(betas_all[iL,iC,:]-beta))
        alphas_WB_up[iC] = alphas_WB_up_all[iL,iC,iB]
        alphas_WB_low[iC] = alphas_WB_low_all[iL,iC,iB]
        alphas_Ref[iC] = alphas_Ref_all[iL,iC,iB]
        alphas_WF[iC] = alphas_WF_all[iL,iC,iB]
    # end iC
    
    plt.plot(chi_list,alphas_WB_up,'-k',linewidth=0.5)
    I = np.argmax(alphas_WB_up)
    alphas_WB_up_max[il] = alphas_WB_up[I]
    chis_WB_up_max[il] = chi_list[I]
#    plt.plot(chi_list[I],alphas_WB_up[I],'o')
    
    plt.plot(chi_list,alphas_WB_up-alphas_WB_low,'-b',linewidth=0.5)
    I = np.argmax(alphas_WB_up-alphas_WB_low)
    alphas_WB_maxWidth[il] = alphas_WB_up[I]-alphas_WB_low[I]
    chis_WB_maxWidth [il] = chi_list[I]
    
    I = np.argmin(np.abs(alphas_WB_up-alphas_Ref))
    alphas_RefWB_up[il] = alphas_WB_up[I]-alphas_Ref[I]
    chis_RefWB_up [il] = chi_list[I]
#    plt.plot(chi_list[I],alphas_WB_up[I]-alphas_WB_low[I],'o')
    
# end iL
    
plt.plot(chis_WB_up_max,alphas_WB_up_max,'k')
plt.plot(chis_WB_maxWidth,alphas_WB_maxWidth,'k')

plt.sca(myAxes['12'])
#plt.plot(LambdaRef_list[ILambdas],chis_WB_up_max,'k')
#plt.plot(LambdaRef_list[ILambdas],chis_WB_maxWidth,'b')
plt.plot(LambdaRef_list[ILambdas],chis_RefWB_up,'r')
