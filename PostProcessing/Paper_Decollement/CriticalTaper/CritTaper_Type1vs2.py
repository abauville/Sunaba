#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 12:59:21 2018

@author: abauville
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 11:21:03 2018

@author: abauville
"""

## Illsutration of stress and fault orientations in accretionary prism
## Based on the critical taper theory
## Nomenclature follows Buiter, 2011, "A review of brittle compressional wedge models"


import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan
import matplotlib
from CritTaper_utils import Taper

nTapers = 1
LambdaRef_list = [0.9]#np.linspace(1e-10,1.0-1e-10,nTapers)

alpha_Ref    = np.zeros(nTapers)
psi_bmin_Ref = np.zeros(nTapers)
psi_bmax_Ref = np.zeros(nTapers)

nW = 100
Weak_list = np.linspace(0.01,0.99,nW)

beta = 0.0

enveloppeRes = 2010

Compute = True

WB_Taper = []

if Compute:
    for iTaper in range(nTapers):
        
        rho_w = 1000.0
        rho = 2500.0
        phiRef   = 30.0*pi/180.0
        LambdaRef=LambdaRef_list[iTaper]
        
        
        
        
        RefTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                         Lambda=LambdaRef, Lambda_b=LambdaRef,
                         rho_w=rho_w, rho=rho)
        RefTaper.computeAlphaVsBeta(n=2010)
        
        for iWeak in range(nW):
    #        alpha_Weak.append(np.zeros(n))
    #        alpha_WeakBase_up.append(np.zeros(n))
    #        alpha_WeakBase_low.append(np.zeros(n))
            
            WeakFac = Weak_list[iWeak]
            LambdaWeak = (1.0-WeakFac) * LambdaRef   + WeakFac
        
            thisTaper = Taper(phi=phiRef, phi_b=phiRef,
                                  Lambda=LambdaRef, Lambda_b=LambdaWeak,
                                  rho_w=rho_w, rho=rho)
            
            
            
            
            thisTaper.computeAlphaVsBeta(n=enveloppeRes)
            WB_Taper.append(thisTaper)
            
            # end iWeak
    # end for iTaper
#    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper.npz",
#             alpha_Ref = alpha_Ref,
#             alpha_Weak = alpha_Weak,
#             alpha_WeakBase_up = alpha_WeakBase_up,
#             alpha_WeakBase_low = alpha_WeakBase_low,
#             alpha_HalfWeak_up = alpha_HalfWeak_up,
#             alpha_HalfWeak_low = alpha_HalfWeak_low,
#             psi_bmin_Ref= psi_bmin_Ref,
#             psi_bmax_Ref= psi_bmax_Ref,
#             psi_bmin_Weak = psi_bmin_Weak,
#             psi_bmax_Weak = psi_bmax_Weak,
#             psi_bmin_WeakBase = psi_bmin_WeakBase,
#             psi_bmax_WeakBase = psi_bmax_WeakBase,
#             psi_bmin_HalfWeak = psi_bmin_HalfWeak,
#             psi_bmax_HalfWeak = psi_bmax_HalfWeak
#             )
    
else: #if Compute   
    daijoubu=1
#    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper.npz");
#    alpha_Ref = loadedData["alpha_Ref"][()]
#    alpha_Weak = loadedData["alpha_Weak"][()]
#    alpha_WeakBase_up = loadedData["alpha_WeakBase_up"][()]
#    alpha_WeakBase_low = loadedData["alpha_WeakBase_low"][()]
#    alpha_HalfWeak_up = loadedData["alpha_HalfWeak_up"][()]
#    alpha_HalfWeak_low = loadedData["alpha_HalfWeak_low"][()]
#    psi_bmin_Ref= loadedData["psi_bmin_Ref"][()]
#    psi_bmax_Ref= loadedData["psi_bmax_Ref"][()]
#    psi_bmin_Weak = loadedData["psi_bmin_Weak"][()]
#    psi_bmax_Weak = loadedData["psi_bmax_Weak"][()]
#    psi_bmin_WeakBase = loadedData["psi_bmin_WeakBase"][()]
#    psi_bmax_WeakBase = loadedData["psi_bmax_WeakBase"][()]
#    psi_bmin_HalfWeak = loadedData["psi_bmin_HalfWeak"][()]
#    psi_bmax_HalfWeak = loadedData["psi_bmax_HalfWeak"][()]
#    
    
nBeta = 100
alphas_Ref = np.zeros(nBeta)
alphas_WB_up = np.zeros(nBeta)
alphas_Diff = np.zeros(nBeta)

betaMinRef = np.min(RefTaper.beta_all)
betaMaxRef = np.max(RefTaper.beta_all)


betas_all = np.zeros((nW,nBeta))
alphas_Diff_all = np.zeros((nW,nBeta))
for iW in range(nW):
    WeakBaseTaper = WB_Taper[iW]
    betaMinWB = np.min(WeakBaseTaper.beta_all)
    betaMaxWB = np.max(WeakBaseTaper.beta_all)

    betaMin = np.max([betaMinRef,betaMinWB])
    betaMax = np.min([betaMaxRef,betaMaxWB])

    
    betas = np.linspace(betaMin,betaMax, nBeta)
    betas_all[iW,:] = betas
    for iB in range(nBeta):
        beta = betas[iB]
        alphas_Ref[iB]  = RefTaper.findAlpha(beta,"upper")
        alphas_Ref[iB] += RefTaper.findAlpha(beta,"lower")
        alphas_Ref[iB] /= 2.0
        
        alphas_WB_up[iB] = WeakBaseTaper.findAlpha(beta,"upper")
    
        alphas_Diff[iB] = alphas_Ref[iB] - alphas_WB_up[iB]
    
    
    alphas_Diff_all[iW,:] = alphas_Diff
    
    
    

Weak_all = Weak_list.copy()
Weak_all = np.matlib.repmat(Weak_all,nBeta,1)
Weak_all = Weak_all.T





plt.figure(1)
plt.clf()
plt.subplot(211)

deg = 180.0/pi
edgeColor = ["r","b"]
faceColor = [np.array([202,231,202])/255,[0,0,0]]
i = 0

for tpr in (WB_Taper[90],RefTaper):
    plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,alpha=1.0,facecolor=faceColor[i])
    plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,facecolor="None",edgecolor=edgeColor[i])
    i+=1

plt.subplot(212)
#plt.plot(betas*deg,alphas_Diff*deg)
#plt.plot(betas*deg,alphas_Ref*deg)
#plt.plot(betas*deg,alphas_WB_up*deg)
plt.pcolor(betas_all*deg, Weak_all, alphas_Diff_all*deg)
#plt.contour(betas_all*deg, Weak_all, alphas_Diff_all*deg, [0.0, 0.00001])
plt.contour(betas_all*deg, Weak_all, alphas_Diff_all*deg)
plt.axis([-50,90,0.0,1.0])
plt.colorbar()




