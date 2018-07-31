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

n = 11
LambdaRef_list = np.linspace(1e-10,1.0-1e-10,n)

alpha_Ref    = np.zeros(n)
psi_bmin_Ref = np.zeros(n)
psi_bmax_Ref = np.zeros(n)

nW = 3
Weak_list = [0.5, 0.2, 0.1]

beta = 0.0
alpha_Weak = np.zeros([nW,n])
alpha_WeakBase_up = np.zeros([nW,n])
alpha_WeakBase_low = np.zeros([nW,n])

alpha_HalfWeak_up = np.zeros([nW,n])
alpha_HalfWeak_low = np.zeros([nW,n])


psi_bmin_Weak = np.zeros([nW,n])
psi_bmax_Weak = np.zeros([nW,n])

psi_bmin_WeakBase = np.zeros([nW,n])
psi_bmax_WeakBase = np.zeros([nW,n])

psi_bmin_HalfWeak = np.zeros([nW,n])
psi_bmax_HalfWeak = np.zeros([nW,n])

Compute = False
if Compute:
    for iTaper in range(n):
        
        rho_w = 0000.0
        rho = 2500.0
        phiRef   = 30.0*pi/180.0
        LambdaRef=LambdaRef_list[iTaper]
        
        
        
        
        RefTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                         Lambda=LambdaRef, Lambda_b=LambdaRef,
                         rho_w=rho_w, rho=rho)
        RefTaper.computeAlphaVsBeta(n=2010)
        alpha_Ref[iTaper]  = RefTaper.findAlpha(beta,"upper")
        alpha_Ref[iTaper] += RefTaper.findAlpha(beta,"lower")
        alpha_Ref[iTaper] /= 2.0
        psi_bmin_Ref[iTaper] = RefTaper.psi_bmin
        psi_bmax_Ref[iTaper] = RefTaper.psi_bmax
        for iWeak in range(3):
    #        alpha_Weak.append(np.zeros(n))
    #        alpha_WeakBase_up.append(np.zeros(n))
    #        alpha_WeakBase_low.append(np.zeros(n))
            
            WeakFac = Weak_list[iWeak]
            LambdaWeak = (1.0-WeakFac) * LambdaRef   + WeakFac
        
            WeakBaseTaper = Taper(phi=phiRef, phi_b=phiRef,
                                  Lambda=LambdaRef, Lambda_b=LambdaWeak,
                                  rho_w=rho_w, rho=rho)
            
            WeakTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                              Lambda=LambdaWeak, Lambda_b=LambdaWeak,
                              rho_w=rho_w, rho=rho)
            
            HalfWeakTaper = Taper(phi=phiRef, phi_b=phiRef,
                                  Lambda=(LambdaRef+LambdaWeak)/2.0, Lambda_b=LambdaWeak,
                                  rho_w=rho_w, rho=rho)
            
            
            WeakBaseTaper.computeAlphaVsBeta(n=2010)
            WeakTaper.computeAlphaVsBeta(n=2010)
            HalfWeakTaper.computeAlphaVsBeta(n=2010)
            
            psi_bmin_Weak[iWeak][iTaper] = WeakTaper.psi_bmin
            psi_bmax_Weak[iWeak][iTaper] = WeakTaper.psi_bmax

            psi_bmin_WeakBase[iWeak][iTaper] = WeakBaseTaper.psi_bmin
            psi_bmax_WeakBase[iWeak][iTaper] = WeakBaseTaper.psi_bmax
            
            psi_bmin_HalfWeak[iWeak][iTaper] = HalfWeakTaper.psi_bmin
            psi_bmax_HalfWeak[iWeak][iTaper] = HalfWeakTaper.psi_bmax
            
            alpha_Weak[iWeak][iTaper]  = WeakTaper.findAlpha(beta,"upper")
            alpha_Weak[iWeak][iTaper] += WeakTaper.findAlpha(beta,"lower")
            alpha_Weak[iWeak][iTaper] /= 2.0
            
            alpha_WeakBase_up[iWeak][iTaper]  = WeakBaseTaper.findAlpha(beta,"upper")
            alpha_WeakBase_low[iWeak][iTaper] = WeakBaseTaper.findAlpha(beta,"lower")
            
            alpha_WeakBase_up[iWeak][iTaper]  = WeakBaseTaper.findAlpha(beta,"upper")
            alpha_WeakBase_low[iWeak][iTaper] = WeakBaseTaper.findAlpha(beta,"lower")
            
            alpha_HalfWeak_up[iWeak][iTaper]  = HalfWeakTaper.findAlpha(beta,"upper")
            alpha_HalfWeak_low[iWeak][iTaper] = HalfWeakTaper.findAlpha(beta,"lower")
            
            # end iWeak
    # end for iTaper
    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper.npz",
             alpha_Ref = alpha_Ref,
             alpha_Weak = alpha_Weak,
             alpha_WeakBase_up = alpha_WeakBase_up,
             alpha_WeakBase_low = alpha_WeakBase_low,
             alpha_HalfWeak_up = alpha_HalfWeak_up,
             alpha_HalfWeak_low = alpha_HalfWeak_low,
             psi_bmin_Ref= psi_bmin_Ref,
             psi_bmax_Ref= psi_bmax_Ref,
             psi_bmin_Weak = psi_bmin_Weak,
             psi_bmax_Weak = psi_bmax_Weak,
             psi_bmin_WeakBase = psi_bmin_WeakBase,
             psi_bmax_WeakBase = psi_bmax_WeakBase,
             psi_bmin_HalfWeak = psi_bmin_HalfWeak,
             psi_bmax_HalfWeak = psi_bmax_HalfWeak
             )
    
else: #if Compute   
    
    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper.npz");
    alpha_Ref = loadedData["alpha_Ref"][()]
    alpha_Weak = loadedData["alpha_Weak"][()]
    alpha_WeakBase_up = loadedData["alpha_WeakBase_up"][()]
    alpha_WeakBase_low = loadedData["alpha_WeakBase_low"][()]
    alpha_HalfWeak_up = loadedData["alpha_HalfWeak_up"][()]
    alpha_HalfWeak_low = loadedData["alpha_HalfWeak_low"][()]
    psi_bmin_Ref= loadedData["psi_bmin_Ref"][()]
    psi_bmax_Ref= loadedData["psi_bmax_Ref"][()]
    psi_bmin_Weak = loadedData["psi_bmin_Weak"][()]
    psi_bmax_Weak = loadedData["psi_bmax_Weak"][()]
    psi_bmin_WeakBase = loadedData["psi_bmin_WeakBase"][()]
    psi_bmax_WeakBase = loadedData["psi_bmax_WeakBase"][()]
    psi_bmin_HalfWeak = loadedData["psi_bmin_HalfWeak"][()]
    psi_bmax_HalfWeak = loadedData["psi_bmax_HalfWeak"][()]
    
    
    
plt.figure(3)
Plot = True
plt.clf()
iWeak = 2
if Plot:
    colors = np.array([ [.4,.4,1.0], [1.0,.4,.4] , [.4,1.0,.4] ])
    
    
    LambdaRef_list_big = np.concatenate([LambdaRef_list, np.flipud(LambdaRef_list)])
    alpha_WeakBase_big = np.concatenate([alpha_WeakBase_up[iWeak], np.flipud(alpha_WeakBase_low[iWeak])])
    
    plt.fill(LambdaRef_list_big, alpha_WeakBase_big*180.0/pi,color=colors[iWeak],alpha=0.2,linestyle="None")
    
    plt.plot(LambdaRef_list, alpha_WeakBase_up[iWeak]*180.0/pi, color=colors[iWeak])
    plt.plot(LambdaRef_list, alpha_WeakBase_low[iWeak]*180.0/pi, color=colors[iWeak])
    
    plt.plot(LambdaRef_list, alpha_HalfWeak_up[iWeak]*180.0/pi, color=colors[iWeak])
    plt.plot(LambdaRef_list, alpha_HalfWeak_low[iWeak]*180.0/pi, color=colors[iWeak])
    
    plt.plot(LambdaRef_list, alpha_Weak[iWeak]*180.0/pi, "--", color=colors[iWeak])
    
    
        
        
    plt.plot(LambdaRef_list, alpha_Ref*180.0/pi, color="k")
    
    plt.axis([0.0,1.0,0.0,30.0])
    
    
    fontdict = {'family': 'Arial',
        'weight': 'bold',
        'size': 16
        }
    #plt.xlabel("$\\lambda$",fontdict=fontdict)
    #plt.ylabel("$\\alpha [°]$",fontdict=fontdict)
    plt.xlabel("$\\lambda$")
    plt.ylabel("$\\alpha [°]$")
    
    
    #    
    matplotlib.rc('font', **fontdict)

#end if Plot





def plotTaper(beta, alpha, psi_b, f0, fL0, fL1, flip=[0,1]):
    
    
    plt.axis("equal")
    plt.plot([0.0, cos(beta)], [0.0, sin(beta)])
    plt.plot([0.0, cos(alpha)*(2.0-cos(alpha))], [0.0, sin(alpha)*(2.0-cos(alpha))])
    
     
    cAngle = 30.0*pi/180.0 # Coulomb Angle
    fAngle0 = psi_b+cAngle+flip[0]*pi
    fAngle1 = psi_b-cAngle+flip[1]*pi
   
    plt.plot([f0[0], f0[0]+cos(fAngle0)*fL0], [f0[1], f0[1]+sin(fAngle0)*fL0])
    plt.plot([f0[0], f0[0]+cos(fAngle1)*fL1], [f0[1], f0[1]+sin(fAngle1)*fL1])



plt.figure(1)
plt.clf()
#iWeak = 2
iTaper = 6
f0 = [0.5,0.0] # origin of the fault
fL0 = 0.3 # Length of Fault 0
fL1 = 0.5 # Length of Fault 1
plt.subplot(3,2,1)
plotTaper(beta,alpha_Ref[iTaper], psi_bmin_Ref[iTaper], f0, fL0, fL1)
plt.axis([0.0,1.0,-0.05,0.4]);
plt.subplot(3,2,2)
fL0 = 0.2 # Length of Fault 0
fL1 = 0.5 # Length of Fault 1
plotTaper(beta,alpha_Weak[iWeak][iTaper], psi_bmin_Weak[iWeak][iTaper], f0, fL0, fL1)

plt.subplot(3,2,3)
plotTaper(beta,alpha_WeakBase_up[iWeak][iTaper], psi_bmax_WeakBase[iWeak][iTaper], f0, fL0, fL1,flip=[0,0])
plt.axis([0.0,1.0,-0.05,0.4]);
plt.subplot(3,2,4)
plotTaper(beta,alpha_WeakBase_low[iWeak][iTaper], psi_bmin_WeakBase[iWeak][iTaper], f0, fL0, fL1)
plt.axis([0.0,1.0,-0.05,0.4]);


plt.subplot(3,2,5)
plotTaper(beta,alpha_HalfWeak_up[iWeak][iTaper], psi_bmax_HalfWeak[iWeak][iTaper], f0, fL0, fL1,flip=[0,0])
plt.axis([0.0,1.0,-0.05,0.4]);
plt.subplot(3,2,6)
plotTaper(beta,alpha_HalfWeak_low[iWeak][iTaper], psi_bmin_HalfWeak[iWeak][iTaper], f0, fL0, fL1)
plt.axis([0.0,1.0,-0.05,0.4]);


#fL0 = 0.2 # Length of Fault 0
#fL1 = 0.5 # Length of Fault 1
#plotTaper(beta,alpha_Weak[iTaper][0], psi_bmin_Weak[iTaper][iWeak], psi_bmin_Weak[iWeak], f0, fL0, fL1)
#plt.axis([0.0,1.0,-0.05,0.4]);


#plt.subplot(3,1,2)
#iTaper = 0
#plt.axis([0.0,1.0,0.0,1.0]);
#plt.axis("equal")
#plt.plot([0.0, cos(beta)], [0.0, sin(beta)])
#plt.plot([0.0, cos(alpha_Ref[iTaper])*(2.0-cos(alpha_Ref[iTaper]))], [0.0, sin(alpha_Ref[iTaper])*(2.0-cos(alpha_Ref[iTaper]))])
#
#f0 = [0.5,0.0] # origin of the fault
#cAngle = 30.0*pi/180.0 # Coulomb Angle
#fAngle0 = psi_bmin_Ref[iTaper]+cAngle
#fAngle1 = psi_bmin_Ref[iTaper]-cAngle+pi
#fL0 = 0.5 # Length of Fault 0
#fL1 = 0.5 # Length of Fault 1
#plt.plot([f0[0], f0[0]+cos(fAngle0)*fL0], [f0[1], f0[1]+sin(fAngle0)*fL0])
#plt.plot([f0[0], f0[0]+cos(fAngle1)*fL1], [f0[1], f0[1]+sin(fAngle1)*fL1])
#    







