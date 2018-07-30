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

from CritTaper_utils import Taper

n = 9
LambdaRef_list = np.linspace(1e-10,1.0-1e-10,n)

alpha_Ref           = np.zeros(n)
alpha_Weak          = np.zeros(n)
alpha_WeakBase_up   = np.zeros(n)
alpha_WeakBase_low  = np.zeros(n)


for iTaper in range(n):
    
    rho_w = 0000.0
    rho = 2500.0
    phiRef   = 30.0*pi/180.0
    LambdaRef=LambdaRef_list[iTaper]
    
    WeakFac = 0.5
    LambdaWeak = (1.0-WeakFac) * LambdaRef   + WeakFac
    
    
    RefTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                     Lambda=LambdaRef, Lambda_b=LambdaRef,
                     rho_w=rho_w, rho=rho)
    
    WeakBaseTaper = Taper(phi=phiRef, phi_b=phiRef,
                          Lambda=LambdaRef, Lambda_b=LambdaWeak,
                          rho_w=rho_w, rho=rho)
    
    WeakTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                      Lambda=LambdaWeak, Lambda_b=LambdaWeak,
                      rho_w=rho_w, rho=rho)
    
    RefTaper.computeAlphaVsBeta(n=2010)
    WeakBaseTaper.computeAlphaVsBeta(n=2010)
    WeakTaper.computeAlphaVsBeta(n=2010)
    
    beta = 0.0
    alpha_Ref[iTaper]  = RefTaper.findAlpha(beta,"upper")
    alpha_Ref[iTaper] += RefTaper.findAlpha(beta,"lower")
    alpha_Ref[iTaper] /= 2.0
    
    alpha_Weak[iTaper]  = WeakTaper.findAlpha(beta,"upper")
    alpha_Weak[iTaper] += WeakTaper.findAlpha(beta,"lower")
    alpha_Weak[iTaper] /= 2.0
    
    alpha_WeakBase_up[iTaper]  = WeakBaseTaper.findAlpha(beta,"upper")
    alpha_WeakBase_low[iTaper] = WeakBaseTaper.findAlpha(beta,"lower")
    
# end iTaper


plt.plot(LambdaRef_list, alpha_Ref*180.0/pi, '-r')
plt.plot(LambdaRef_list, alpha_Weak*180.0/pi, '-b')
plt.plot(LambdaRef_list, alpha_WeakBase_up*180.0/pi, '-g')
plt.plot(LambdaRef_list, alpha_WeakBase_low*180.0/pi, '-g')
    
    




    

    
#    plt.clf()
#    plt.xlabel("base angle []")
#    plt.ylabel("surface angle []")
#    plt.axis([0,30,0,50])
#    plt.plot(RefTaper.beta_all/pi*180.0,RefTaper.alpha_all/pi*180.0,color='r')
#    plt.plot(WeakBaseTaper.beta_all/pi*180.0,WeakBaseTaper.alpha_all/pi*180.0,color='g')
#    plt.plot(WeakTaper.beta_all/pi*180.0,WeakTaper.alpha_all/pi*180.0,color='b')

#
#
#
#beta = 0.0
#alpha = findAlpha(beta,"upper",alpha_all,beta_all)
#print("alpha = " + str(alpha*180.0/pi) + "°")
#alpha = findAlpha(beta,"lower",alpha_all,beta_all)
#print("alpha = " + str(alpha*180.0/pi) + "°")
#
#
#
#
##plt.figure(1)
##plt.clf()
##plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.7,0.7,0.75,.5])
#
#
###Lambda=0.5
##Lambda_b=0.76
##alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
##plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.75,0.7,0.7,.5])
##
##Lambda_b=0.7414
##alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
##plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.7,0.75,0.7,.5])
#plt.figure(2)
#plt.clf()
#plt.xlabel("base angle []")
#plt.ylabel("surface angle []")
#plt.axis([0,30,0,50])
#
#
##plt.plot(beta_all/pi*180.0,(alpha_all+beta_all)/pi*180.0)
#plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,color='r')
#
#
#
#
#
#
#
##psi_b = pi/4.0 - phi/2.0# Coulomb orientation
#WeakFac = 0.50
#
#phi = phi0
#phi_b = phi_b0
#Lambda0 = Lambda
#Lambda_b = Lambda_b0
#
#Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
#
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b,rho_w=rho_w)
#
#plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'-y')
#
#
##phi = (1.0-WeakFac)*phi
#
#phi = phi0
#phi_b = phi_b0
#Lambda = Lambda0
#Lambda_b = Lambda_b0
#
#Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
#Lambda = (1.0-WeakFac)*Lambda0+WeakFac
##phi_b = (1.0-WeakFac)*phi_b
#
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b,rho_w=rho_w)
#
#plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'--y')
#
#
#
#
##
##phi = phi0
##phi_b = phi_b0
##Lambda = Lambda0
##Lambda_b = Lambda_b0
###WeakFac = 1.0-np.arcsin(1.0-0.5)
##
##b = phi
##A = 1.0/b*np.arctan((1.0-WeakFac)*tan(b))
##WeakFac = 1.0-A
###Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
###Lambda = (1.0-WeakFac)*Lambda0+WeakFac
##phi_b = (1.0-WeakFac)*phi - 1e-6
##
##alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
##
##plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'-k')
##
##
##hi = phi0
##phi_b = phi_b0
##Lambda = Lambda0
##Lambda_b = Lambda_b0
###WeakFac = 1.0-np.arcsin(1.0-0.5)
###WeakFac = 0.47
###Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
###Lambda = (1.0-WeakFac)*Lambda0+WeakFac
##
##phi = (1.0-WeakFac)*phi
##phi_b = phi - 1e-6
##
##alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
##
##plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'--k')
##
##
#




