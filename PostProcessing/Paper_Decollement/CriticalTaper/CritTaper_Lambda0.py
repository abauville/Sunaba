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

# =============================================================================
##                Function to compute the graph alpha vs beta                ##

# Some trigo functions
def csc(x):
    return 1.0/sin(x)
def sec(x):
    return 1.0/cos(x)
def cot(x):
    return 1.0/tan(x)

# F function from Dahlen, 1984, eq (21)
def F_function(psi,phi):
    return tan(2.0*psi) / (csc(phi)*sec(2.0*psi) - 1.0)
    
def computeAlphaVsBeta(phi,phi_b,n,Lambda=0.0, Lambda_b=0.0, rho_w=1000, rho=2500):
    # Computes the graph of surface vs base angle based on the graphical method described in
    # Dahlen, 1984, Noncohesive Critical Coulomb Wedges (paragraph Stable vs Unstable Wedges p. 10,128).
    # 
    # Returns a 2*N matrix with [alpha;beta]; N=2*n
    # unit is radian
    # phi = internal friction angle
    # mu_b = basal coefficient of friction = ()
    # alpha = surface angle
    # beta = base angle
    
    # psi   = stress orientation
    # psi_b = stress orientation with respect to the base direction
    # psi_0 = stress orientation with respect to the surface direction
    
    alphaP_list = np.linspace(-phi+1e-4,phi-1e-4,n) #  note: alphaMax = phi
    if np.where(alphaP_list==0)[0].shape[0] > 0: # if value 0 is found
        # then increment n by one and reconstruct alphaP_list to avoid 0 (there is no solution for the value 0)
        n=n+1
        alphaP_list = np.linspace(-phi+1e-6,phi-1e-6,n) #  note: alphaMax = phi
    
    
    
    mu_bP = tan(phi_b) * (1.0-Lambda_b)/(1.0-Lambda)
    if Lambda>0.0:
        alpha_list = arctan(tan(alphaP_list) * (1.0-Lambda)/(1.0-rho_w/rho))
#        alpha_list = arctan(tan(alphaP_list) * (1.0-Lambda))
    else:
        alpha_list = alphaP_list
                
    
    psi = np.linspace(-pi/2.0,pi/2.0,n) 
    psi_0    = np.zeros(n)
    psi_02    = np.zeros(n)

    beta_Left = np.zeros(n)
    beta_Right = np.zeros(n)
    
    F = F_function(psi,phi)


    

    for i in range(n):
        alphaP = alphaP_list[i]
        alpha  = alpha_list[i]
        # Find phi_b
        I = np.where((F[:-1]-mu_bP)/(F[1:]-mu_bP) < 0.0 )[0]            
        psi_bmin = (psi[I][0]+psi[I+1][0])/2.0
        psi_bmax = (psi[I][1]+psi[I+1][1])/2.0
        # Find phi_0
        I = np.where((F[:-1]-tan(alphaP))/(F[1:]-tan(alphaP)) < 0.0 )[0][0]        
        psi_0 = (psi[I]+psi[I+1])/2.0
        I = np.where((F[:-1]-tan(alphaP))/(F[1:]-tan(alphaP)) < 0.0 )[0][1]
        psi_02 = (psi[I]+psi[I+1])/2.0


        ApBmin = psi_bmin-psi_0 # alpha + beta = critical taper
        ApBmin2 = psi_bmin-psi_02 # alpha + beta = critical taper
        ApBmax = psi_bmax-psi_0 # alpha + beta = critical taper
        ApBmax2 = psi_bmax-psi_02 # alpha + beta = critical taper
        
        betaMin = ApBmin-alpha
        betaMin2 = ApBmin2-alpha
        betaMax = ApBmax-alpha
        betaMax2 = ApBmax2-alpha
    
        if alphaP>0.0:
            beta_Left[i] = np.max([betaMin,betaMax2])
            beta_Right[i] = betaMax
        else:
            beta_Left[i] = betaMin2
            beta_Right[i] = np.min([betaMin,betaMax2])
    
    beta_all = np.concatenate([beta_Left,np.flipud(beta_Right)])
#    alphaP_all = np.concatenate([alphaP_list,np.flipud(alphaP_list)])
    alpha_all = np.concatenate([alpha_list,np.flipud(alpha_list)])
    return (alpha_all,beta_all,psi_bmin,psi_bmax)


##                Function to compute the graph alpha vs beta                ##
# =============================================================================




n = 2010



##### Value to reproduce Dahlen, Fig 12
#phi = np.arctan(1.1)
#phi_b = np.arctan(0.85)#30.0*pi/180.0
#Lambda=0.8
#Lambda_b=0.8
#print("phi_b = " + str(phi_b*180.0/pi) + "°")
#


#### Value to reproduce Dahlen, Fig 9
#phi   = 30.0*pi/180.0
#phi_b = 10.0*pi/180.0
#Lambda=0.0
#Lambda_b=0.0
#print("phi_b = " + str(phi_b*180.0/pi) + "°")



#### Choose phi_b such that psi_b is equal to the given psi_b
#phi   = 30.0*pi/180.0
#psi_b = pi/4.0 - phi/2.0# Coulomb orientation
#phi_b = np.arctan(F_function(psi_b,phi))-0.0001
#Lambda=0.0
#Lambda_b=0.0


### Choose phi_b such that psi_b is equal to the given psi_b
rho_w = 1000.0
phi   = 30.0*pi/180.0
#psi_b = pi/4.0 - phi/2.0# Coulomb orientation
phi_b = 30.0*pi/180.0 - 1e-6
Lambda=0.864
Lambda_b=Lambda

phi0 = phi
phi_b0 = phi_b
Lambda0 = Lambda
Lambda_b0 = Lambda_b

alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b,rho_w=rho_w)
        
#plt.figure(1)
#plt.clf()
#plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.7,0.7,0.75,.5])


##Lambda=0.5
#Lambda_b=0.76
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
#plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.75,0.7,0.7,.5])
#
#Lambda_b=0.7414
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
#plt.fill(beta_all/pi*180.0,alpha_all/pi*180.0,color=[0.7,0.75,0.7,.5])
plt.figure(2)
plt.clf()
plt.xlabel("base angle []")
plt.ylabel("surface angle []")
plt.axis([0,30,0,50])


#plt.plot(beta_all/pi*180.0,(alpha_all+beta_all)/pi*180.0)
plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,color='r')







#psi_b = pi/4.0 - phi/2.0# Coulomb orientation
WeakFac = 0.50

phi = phi0
phi_b = phi_b0
Lambda0 = Lambda
Lambda_b = Lambda_b0

Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac

alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b,rho_w=rho_w)

plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'-y')


#phi = (1.0-WeakFac)*phi

phi = phi0
phi_b = phi_b0
Lambda = Lambda0
Lambda_b = Lambda_b0

Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
Lambda = (1.0-WeakFac)*Lambda0+WeakFac
#phi_b = (1.0-WeakFac)*phi_b

alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b,rho_w=rho_w)

plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'--y')




#
#phi = phi0
#phi_b = phi_b0
#Lambda = Lambda0
#Lambda_b = Lambda_b0
##WeakFac = 1.0-np.arcsin(1.0-0.5)
#
#b = phi
#A = 1.0/b*np.arctan((1.0-WeakFac)*tan(b))
#WeakFac = 1.0-A
##Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
##Lambda = (1.0-WeakFac)*Lambda0+WeakFac
#phi_b = (1.0-WeakFac)*phi - 1e-6
#
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
#
#plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'-k')
#
#
#hi = phi0
#phi_b = phi_b0
#Lambda = Lambda0
#Lambda_b = Lambda_b0
##WeakFac = 1.0-np.arcsin(1.0-0.5)
##WeakFac = 0.47
##Lambda_b = (1.0-WeakFac)*Lambda_b0+WeakFac
##Lambda = (1.0-WeakFac)*Lambda0+WeakFac
#
#phi = (1.0-WeakFac)*phi
#phi_b = phi - 1e-6
#
#alpha_all, beta_all, psi_bmin, psi_bmax = computeAlphaVsBeta(phi,phi_b,n,Lambda,Lambda_b)
#
#plt.plot(beta_all/pi*180.0,alpha_all/pi*180.0,'--k')
#
#





