#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 20:49:33 2018

@author: abauville
"""

# Critical Taper utilities
import numpy as np
from numpy import sin, cos, tan, pi, arctan
import matplotlib.pyplot as plt
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
    


class Taper():
    def __init__(self, phi=30.0*np.pi/180.0, phi_b=30.0*np.pi/180.0-1e-6,
                 Lambda=0.0, Lambda_b=0.0,
                 rho=2500.0, rho_w=1000.0):
        self.phi = phi
        self.phi_b = phi_b
        self.Lambda = Lambda
        self.Lambda_b = Lambda_b
        
        self.rho = rho
        self.rho_w = rho_w
        
        self.alpha_all = []
        self.beta_all = []
        self.psi_bmin = 0
        self.psi_bmax = 0
        
    # =========================================================================        
    #               Function to compute the graph alpha vs beta              
    def computeAlphaVsBeta(self, n=1001):
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
        
        phi = self.phi
        phi_b = self.phi_b
        Lambda = self.Lambda
        Lambda_b = self.Lambda_b
        
        rho_w = self.rho_w
        rho = self.rho
        
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
        
        self.beta_all = np.concatenate([beta_Left,np.flipud(beta_Right)])
        self.alpha_all = np.concatenate([alpha_list,np.flipud(alpha_list)])
        self.psi_bmin = psi_bmin
        self.psi_bmax = psi_bmax


    #               Function to compute the graph alpha vs beta              
    # =========================================================================
    
    
    # =========================================================================    
    #                      Find alpha for a given beta                       
    def findAlpha(self, beta , enveloppe, tol=1e-2):
        # enveloppe can be "upper" or "lower"
        # Find alpha for a given beta        
#        if beta>pi:
#            print("Warning: beta should be in radians")
        
        IBetaMin = np.argmin(self.beta_all)
        IBetaMax = np.argmax(self.beta_all)
    
        I0 = min(IBetaMin, IBetaMax)
        I1 = max(IBetaMin, IBetaMax)
        beta_env1 = self.beta_all[I0:I1+1]
        beta_env2 = np.concatenate((self.beta_all[I1::],self.beta_all[0:I0+1]))
    
        alpha_env1 = self.alpha_all[I0:I1+1]
        alpha_env2 = np.concatenate((self.alpha_all[I1::],self.alpha_all[0:I0+1]))
        Ienv1 = np.argmin(np.abs(beta_env1-beta))
        Ienv2 = np.argmin(np.abs(beta_env2-beta))
        
        if enveloppe=="upper":
            alpha = np.max([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
        elif enveloppe=="lower":
            alpha = np.min([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
        elif enveloppe=="average":
            alpha  = ( alpha_env1[Ienv1]+alpha_env2[Ienv2] )/2.0# min or max to get the lower or upper enveloppe of stability            
        else:
            print("unknown enveloppe")
#        alphaUp = np.max([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
#        alphaLow = np.min([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
#        plt.clf()
#        plt.plot(self.beta_all,self.alpha_all,"k")
#        plt.plot(beta_env1,alpha_env1,"--r")
#        plt.plot(beta_env2,alpha_env2,"--b") 
#        plt.plot(beta,alphaUp,"or")
#        plt.plot(beta,alphaLow,"ob")
        return alpha
    #                      Find alpha for a given beta                          
    # =========================================================================