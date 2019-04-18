#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:45:32 2018

@author: abauville
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, tan, pi, arccos, arcsin, arctan
from CritTaper_utils import Taper
deg = pi/180.0


# Attempt at analytical construction of the alpha vs beta plot
n = 1000

rho_w = 1000.0
rho = 2500.0

phi = 30.0 * deg
phi_b = 20.0* deg


Lambda_b = 0.9
Lambda = 0.6

Lambda_hydro = rho_w/rho

Lambda_ov   = 1.0 - (1.0-Lambda  )/(1.0-Lambda_hydro)
Lambda_b_ov = 1.0 - (1.0-Lambda_b)/(1.0-Lambda_hydro)

chi = 1.0-(1.0-Lambda_b)/(1.0-Lambda)



mu = tan(phi)
mu_b = tan(phi_b)
mu_b_p = mu_b*((1.0-Lambda_b)/(1.0-Lambda))
phi_b_p = arctan(mu_b_p)

alpha_m = np.arctan( (1.0-Lambda_b_ov)*np.tan(phi_b))
beta_m  = -alpha_m
beta_c  = np.pi/4.0-phi_b_p/2.0

alpha_max = np.arctan((1.0-Lambda_ov)*mu)


alpha = np.linspace(-alpha_max,alpha_max,n)# * deg

alpha_up_left   = np.linspace( alpha_m  , alpha_max,n)# * deg
alpha_up_right  = np.linspace( alpha_max,-alpha_m  ,n)# * deg
alpha_bot_left  = np.linspace( alpha_m  ,-alpha_max,n)# * deg
alpha_bot_right = np.linspace(-alpha_max,-alpha_m  ,n)# * deg

alpha_list = [alpha_m,alpha_max,-alpha_m,-alpha_max,alpha_m]
stepRef = 0.01


fac_psi_b = [-1.0,-1.0, 1.0, 1.0]
fac_psi_0 = [-1.0, 1.0,-1.0, 1.0]
fac_psi_b_pi = [0.5,0.5,1.0,0.0]
fac_psi_0_pi = [0.5,0.0,0.5,0.0]



alpha_all = np.array([])
beta_all = np.array([])

alpha_up_all = np.array([])
beta_up_all = np.array([])

alpha_bot_all = np.array([])
beta_bot_all = np.array([])
psi_b_all = np.array([])
psi_0_all = np.array([])
for i in range(len(alpha_list)-1):
#for i in range(1):
    if alpha_list[i]<alpha_list[i+1]:
        step = np.abs(stepRef)
    else:
        step = -np.abs(stepRef)
    
    if np.abs(1.0-np.abs(alpha_list[i]/alpha_max))<1e-6: 
        alpha1 = alpha_list[i]
        alpha0 = alpha_list[i+1]
    else:
        alpha0 = alpha_list[i]
        alpha1 = alpha_list[i+1]
            
        
    stepRef = step
#    alpha = np.arange(alpha_list[i],alpha_list[i+1],step)
    print(alpha1/alpha_max)
    
    alpha = np.array([])
    new_alpha = alpha_list[i]
#    print(step)
    while np.abs(new_alpha-alpha_list[i])<np.abs(alpha_list[i+1]-alpha_list[i]):
        alpha = np.append(alpha,new_alpha)
        if alpha0 == alpha_list[i]:
            step = stepRef*(1.05-np.abs((new_alpha+step-alpha0)/(alpha1-alpha0)))
        else:
            step = stepRef*(1.05-np.abs((new_alpha-alpha0)/(alpha1-alpha0)))
        new_alpha = new_alpha + step
        
        
    
        
#    print(alpha)
    
    alpha_p = arctan( 1.0/(1.0-Lambda_ov)*tan(alpha) )
    
    
    psi_b = fac_psi_b[i] * 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p + fac_psi_b_pi[i] * pi
    psi_0 = fac_psi_0[i] * 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p + fac_psi_0_pi[i] * pi
    taperAngle = psi_b-psi_0
    beta =  taperAngle-alpha

    
    if i<2:
        alpha_up_all = np.append(alpha_up_all,alpha)
        beta_up_all = np.append(beta_up_all,beta)
    elif i:
        alpha_bot_all = np.append(alpha_bot_all,alpha)
        beta_bot_all = np.append(beta_bot_all,beta)
    alpha_all = np.append(alpha_all,alpha)
    beta_all = np.append(beta_all,beta)
    psi_b_all = np.append(psi_b_all,psi_b)
    psi_0_all = np.append(psi_0_all,psi_0)
   
alpha_all = np.append(alpha_all,alpha_all[0])
beta_all = np.append(beta_all,beta_all[0])

alpha_bot_all = np.append(alpha_bot_all,alpha_bot_all[0])
beta_bot_all = np.append(beta_bot_all,beta_bot_all[0])

alpha_bot_all = np.flipud(alpha_bot_all)
beta_bot_all = np.flipud(beta_bot_all)
  
plt.cla()
plt.plot(beta_all/deg,alpha_all/deg,'-')
plt.plot(beta_up_all/deg,alpha_up_all/deg,'.r')
plt.plot(beta_bot_all/deg,alpha_bot_all/deg,'.b')

plt.clf()
plt.plot(psi_b_all/deg,'o')

#
#psi_b = 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p
#psi_0 = 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p
#taperAngle = psi_b-psi_0
#beta =  taperAngle-alpha
#
#
#psi_b = 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p + pi
#psi_0 = pi/2.0 - 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p
#taperAngle = psi_b-psi_0
#beta1 =  taperAngle-alpha
#
#psi_b = pi/2.0 - 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p
#psi_0 =  0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p
#taperAngle = psi_b-psi_0
#beta2 =  taperAngle-alpha
#
#psi_b = pi/2.0 - 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p
#psi_0 = pi/2.0 - 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p
#taperAngle = psi_b-psi_0
#beta3 =  taperAngle-alpha
#
#
##plt.plot(beta_previous*1.0/deg,alpha*1.0/deg,'-b')
#plt.plot(beta1*1.0/deg,alpha*1.0/deg,'--k')
#plt.plot(beta2*1.0/deg,alpha*1.0/deg,'--k')
#plt.plot(beta3*1.0/deg,alpha*1.0/deg,'--k')
#
#
#plt.plot(beta_m/deg, alpha_m/deg,'or')
#plt.plot(beta_c/deg, 0.0/deg, 'ob')
#plt.plot((beta_m+2.0*(beta_c-beta_m))/deg, -alpha_m/deg,'om')
#
#


#plt.xlim([-8.0,-4.0])
#plt.ylim([4.0,8.0])

#plt.plot(beta*1.0/deg+2.0*((np.pi/4.0-phi_b_p/2.0)-beta[I])/deg,alpha*1.0/deg,'-g')