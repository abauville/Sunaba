#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:45:32 2018

@author: abauville
"""

import sys
sys.path.insert(0, '../../Utils/')
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, tan, pi, arccos, arcsin, arctan
from CritTaper_utils import Taper
import time
deg = pi/180.0


# Attempt at analytical construction of the alpha vs beta plot
n = 1000
alpha = np.linspace(-45.0,45.0,n) * deg
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

alpha_p = arctan( (1.0-rho_w/rho)/(1.0-Lambda)*tan(alpha) )

psi_b = 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p
psi_0 = 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p

taperAngle = psi_b-psi_0

beta =  taperAngle-alpha
beta2= pi-(taperAngle-alpha)-0.5

beta_center  = 0.5* (np.max(beta ) + np.min(beta ))
alpha_center = 0.5* (np.max(alpha) + np.min(alpha))

plt.figure(100)
plt.clf()



tic = time.time()
for i in range(1000):
    tpr = Taper(phi=phi, phi_b=phi_b,
                Lambda=Lambda, Lambda_b=Lambda_b,
                rho_w=rho_w, rho=rho)
    tpr.computeAlphaVsBeta()
print("time = %.2f s" % (time.time()-tic))
plt.plot(tpr.beta_all/deg,tpr.alpha_all/deg,'-r')
#
#mueff = np.tan(phi_b)*((1.0-Lambda_b)/(1.0-Lambda))
#beta_Coulomb2 = pi/4.0-(np.arctan(mueff))/2.0
#t
#I = np.argmin(abs(alpha-0))
#beta0 = beta[I]
#plt.plot(beta*1.0/deg,alpha*1.0/deg,'-b')
##plt.plot(beta2*1.0/deg*.4+25.0,alpha*1.0/deg,'-b')
##plt.plot((np.pi/4.0-phi_b_p/2.0)/deg,.0,'ob')
#plt.plot(beta*1.0/deg+2.0*((np.pi/4.0-phi_b_p/2.0)-beta[I])/deg,alpha*1.0/deg,'-g')
#plt.axis('equal')
#
#
#alpha_m = np.arctan( (1.0-Lambda_b_ov)*np.tan(phi_b))
##alpha_m = np.arctan((1.0-chi)*Lambda_ov*np.tan(phi) )
##alpha_m = -np.arctan((chi)*(1.0-Lambda_ov)*np.tan(phi) )
#beta_m  = -alpha_m
#
#
#
#beta_previous = beta
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
##beta_c = np.pi/4.0-np.arctan((1.0-chi)*np.tan(phi))/2.0
#beta_c = np.pi/4.0-phi_b_p/2.0
#
#plt.cla()
#
#plt.plot(beta_previous*1.0/deg,alpha*1.0/deg,'--b')
#plt.plot(beta1*1.0/deg,alpha*1.0/deg,'--k')
#plt.plot(beta2*1.0/deg,alpha*1.0/deg,'--r')
#plt.plot(beta3*1.0/deg,alpha*1.0/deg,'--g')
#
#
#plt.plot(beta_m/deg, alpha_m/deg,'or')
#plt.plot(beta_c/deg, 0.0/deg, 'ob')
#plt.plot((beta_m+2.0*(beta_c-beta_m))/deg, -alpha_m/deg,'om')
#



#plt.xlim([-8.0,-4.0])
#plt.ylim([4.0,8.0])

#plt.plot(beta*1.0/deg+2.0*((np.pi/4.0-phi_b_p/2.0)-beta[I])/deg,alpha*1.0/deg,'-g')