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

phi = 30.0 * deg
phi_b = 0.01* deg
Lambda_b = 0.4
Lambda = 0.4

n = 1000.0
alpha = np.linspace(-45.0,45.0,n) * deg
rho_w = 1000.0
rho = 2500.0


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

plt.cla()




tpr = Taper(phi=phi, phi_b=phi_b,
                Lambda=Lambda, Lambda_b=Lambda_b,
                rho_w=rho_w, rho=rho)
tpr.computeAlphaVsBeta(n=2010)

plt.plot(tpr.beta_all/deg,tpr.alpha_all/deg,'-r')

mueff = np.tan(phi_b)*((1.0-Lambda_b)/(1.0-Lambda))
beta_Coulomb2 = pi/4.0-(np.arctan(mueff))/2.0


plt.plot(beta*1.0/deg,alpha*1.0/deg,'-b')
plt.plot(beta2*1.0/deg,alpha*1.0/deg,'-b')
plt.plot(beta*1.0/deg+(pi/2.0-2.0*phi_b_p)/deg,alpha*1.0/deg,'-g')