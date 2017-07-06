#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 17:12:02 2017

@author: abauville
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi
from matplotlib.pyplot import plot

Tauxx = -1.0
Tauxy = .5
Tauyx = Tauxy
Tauyy = -Tauxx

P = 1.0



Sigmaxx = P+Tauxx
Sigmayy = P+Tauyy

Sigma1 = P + sqrt(P**2-(Sigmaxx*Sigmayy-Tauxy*Tauyx))
Sigma3 = P - sqrt(P**2-(Sigmaxx*Sigmayy-Tauxy*Tauyx))

TauII = sqrt(0.5*(Tauxx**2 + Tauyy**2 + Tauxy**2 + Tauyx**2))


nPoints = 100
phi = np.linspace(0,2*pi,nPoints)

StressRatio = Tauxx/Tauxy 
if (Tauxy>0.0):
    psi = 2*np.arctan(-StressRatio+sqrt(StressRatio**2+1))
else:
	psi = 2*np.arctan(-StressRatio-sqrt(StressRatio**2+1))


plt.hold(True)
plot(P+TauII*np.cos(phi), TauII*np.sin(phi))
plot(Sigma1,0,"xr")
plot(Sigma3,0,"xb")
plot(P+Tauxx, Tauxy, "xg")
plot(P+Tauyy, Tauyx, "xm")

plot([P,P+TauII*np.cos(psi)],[0,TauII*np.sin(psi)])

plt.axis("equal")





