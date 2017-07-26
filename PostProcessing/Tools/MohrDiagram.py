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

Tauxx = -.5
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
    psi = np.arctan(-StressRatio+sqrt(StressRatio**2+1))
else:
	psi = np.arctan(-StressRatio-sqrt(StressRatio**2+1))



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')

# Eliminate upper and right axes
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


# linewidth
Axis_width = 1.5
LabelFontSize = 15
ax.spines['left'].set_linewidth(Axis_width)
ax.spines['bottom'].set_linewidth(Axis_width)
ax.xaxis.set_tick_params(width=Axis_width)
ax.yaxis.set_tick_params(width=Axis_width)
ax.xaxis.set_tick_params(which='minor', left='on')

ax.set_xlabel('$\sigma_n$',size=LabelFontSize, x=1, va='bottom')
ax.set_ylabel('$\\tau$',size=LabelFontSize, y=1, ha='left', va='bottom', rotation=0)


xPad = abs(Tauxx)/10.0

#plt.xticks([P, Sigma3-xPad, Sigma1+xPad], ['P', '$\sigma_3$', '$\sigma_1$'], size=LabelFontSize)


plt.hold(True)
plot(P+TauII*np.cos(phi), TauII*np.sin(phi), linewidth=1.5, color='k')
plot(Sigma1,0,"xr")
plot(Sigma3,0,"xb")
plot(P+Tauxx, Tauxy, "xg")
plot(P+Tauyy, Tauyx, "xm")

plot([P,P+TauII*np.cos(2*psi)],[0,TauII*np.sin(2*psi)])

plt.axis("equal")

plt.show()




