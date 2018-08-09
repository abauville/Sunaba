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



Compute = True
if Compute:

    
    rho_w = 0000.0
    rho = 2500.0
    phiRef   = 30.0*pi/180.0
    LambdaRef=0.60
    
    
    
    
    RefTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                     Lambda=LambdaRef, Lambda_b=LambdaRef,
                     rho_w=rho_w, rho=rho)
    RefTaper.computeAlphaVsBeta(n=2010)
#        alpha_Weak.append(np.zeros(n))
#        alpha_WeakBase_up.append(np.zeros(n))
#        alpha_WeakBase_low.append(np.zeros(n))
        
    WeakFac = 0.1
    Weak = int(100.0*WeakFac)
    beta = 0.0
    LambdaWeak = (1.0-WeakFac) * LambdaRef   + WeakFac

    WeakBaseTaper = Taper(phi=phiRef, phi_b=phiRef,
                          Lambda=LambdaRef, Lambda_b=LambdaWeak,
                          rho_w=rho_w, rho=rho)
    
    WeakTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                      Lambda=LambdaWeak, Lambda_b=LambdaWeak,
                      rho_w=rho_w, rho=rho)
    
    Ratio = 0.6
    HalfWeakTaper = Taper(phi=phiRef, phi_b=phiRef,
                          Lambda=(1.0-Ratio)*LambdaRef+Ratio*LambdaWeak, Lambda_b=LambdaWeak,
                          rho_w=rho_w, rho=rho)
    
    
    WeakBaseTaper.computeAlphaVsBeta(n=2010)
    WeakTaper.computeAlphaVsBeta(n=2010)
    HalfWeakTaper.computeAlphaVsBeta(n=2010)
        
    
    alpha_Ref  = RefTaper.findAlpha(beta,"upper")
    alpha_Ref += RefTaper.findAlpha(beta,"lower")
    alpha_Ref /= 2.0
        
    
    alpha_Weak  = WeakTaper.findAlpha(beta,"upper")
    alpha_Weak += WeakTaper.findAlpha(beta,"lower")
    alpha_Weak /= 2.0
    
    alpha_WeakBase_up  = WeakBaseTaper.findAlpha(beta,"upper")
    alpha_WeakBase_low = WeakBaseTaper.findAlpha(beta,"lower")
    
    alpha_HalfWeak_up  = HalfWeakTaper.findAlpha(beta,"upper")
    alpha_HalfWeak_low = HalfWeakTaper.findAlpha(beta,"lower")
        # end iWeak

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
    daijoubu = 1
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
    
    
plt.figure(1)
plt.clf()
colors = [ [0.2,0.3,0.8], [0.2,0.8,0.3], [1.0,0.1,0.2] ]
i=0
taper = WeakBaseTaper
labels=["$\\lambda_w, \\lambda_w$", "$\\lambda_{ref}, \\lambda_w$", "$\\lambda_{ref}, \\lambda_{ref}$"]
plt.fill(taper.beta_all *180.0/pi,taper.alpha_all*180.0/pi, color = colors[1],alpha=0.3)
for taper in (WeakTaper, WeakBaseTaper, RefTaper):
    plt.plot(taper.beta_all *180.0/pi,taper.alpha_all*180.0/pi, color = colors[i],label=labels[i])
    i+=1
    
plt.plot([0.0,0.0],[-15.0,15.0],"--k")
plt.axis([-15.0,80.0,-15.0,15.0])
plt.ylabel("$\\alpha$ [ ]")
plt.xlabel("$\\beta$ [ ]")
plt.legend()

Lambda = int(LambdaRef * 100.0)
superDirList = ["Hc0.062_Lambda%i" % Lambda]

plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/AlphaVsBeta_Weak%i_%s.png" % (Weak, superDirList[0]))




## Load Surface Data of Sim

#Weak = 50

#Lambda = 60

#    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/SurfaceAngle.npz");
loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/SurfaceAngle_Weak%i_%s.npz" % (Weak,superDirList[0]));
strainFront = loadedData["strainFront"][()]
strainBase = loadedData["strainBase"][()]
xFront = loadedData["xFront"][()]
xBase = loadedData["xBase"][()]
timeList = loadedData["timeList"][()]
slope = loadedData["slope"][()]
    
#end Compute
    
    
    
cdictRed = {'red':  ((0.0 , 1.0, 0.0),
                    (1.0 , 1.0, 0.0)),

          'green': ((0.0 , 1.0, 1.0),
                    (1.0 , 0.0, 1.0)),
 
          'blue':  ((0.0 , 1.0, 1.0),
                    (1.0 , 0.0, 0.0))
    }
cdictGray = {'red':  ((0.0 , 1.0, 0.0),
                     (1.0 , 0.0, 0.0)),

            'green': ((0.0 , 1.0, 1.0),
                      (1.0 , 0.0, 1.0)),

            'blue':  ((0.0 , 1.0, 1.0),
                      (1.0 , 0.0, 0.0))
    }


CMAP = np.zeros((256,4));
CMAP[:,0] = np.linspace(1.0,1.0,256);
CMAP[:,1] = np.linspace(1.0,0.0,256);
CMAP[:,2] = np.linspace(1.0,0.0,256);
CMAP[:,3] = np.linspace(0.25,1.0,256);
CMAP[0,3] = 0.0;


myCMAP =  matplotlib.colors.ListedColormap(CMAP,"RedTransparent",256);
#CMAP = LinearSegmentedColormap('RedTransparent', cdictRed)
plt.register_cmap(cmap=myCMAP)

CMAP = np.zeros((256,4));
CMAP[:,0] = np.linspace(1.0,0.0,256);
CMAP[:,1] = np.linspace(1.0,0.0,256);
CMAP[:,2] = np.linspace(1.0,0.0,256);
CMAP[:,3] = np.linspace(0.25,1.0,256);
CMAP[0,3] = 0.0;


myCMAP =  matplotlib.colors.ListedColormap(CMAP,"GrayTransparent",256);
#CMAP = LinearSegmentedColormap('RedTransparent', cdictRed)
plt.register_cmap(cmap=myCMAP)

plt.figure(2)
plt.clf()
plt.set_cmap("GrayTransparent")

Hsed = 1.0e3
#dx = (Setup.Grid.xmax-Setup.Grid.xmin)/Setup.Grid.nxC / Hsed

iSim = 0
yr = 3600.0*24.0*365.25
kyr = 1e3*yr
#for iSim in range(iSim0,nSim):
#plt.plot(timeList[iSim]/kyr,slope[iSim]*180.0/np.pi,'.')
smoothSlope = np.zeros(slope[iSim].shape)
smoothWindowSize = 2
for i in range(smoothSlope.size):
    i0 = i-smoothWindowSize
    i1 = i+smoothWindowSize+1
    
    if i0<0:
        i0 = 0
    if i1>smoothSlope.size:
        i1 = smoothSlope.size
    
    smoothSlope[i] = np.mean(slope[iSim][i0:i1])

cm = 0.01
km = 1e3
vBack = 10.0 * cm/yr
Hsed = 2.0 * km
shortList = timeList[0]*vBack/Hsed

iBegin = np.argmin(timeList[0])
iBegin = 6
#plt.plot(shortList[iBegin::]/kyr,slope[iSim][iBegin::]*180.0/np.pi,'.')
plt.fill([0, shortList[-1], shortList[-1], 0],
         np.array([alpha_WeakBase_low, alpha_WeakBase_low, alpha_WeakBase_up, alpha_WeakBase_up])*180.0/pi
         ,color=colors[1],alpha=0.3)
#plt.plot(shortList[iBegin::],smoothSlope[iBegin::]*180.0/np.pi,'-k',linewidth=3)
plt.plot(shortList[iBegin::],slope[0][iBegin::]*180.0/np.pi,'-k',linewidth=3)
plt.plot([0, shortList[-1]], np.array([alpha_Ref, alpha_Ref])*180.0/pi,"-",color=colors[2])
plt.plot([0, shortList[-1]], np.array([alpha_Weak, alpha_Weak])*180.0/pi,"-",color=colors[0])
plt.plot([0, shortList[-1]], np.array([alpha_WeakBase_up, alpha_WeakBase_up])*180.0/pi,"-",color=colors[1])
plt.plot([0, shortList[-1]], np.array([alpha_WeakBase_low, alpha_WeakBase_low])*180.0/pi,"-",color=colors[1])

#plt.axis([0.0, shortList[-1], 0.0, 5])
#plt.axis([0.0, 17.0, 3.0, 12.0])
plt.axis([0.0, 5.0, 0.0, 4.0])


# set a font dict
fontdict = {'family': 'Montserrat',
'weight': 'bold',
'size': 18
}

plt.xlabel("shortening",fontdict=fontdict)
plt.ylabel("$\\alpha$ [ ]",fontdict=fontdict)


matplotlib.rc('font', **fontdict)
plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/SurfaceAngle_Weak%i_%s.png" % (Weak, superDirList[0]))