#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 12:59:21 2018

@author: abauville
"""

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
import numpy.matlib as matlib
from CritTaper_utils import Taper
nW = 50
nBeta = 50
nLambda = 11
LambdaRef_list =np.linspace(1e-10,1.0-1e-10,nLambda)

alpha_Ref    = np.zeros(nLambda)
psi_bmin_Ref = np.zeros(nLambda)
psi_bmax_Ref = np.zeros(nLambda)


Weak_list = np.linspace(0.01,0.99,nW)

beta = 0.0

enveloppeRes = 1501



Taper_WB = []
Taper_WF = [] # Fully weakened
Taper_Ref = []

deg = 180.0/pi


alphas_Ref = np.zeros(nBeta)
alphas_WB_up = np.zeros(nBeta)
alphas_WB_low = np.zeros(nBeta)
alphas_Diff = np.zeros(nBeta)

betas_all = np.zeros((nLambda,nW,nBeta))
Weak_all = np.zeros((nLambda,nW,nBeta))
#alphas_Diff_all = np.zeros((nLambda,nW,nBeta))
alphas_Ref_all = np.zeros((nLambda,nW,nBeta))
alphas_WB_up_all = np.zeros((nLambda,nW,nBeta))
alphas_WB_low_all = np.zeros((nLambda,nW,nBeta))


Weak_small = Weak_list.copy()
Weak_small = matlib.repmat(Weak_small,nBeta,1)
Weak_small = Weak_small.T


Compute = False
if Compute:
    Counter = 0
    maxCounter = nLambda*nW
    for iTaper in range(nLambda):
        Weak_all[iTaper,:,:] = Weak_small
        
        if Counter%10==0:
            print("Counter = %i/%i" % (Counter,maxCounter))
        
        rho_w = 1000.0
        rho = 2500.0
        phiRef   = 30.0*pi/180.0
        LambdaRef=LambdaRef_list[iTaper]
        
        
        
        ## ============= RefTaper =================    
        thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                         Lambda=LambdaRef, Lambda_b=LambdaRef,
                         rho_w=rho_w, rho=rho)
        thisTaper.computeAlphaVsBeta(n=2010)
        
        betaMinRef = np.min(thisTaper.beta_all)
        betaMaxRef = np.max(thisTaper.beta_all)
        
        Taper_Ref.append(thisTaper)
        ## ========================================
        
        for iWeak in range(nW):

            
            WeakFac = Weak_list[iWeak]
            LambdaWeak = (1.0-WeakFac) * LambdaRef   + WeakFac
        
            thisTaper = Taper(phi=phiRef, phi_b=phiRef,
                                  Lambda=LambdaRef, Lambda_b=LambdaWeak,
                                  rho_w=rho_w, rho=rho)
            
            
            
            
            thisTaper.computeAlphaVsBeta(n=enveloppeRes)
            
            betaMinWB = np.min(thisTaper.beta_all)
            betaMaxWB = np.max(thisTaper.beta_all)
                               
            Taper_WB.append(thisTaper)
            
            betaMin = np.max([betaMinRef,betaMinWB])
            betaMax = np.min([betaMaxRef,betaMaxWB])
            

            betas_all[iTaper,iWeak,:] = np.linspace(betaMin,betaMax, nBeta)
            
            
            ## Fully weakened
            thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                              Lambda=LambdaWeak, Lambda_b=LambdaWeak,
                              rho_w=rho_w, rho=rho)            
            thisTaper.computeAlphaVsBeta(n=enveloppeRes)
            Taper_WF.append(thisTaper)
            
            
            for iB in range(nBeta):
                beta = betas_all[iTaper,iWeak,iB]
                alphas_Ref[iB]  = Taper_Ref[iTaper].findAlpha(beta,"upper")
                alphas_Ref[iB] += Taper_Ref[iTaper].findAlpha(beta,"lower")
                alphas_Ref[iB] /= 2.0
                
                alphas_WB_up[iB] = Taper_WB[iTaper*nW+iWeak].findAlpha(beta,"upper",tol=1e-3)
                alphas_WB_low[iB] = Taper_WB[iTaper*nW+iWeak].findAlpha(beta,"lower",tol=1e-3)
            # end iB
                
#                alphas_Diff[iB] = alphas_Ref[iB] - alphas_WB_up[iB]
            
            alphas_Ref_all[iTaper,iWeak,:] = alphas_Ref
            alphas_WB_up_all[iTaper,iWeak,:] = alphas_WB_up
            alphas_WB_low_all[iTaper,iWeak,:] = alphas_WB_low            
#            alphas_Diff_all[iTaper,iW,:] = alphas_Diff
            
            Counter+=1
            # end iWeak
    # end for iTaper
    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz",
             betas_all = betas_all,
             alphas_Ref_all = alphas_Ref_all,
             alphas_WB_up_all = alphas_WB_up_all,
             alphas_WB_low_all = alphas_WB_low_all,
             Weak_all = Weak_all,
             Taper_Ref = Taper_Ref,
             Taper_WB = Taper_WB,
             Taper_WF = Taper_WF
             )
    
else: #if Compute   
    daijoubu=1
    loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz");
    betas_all = loadedData["betas_all"][()]
    alphas_Ref_all = loadedData["alphas_Ref_all"][()]
    alphas_WB_up_all = loadedData["alphas_WB_up_all"][()]
    alphas_WB_low_all = loadedData["alphas_WB_low_all"][()]
    Weak_all = loadedData["Weak_all"][()]
    Taper_Ref = loadedData["Taper_Ref"][()]
    Taper_WB = loadedData["Taper_WB"][()]
    Taper_WF = loadedData["Taper_WF"][()]




plt.figure(1)
plt.clf()
#plt.subplot(212)
#
#plt.subplot(211)


edgeColor = ["r","r","r"]
edgeColorWeak = ["g","m"]
faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
linestyle = ["-","-"]
i = 0

#for tpr in (Taper_WB[0],RefTaper):
#    plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,alpha=1.0,facecolor=faceColor[i])
#    plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,facecolor="None",edgecolor=edgeColor[i])
#    i+=1


#plt.plot(betas*deg,alphas_Diff*deg)
#plt.plot(betas*deg,alphas_Ref*deg)
#plt.plot(betas*deg,alphas_WB_up*deg)
#plt.pcolor(betas_all*deg, Weak_all, alphas_Diff_all*deg)

#plt.contour(betas_all*deg, Weak_all, alphas_Diff_all*deg)


#plt.subplot(212)
iCount = 0
#colors = ["r","g","b","y","m"]
colors = np.random.rand(nLambda,4)
colors[:,-1] = 1.0

ax11 = plt.subplot(231)
ax12 = plt.subplot(232)
ax13 = plt.subplot(233)
#
#ax21 = plt.subplot(434)
#ax22 = plt.subplot(435)
#ax23 = plt.subplot(436)
#
#ax31 = plt.subplot(437)
#ax32 = plt.subplot(438)
#ax33 = plt.subplot(439)

ax21 = plt.subplot(223)
ax22 = plt.subplot(224)

chiList = [0.25, 0.8]
#axList = [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33]
axList = [ax11, ax12, ax13]
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all
AxCount = 0
for iTaper in [0,6,9]:
    
    betas = betas_all[iTaper,:,:]
    alphas_diff = alphas_diff_all[iTaper,:,:]
    alphas_Ref = alphas_Ref_all[iTaper,:,:]
    alphas_WB_up = alphas_WB_up_all[iTaper,:,:]
    alphas_WB_low = alphas_WB_low_all[iTaper,:,:]
    chis = Weak_all[iTaper,:,:]
    
    beta_outline = np.concatenate((betas[0,:],betas[1:-2,-1],betas[-1,-1::-1],betas[-2::-1,0]))
    chi_outline = np.concatenate((chis[0,:],chis[1:-2,-1],chis[-1,-1::-1],chis[-2::-1,0]))

    for iSub in range(len(chiList)):
        plt.sca(axList[AxCount])
        
        I = np.argmin(abs(Weak_list-chiList[iSub]))
        i=0
        for tpr in (Taper_WB[iTaper*nW+I],Taper_Ref[iTaper]):
#            plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,alpha=0.5,facecolor=faceColor[i])
#            plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,facecolor="None",edgecolor=edgeColor[i])
            
            if i==1:
                color = edgeColor[i]
            else:
                color = edgeColorWeak[iSub]
                
            plt.fill(tpr.beta_all*deg, (tpr.alpha_all+tpr.beta_all)*deg,alpha=0.08,facecolor=color)
            plt.fill(tpr.beta_all*deg, (tpr.alpha_all+tpr.beta_all)*deg,facecolor="None",edgecolor=color,linestyle=linestyle[iSub])
            i+=1
        
        x0 = 30.0-70.0
        x1 = 30.0+70.0
#        y0 = -45.0
#        y1 = 45.0
        y0 = 0.0
        y1 = 90.0
        plt.axis([x0,x1,y0,y1])
        
        plt.plot([30.0,30.0],[y0,y1],':k',linewidth=0.5)
#        plt.plot([0.0,0.0],[y0,y1],':k',linewidth=0.5)
        plt.plot([x0,x1],[30.0,30.0],':k',linewidth=0.5)        
#        if iSub==1:
        Lambda = LambdaRef_list[iTaper]
        chi = chiList[iSub]
        plt.text(x0+0.025*(x1-x0),y1-0.05*(y1-y0),"$\\lambda=$%i,  $\\chi=$%i " % (int(Lambda*100.0), int(chi*100.0)))
    AxCount+=1
    
    
    
    plt.sca(ax21)
#    plt.pcolor(betas*deg, chis, alphas_diff*deg)
#    plt.pcolor(betas*deg, chis, (alphas_WB_up)*deg)
#    cbar = plt.colorbar()
    
#    plt.pcolor(betas*deg, chis, (alphas_Ref)*deg)
#    plt.pcolor(betas*deg, chis, (alphas_WB_up-alphas_WB_low)*deg)
#    cbar = plt.colorbar()
    plt.contour(chis, betas*deg, alphas_diff*deg, [0.0,1e10],colors=[colors[iTaper,:]])
#    plt.plot (beta_outline*deg,chi_outline,color=colors[iTaper,:],linestyle='--',linewidth = 0.5)
#    plt.plot([0.0,0.0],[0.0,1.0],"--k")
#    plt.plot([8.0,8.0],[0.0,1.0],"--k")
    
    for iSub in range(len(chiList)):
        plt.plot([chiList[iSub],chiList[iSub]],[-360.0,+360.0],"--k",linewidth=1)
    
    
    
plt.sca(ax11)
plt.ylabel("taper angle")
plt.xlabel("beta")
    
plt.sca(ax21)
plt.axis([0.0,1.0,-50,70])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.ylabel("beta")
plt.xlabel("chi")

plt.sca(ax22)
# Note : for this graph
# 1. Extract the beta values from the contour plot on ax21
# 2. tpr.findAlpha for those beta values
# 3. Plot the taper
# It looks like the taper angle of transition is independent of Lambda

plt.axis([0.0,1.0,-50,70])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.ylabel("taper")
plt.xlabel("chi")

#iCount+=1




