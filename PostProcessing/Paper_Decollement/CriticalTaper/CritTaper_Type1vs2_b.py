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
import matplotlib
nW = 40
nBeta = 40
nLambda = 11
LambdaRef_list =np.linspace(0,1.0,nLambda)
LambdaRef_list[ 0] += 1e-10
LambdaRef_list[-1] -= 1e-10

alpha_Ref    = np.zeros(nLambda)
psi_bmin_Ref = np.zeros(nLambda)
psi_bmax_Ref = np.zeros(nLambda)


Weak_list = np.linspace(0.01,0.99,nW)

beta = 0.0

enveloppeRes = 2001



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

Lambdas_Ref_all = np.zeros((nLambda,nW,nBeta))


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
                alphas_Ref[iB]  = Taper_Ref[iTaper].findAlpha(beta,"average")
                
                alphas_WB_up[iB] = Taper_WB[iTaper*nW+iWeak].findAlpha(beta,"upper",tol=1e-3)
                alphas_WB_low[iB] = Taper_WB[iTaper*nW+iWeak].findAlpha(beta,"lower",tol=1e-3)
                
                Lambdas_Ref_all[iTaper,iWeak,iB] = LambdaRef
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
             Lambdas_Ref_all = Lambdas_Ref_all,
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
    Lambdas_Ref_all = loadedData["Lambdas_Ref_all"][()]
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

ax11 = plt.subplot(331)
ax12 = plt.subplot(332)
ax13 = plt.subplot(333)
#
#ax21 = plt.subplot(434)
#ax22 = plt.subplot(435)
#ax23 = plt.subplot(436)
#
#ax31 = plt.subplot(437)
#ax32 = plt.subplot(438)
#ax33 = plt.subplot(439)

ax21 = plt.subplot(334)
ax22 = plt.subplot(335)
ax23 = plt.subplot(336)

ax31 = plt.subplot(337)
ax32 = plt.subplot(338)
#ax33 = plt.subplot(339)


chiList = [0.25, 0.99]
#axList = [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33]
axList = [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32]
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all
AxCount = 0
AxCount2=0

for iTaper in range(nLambda):
    
    betas = betas_all[iTaper,:,:]
    alphas_diff = alphas_diff_all[iTaper,:,:]
    alphas_Ref = alphas_Ref_all[iTaper,:,:]
    alphas_WB_up = alphas_WB_up_all[iTaper,:,:]
    alphas_WB_low = alphas_WB_low_all[iTaper,:,:]
    chis = Weak_all[iTaper,:,:]
    
    beta_outline = np.concatenate((betas[0,:],betas[1:-2,-1],betas[-1,-1::-1],betas[-2::-1,0]))
    chi_outline = np.concatenate((chis[0,:],chis[1:-2,-1],chis[-1,-1::-1],chis[-2::-1,0]))

    if iTaper in [0,6,9]:
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
        
    
        plt.sca(axList[AxCount+3])
    
#        plt.sca(ax21)
        plt.pcolor(betas*deg, chis, alphas_diff*deg,vmin=-20.0,vmax=20.0)
        plt.set_cmap("seismic")
        plt.axis([-40,100,0.0,1.0])
        AxCount+=1
#        plt.colorbar()
#        plt.contour(chis, betas*deg, alphas_diff*deg, [0.0,1e10],colors=[[0.0,0.0,0.0,1.0]])
        
    if iTaper in range(0,nLambda,1):
        plt.sca(ax31)
        CS = plt.contour(betas*deg, chis, alphas_diff*deg, [0.0,1e10],colors=[colors[iTaper,:]])
        
    
#        for iSub in range(len(chiList)):
#            plt.plot([chiList[iSub],chiList[iSub]],[-360.0,+360.0],"--k",linewidth=1)
        plt.plot([0.0,0.0],[0.0,1.0],"--k",linewidth=1)
            
            
        plt.sca(ax32)
        # 1. Extract the beta values from the contour plot on ax21
        beta_contour = CS.allsegs[0][0][:,0]/180.0*pi
        chi_contour  = CS.allsegs[0][0][:,1]    
        # 2. tpr.findAlpha for those beta values
        alpha_contour = np.zeros(beta_contour.shape)
        for iB in range(len(beta_contour)):
            alpha_contour[iB]   = Taper_Ref[iTaper].findAlpha(beta_contour[iB],"average")
            
        plt.plot((beta_contour+alpha_contour)*180.0/pi,chi_contour,".")
        # 3. Plot the taper
        # It looks like the taper angle of transition is independent of Lambda
        
        if iTaper in [0,6,9]:
            
            for iC in range(len(chiList)):
                I = np.argmin(abs(chi_contour-chiList[iC]))
                plt.sca(axList[AxCount2])
                plt.plot(beta_contour[I]*180.0/pi, (beta_contour[I]+alpha_contour[I])*180.0/pi,"ok",markerFaceColor="None")
                plt.plot([-180.0,180.0], np.ones(2)*(beta_contour[I]+alpha_contour[I])*180.0/pi,"--k",color=[.5,.5,.5])
                plt.sca(axList[AxCount2+3])
                plt.plot(beta_contour[I]*180.0/pi, chi_contour[I],"ok",markerFaceColor="None")
                plt.plot([-180.0,180.0], np.ones(2)*(beta_contour[I]+alpha_contour[I])*180.0/pi,"--k")
            AxCount2+=1
    

    
plt.sca(ax11)
plt.ylabel("taper angle")
plt.xlabel("beta")
#plt.box("off")

plt.sca(ax21)
#plt.axis([0.0,1.0,-50,70])

#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("beta")
plt.ylabel("chi")
    
plt.sca(ax31)
plt.axis([-45,25,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("beta")
plt.ylabel("chi")


plt.sca(ax32)
plt.axis([0.0,25.0,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("taper angle")
plt.ylabel("chi")
#iCount+=1

for ax in axList:
        
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


ax12.axes.get_yaxis().set_ticklabels([])
ax13.axes.get_yaxis().set_ticklabels([])


plt.figure(2)
plt.clf()
n = np.int(nLambda/2.0)
plt.axis([-45,105.0,0.0,100.0])
for iPlot in range(2):
    if iPlot==0:
        I = np.arange(0,n+1)
    else:
        I = np.arange(n,nLambda)
        
    
    betas = betas_all[I,0,:]
    alphas_Ref = alphas_Ref_all[I,0,:]
    taper_angles = betas+alphas_Ref
    Lambdas = Lambdas_Ref_all[I,0,:]
        
    
#    plt.contour(betas*180.0/pi,Lambdas,taper_angles*180.0/pi,colors=[[.8,.8,.8,1.0]],levels=np.arange(0.0,60.0,2.0))
    CS = plt.contour(betas*180.0/pi,Lambdas*100.0,taper_angles*180.0/pi,colors="k",levels=np.arange(0.0,60.0,10.0))
    fmt = {}
    if iPlot==0:
#        ax.clabel(CS, CS.levels, inline=True)
        Ilevels = np.arange(0,CS.levels.size)
        for angle in CS.levels:
            if plt.rcParams["text.usetex"]:
                fmt[angle] = r'%.0f \°' % angle
            else:
                fmt[angle] = '%.0f °' % angle
    else:
        tapers_contour = (beta_contour+alpha_contour)*180.0/pi
        chiLevels = np.ones(CS.levels.shape)
        Ilevels = []
        i = 0
        
        for angle in CS.levels:
            chiLevels[i] = chi_contour[np.argmin(abs(tapers_contour-angle))]
            if chiLevels[i]<chi_contour[-1]:
                Ilevels.append(i)
                if plt.rcParams["text.usetex"]:
                    fmt[angle] = r'%.0f \%%' % (chiLevels[i] * 100.0)
                else:
                    fmt[angle] = '%.0f %%' % (chiLevels[i] * 100.0)
            i+=1
    
    textColor = ["r","b"]    
    cLabel = ax.clabel(CS, CS.levels[Ilevels], inline=True, fmt=fmt, fontsize=16)
    for thisText in cLabel:
        thisText._color = textColor[iPlot]
        thisText._fontproperties._family = 'Montserrat'

# set a font dict
font = {'family': 'Montserrat',
'weight': 'bold',
'size': np.int(16)
}
plt.xlabel("$\\beta$ [°]",fontdict=font)
plt.ylabel("$\\lambda$ [%]",fontdict=font)

xTicks = np.arange(-45.0,105.0,15.0)
xTickLabels = []
for i in range(xTicks.size):
    xTickLabels.append("%.0f" % xTicks[i])
plt.xticks(xTicks, xTickLabels)

#matplotlib.rc('font', **font)
#plt.contour()