#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import CritTaper_dataMaker

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

deg = 180.0/pi

plt.figure(5)
plt.clf()
#plt.subplot(212)
#
#plt.subplot(211)

nLambda = 0
nChi = 0
nBeta = 9

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)


alphas_diff_all = alphas_Ref_all - alphas_WB_up_all


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
#plt.pcolor(betas_all*deg, chis_all, alphas_Diff_all*deg)

#plt.contour(betas_all*deg, chis_all, alphas_Diff_all*deg)


#plt.subplot(212)
iCount = 0
#colors = ["r","g","b","y","m"]
colors = np.random.rand(nLambda,4)
colors[:,-1] = 1.0

ax11 = plt.subplot(131)
ax12 = plt.subplot(132)
ax13 = plt.subplot(133)


chiList = [0.25, 0.99]
#axList = [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33]
axList = [ax11, ax12, ax13]
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all
AxCount = 0
AxCount2=0

for iTaper in range(nLambda):
    
    betas = betas_all[iTaper,:,:]
    alphas_diff = alphas_diff_all[iTaper,:,:]
    alphas_Ref = alphas_Ref_all[iTaper,:,:]
    alphas_WB_up = alphas_WB_up_all[iTaper,:,:]
    alphas_WB_low = alphas_WB_low_all[iTaper,:,:]
    chis = chis_all[iTaper,:,:]
    
    beta_outline = np.concatenate((betas[0,:],betas[1:-2,-1],betas[-1,-1::-1],betas[-2::-1,0]))
    chi_outline = np.concatenate((chis[0,:],chis[1:-2,-1],chis[-1,-1::-1],chis[-2::-1,0]))

    if iTaper in [nLambda-10]:
#        for iSub in range(len(chiList)):
#            plt.sca(axList[AxCount])
#            
#            I = np.argmin(abs(chi_list-chiList[iSub]))
#            i=0
#            for tpr in (Taper_WB[iTaper*nChi+I],Taper_Ref[iTaper]):
#    #            plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,alpha=0.5,facecolor=faceColor[i])
#    #            plt.fill(tpr.beta_all*deg, tpr.alpha_all*deg,facecolor="None",edgecolor=edgeColor[i])
#                
#                if i==1:
#                    color = edgeColor[i]
#                else:
#                    color = edgeColorWeak[iSub]
#                    
#                plt.fill(tpr.beta_all*deg, (tpr.alpha_all+tpr.beta_all)*deg,alpha=0.08,facecolor=color)
#                plt.fill(tpr.beta_all*deg, (tpr.alpha_all+tpr.beta_all)*deg,facecolor="None",edgecolor=color,linestyle=linestyle[iSub])
#                i+=1
#            
#            x0 = 30.0-70.0
#            x1 = 30.0+70.0
#    #        y0 = -45.0
#    #        y1 = 45.0
#            y0 = 0.0
#            y1 = 90.0
#            plt.axis([x0,x1,y0,y1])
#            
#            plt.plot([30.0,30.0],[y0,y1],':k',linewidth=0.5)
#    #        plt.plot([0.0,0.0],[y0,y1],':k',linewidth=0.5)
#            plt.plot([x0,x1],[30.0,30.0],':k',linewidth=0.5)        
#    #        if iSub==1:
#            Lambda = LambdaRef_list[iTaper]
#            chi = chiList[iSub]
#            plt.text(x0+0.025*(x1-x0),y1-0.05*(y1-y0),"$\\lambda=$%i,  $\\chi=$%i " % (int(Lambda*100.0), int(chi*100.0)))
#        
#    
#        plt.sca(axList[AxCount+3])
        
    
        plt.sca(ax11)
        plt.pcolor(betas*deg, chis, alphas_diff*deg,vmin=-20.0,vmax=20.0)
        plt.set_cmap("seismic")
        plt.axis([-40,100,0.0,1.0])
        AxCount+=1
#        plt.colorbar()
#        plt.contour(chis, betas*deg, alphas_diff*deg, [0.0,1e10],colors=[[0.0,0.0,0.0,1.0]])
        
    if iTaper in range(0,nLambda,10):
        plt.sca(ax12)
        CS = plt.contour(betas*deg, chis, alphas_diff*deg, [0.0,1e10],colors=[colors[iTaper,:]])
#        plt.plot([0.0,0.0],[0.0,1.0],"--k",linewidth=1)
            
            
        plt.sca(ax13)
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
        
#        if iTaper in [0,6,9]:
#            
#            for iC in range(len(chiList)):
#                I = np.argmin(abs(chi_contour-chiList[iC]))
#                plt.sca(axList[AxCount2])
#                plt.plot(beta_contour[I]*180.0/pi, (beta_contour[I]+alpha_contour[I])*180.0/pi,"ok",markerFaceColor="None")
#                plt.plot([-180.0,180.0], np.ones(2)*(beta_contour[I]+alpha_contour[I])*180.0/pi,"--k",color=[.5,.5,.5])
#                plt.sca(axList[AxCount2+3])
#                plt.plot(beta_contour[I]*180.0/pi, chi_contour[I],"ok",markerFaceColor="None")
#                plt.plot([-180.0,180.0], np.ones(2)*(beta_contour[I]+alpha_contour[I])*180.0/pi,"--k")
#            AxCount2+=1
    

    

plt.sca(ax11)
plt.axis([-15.0,75.0,0.0,1.0])

#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("$\\beta$ [°]")
plt.ylabel("$\\chi$")
    
plt.sca(ax12)
plt.axis([-45,25,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("$\\beta$  [°]")


plt.sca(ax13)
plt.axis([0.0,25.0,0.0,1.0])
#cbar = plt.colorbar()
#cbar.set_label("alpha Diff")
plt.xlabel("$\\alpha+\\beta$  [°]")
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
