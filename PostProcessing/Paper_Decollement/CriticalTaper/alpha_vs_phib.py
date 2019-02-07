#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 14:05:06 2019

@author: abauville

Plot of alpha vs phi b 
or variants
"""
from CritTaper_utils import Taper
import Figz_Utils
import CritTaper_dataMaker
import numpy as np
import matplotlib.pyplot as plt
from numpy import array as arr

deg = np.pi/180.0
fig = Figz_Utils.Figure(100,height=29.0,width=16.0)
graphAxes = Figz_Utils.makeAxes(fig,2,1,aspectRatio=1.0,rightMarginPad = 0.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
#graphAxes['12'].axis('off')

#drawAxes = Figz_Utils.makeAxes(fig,2,2,aspectRatio=0.47)

#(nChi, nBeta, nLambda, LambdaRef_list, 
# chi_list, betas_all, alphas_Ref_all, 
# alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=np.linspace(0.0,30.0,13.0)*np.pi/180.0, nChi=61, nLambda=61,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)

Compute = False
if Compute:
    ## Create taper and get data
    rho_w = 1000.0
    rho = 2500.0
    phiRef   = 30.0*np.pi/180.0
    nLambda = 200
    LambdaRef=np.linspace(.4,1.0,nLambda)
    LambdaRef[-1]=.995
    
    Lambda_hydro = 0.4
    
    Lambda_ov = (LambdaRef-Lambda_hydro)/(1.0-Lambda_hydro)
    
    LambdaWeak = (1.0-chi) * LambdaRef   + chi
    
    ## ============= RefTaper =================    
    nPhi_b = 200
    phi_b = np.linspace(1e-4,phiRef,nPhi_b)
    beta = 0.0
    
    chi = np.linspace(100.0,1.0,nPhi_b)    
    
    alpha_up  = np.zeros([nPhi_b,nLambda])
    alpha_low = np.zeros([nPhi_b,nLambda])
    
    iTpr = 0
    for iLambda in range(nLambda):
        print("iL = %i/%i" % (iLambda,nLambda))
        for iPhi_b in range(nPhi_b):
            this_phi_b = phi_b[iPhi_b]
            Lambda = LambdaRef[iLambda]
            tpr = Taper(phi=phiRef, phi_b=this_phi_b,
                        Lambda=Lambda, Lambda_b=Lambda+1e-6,
                        rho_w=rho_w, rho=rho)
            tpr.computeAlphaVsBeta(n=4010)
            
        #    betaMinRef = np.min(tpr.beta_all)
        #    betaMaxRef = np.max(tpr.beta_all)
            
            alpha_up [iPhi_b,iLambda] = tpr.findAlpha(beta,"upper",tol=1e-3)
            alpha_low[iPhi_b,iLambda] = tpr.findAlpha(beta,"lower",tol=1e-3)
        
    

    
    #   Save stuff
    # ============================================  
    np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/alpha_phi_b.npz",
         alpha_up = alpha_up,
         alpha_low= alpha_low,
         phi_b = phi_b,
         Lambda_Ref = Lambda_Ref,
         phiRef = phiRef)
     
else:
#    #   Load stuff
#    # ============================================  
#    loadedData  = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/alpha_phi_b.npz")
#    alpha_up    = loadedData["alpha_up"][()]
#    alpha_low   = loadedData["alpha_low"][()]
#    phi_b       = loadedData["phi_b"][()]
#    Lambda_Ref  = loadedData["LambdaRef"][()]
#    phiRef  = loadedData["phiRef"][()]
    
    nPhi_b = phi_b.shape[0]
    Lambda_hydro = 0.4
    Lambda_ov = (LambdaRef-Lambda_hydro)/(1.0-Lambda_hydro)
    chi = np.linspace(100.0,1.0,nPhi_b)    
## Plotting  
#alpha_repose = (1.0-Lambda_ov)*np.tan(phiRef)
alpha_repose = np.arctan((1.0-Lambda_ov)*np.tan(phiRef))
alpha_p_repose = alpha_repose/(1.0-Lambda_ov)      
alpha_p_up = alpha_up/(1.0-Lambda_ov)            
alpha_p_low= alpha_low/(1.0-Lambda_ov)
            

#plt.plot(100.0-chi,alpha_up /deg,'b')
#plt.plot(100.0-chi,alpha_up/alpha_repose,'b')
#plt.plot(100.0-chi,alpha_low/deg,'b')
#plt.plot([.0,100.0],arr([alpha_up[-1],alpha_up[-1]])/deg,':b')

#plt.plot(100.0-chi,alpha_p_up/alpha_p_repose*1.0,'k')
#plt.plot(100.0-chi,alpha_p_low/alpha_p_repose*1.0,'k')


#plt.plot(100.0-chi,alpha_p_up/alpha_p_repose*30.0-20.0,'b')
#plt.plot(100.0-chi,alpha_p_low/alpha_p_repose*1.0,'k')

plt.figure(1)
plt.cla()

plotLambda_ov = 0.33
I = np.argmin(np.abs(plotLambda_ov-Lambda_ov))
plotLambda_ov = Lambda_ov[I]
plotAlpha_repose = np.arctan((1.0-plotLambda_ov)*np.tan(phiRef))



#plt.plot((1.0-chi/100.0),(alpha_up[:,I]/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),(alpha_up[:,I]/alpha_up[-1,I]),'k')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'k')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'k')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'b')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'b')

#plt.plot((1.0-chi/100.0),alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'r')

F = alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose)
If = np.argmin(F)
plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
F = (alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose
If = np.argmax(F)
plt.plot([1.0-chi[If]/100.0,1.0-chi[If]/100.0],[-1.0,1.0],'--r')
plt.plot((1.0-chi/100.0),((alpha_up[:,I])/plotAlpha_repose),'k')

F = (alpha_up[:,I])/plotAlpha_repose
If = np.argmax(F)
plt.plot([1.0-chi[If]/100.0,1.0-chi[If]/100.0],[-1.0,1.0],':b')


#plt.plot((1.0-chi/100.0),(alpha_up[:,I]/plotAlpha_repose*(alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')

#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I]-alpha_up[-1,I])/alpha_up[-1,I]),'k')

#plt.plot((1.0-chi/100.0),alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),alpha_up[:,I]/plotAlpha_repose*((alpha_up[:,I]-alpha_up[-1,I])/plotAlpha_repose),'k')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I])/alpha_up[-1,I]),'k')
#plt.plot([1.0-chi[If]/100.0,1.0-chi[If]/100.0],[-1.0,1.0],'--r')

#plt.plot([.0,1.0],[((alpha_up[If,I])/alpha_up[-1,I]),((alpha_up[If,I])/alpha_up[-1,I])],'--r')

#plt.plot((1.0-chi/100.0),((alpha_up[:,I])*180.0/np.pi),'b')
#plt.plot([1.0-chi[If]/100.0,1.0-chi[If]/100.0],[-1.0,30.0],'--r')

#plt.plot((1.0-chi/100.0),alpha_up[:,I]/(plotAlpha_repose-alpha_up[-1,I])*((alpha_up[:,I]-alpha_up[-1,I])/(plotAlpha_repose-alpha_up[-1,I])),'k')

#plt.plot((1.0-chi/100.0),(1.0-chi/100.0)*((alpha_up[:,I])/alpha_up[-1,I]),'r')
#plt.plot((1.0-chi/100.0),(1.0-chi/100.0)*((alpha_up[:,I])-alpha_up[-1,I]),'r')
#plt.plot((1.0-chi/100.0),((alpha_up[:,I])-alpha_up[-1,I])/alpha_up[-1,I],'b')


plt.axis([.0,1.0,-1.0,1.0])

plt.sca(graphAxes['11'])
plt.cla()
#plt.contourf(Lambda_ov*100.0,chi,((alpha_up-alpha_up[-1,:])/(alpha_repose)),np.linspace(-1.0,1.0,20))
plt.contourf(Lambda_ov*100.0,chi,(alpha_up/alpha_repose[0]),np.linspace(.0,1.0,20))
plt.contour(Lambda_ov*100.0,chi,(alpha_up/alpha_repose[0]),[-1.0, .1, .5])

plt.contourf(Lambda_ov*100.0,chi,(alpha_up)*180.0/np.pi,np.linspace(.0,30.0,20))

plt.colorbar()
plt.ylim([100.0,.0])
plt.set_cmap('seismic')

plt.sca(graphAxes['21'])
plt.cla()

chi2D,dum = np.meshgrid(chi,chi)

##plt.plot(100.0-chi,(alpha_p_up-alpha_p_up[-1]) /deg,'k')
##plt.plot(100.0-chi,(alpha_p_up-alpha_p_up[-1])/alpha_p_repose,'k')
##plt.plot(chi,(alpha_p_up-alpha_p_up[-1])/(alpha_p_repose-alpha_p_up[-1]),'k')
#plt.plot(phi_b,(alpha_p_up-alpha_p_up[-1])/(alpha_p_repose),'k')
#
#
##plt.plot(100.0-percent,(alpha_p_up-alpha_p_up[-1]) /deg,'k')
##plt.plot(100.0-percent,alpha_p_low/deg,'k')
#
##plt.plot(100.0-percent,alpha_p_up /deg,'k')
##plt.plot(100.0-percent,alpha_p_low/deg,'k')
##plt.plot([.0,100.0],arr([alpha_p_up[-1],alpha_p_up[-1]])/deg,':k')
#
#
##plt.plot(phi_b/deg,alpha_up/deg)
##plt.plot(phi_b/deg,alpha_low/deg)
#plt.xlabel('$\\chi$ [°]')
##plt.xlabel('$\\phi_b$ [°]')
#plt.ylabel("$\\alpha'$ [°]")
#plt.ylim([-1.0,1.0])
#    
#
#plt.sca(graphAxes['11'])
#plt.cla()
#
##plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1]),20)
##plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up-alpha_p_up[-1]),40)
#plt.contourf(Lambda_ov*100.0,chi,((alpha_p_up)/(alpha_p_repose)),np.linspace(0.98,1.0,100))
#plt.contourf(Lambda_ov*100.0,100.0-chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose)),arr([-1,.0,.5,1.0]))
#plt.contourf(Lambda_ov*100.0,100.0-chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose)),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,100.0-chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose)),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,100.0-chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose-alpha_p_up[-1])),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,np.log10(np.abs(((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose-alpha_p_up[-1])))),np.linspace(-5,5.0,100))

#plt.contourf(Lambda_ov*100.0,chi,(((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose-alpha_p_up[-1]))),np.linspace(.0,1.0,100))

#plt.contourf(Lambda_ov*100.0,chi,alpha_up*(alpha_p_up-alpha_p_up[-1,:])/(alpha_p_up[-1]),np.linspace(-.1,.1,20))
#plt.contourf(Lambda_ov*100.0,chi,alpha_up*(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-.1,.1,20))
#plt.contourf(Lambda_ov*100.0,chi,alpha_up/(alpha_p_up[-1])*(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-.6,.6,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,alpha_up/(alpha_p_up[-1])*(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-.6,.6,20))
#plt.contourf(Lambda_ov*100.0,chi,alpha_up/(alpha_p_up[-1])*(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-.6,.6,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-1.0,1.0,20))

#plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up-alpha_p_up[-1,:]),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(1.0-chi2D.T/100.0)*(alpha_p_up-alpha_p_up[-1,:]),np.linspace(-.1,.1,200))
#plt.contourf(Lambda_ov*100.0,chi,(1.0-chi2D.T/100.0)*(alpha_up-alpha_up[-1,:]),np.linspace(-.2,.2,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up-alpha_p_up[-1,:]),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up-alpha_p_up[-1,:])/(alpha_p_repose),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_repose),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,:]),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_repose)*(1.0-(-alpha_repose+alpha_up[-1,:])/(alpha_repose)),np.linspace(-1.0,1.0,20))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,0]),np.linspace(-1.0,1.0,100))
#plt.contour(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,0]),[-2.0, -.5, -.25])
#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,0]),np.linspace(-1.0,1.0,20))

#plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,:])*(alpha_up)/(alpha_repose[0]),np.linspace(-1.0,1.0,20))
plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1,:])/(alpha_up[-1,:])*(alpha_up)/(alpha_up[-1,0]),np.linspace(-1.0,1.0,50))
#plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up),np.linspace(-1,1.0,20))

plt.colorbar()
plt.ylim([100.0,.0])
plt.set_cmap('seismic')
#
#plt.sca(graphAxes['21'])
#plt.cla()
##plt.contourf(Lambda_ov*100.0,chi,(alpha_up-alpha_up[-1]),20)
##plt.contourf(Lambda_ov*100.0,chi,(alpha_p_up-alpha_p_up[-1]),20)
#plt.contourf(Lambda_ov*100.0,chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose)),np.linspace(-1.0,1.0,30))
#plt.plot([33.0,33.0],[.0,100.0],'--',linewidth=0.5)
#plt.plot([66.0,66.0],[.0,100.0],'--',linewidth=0.5)
##plt.contourf(Lambda_ov*100.0,chi,((alpha_p_up-alpha_p_up[-1])/(alpha_p_repose-alpha_p_up[-1])),np.linspace(.98,1.0,200))
#plt.colorbar()
#plt.ylim([100.0,.0])
#plt.set_cmap('seismic')
#




