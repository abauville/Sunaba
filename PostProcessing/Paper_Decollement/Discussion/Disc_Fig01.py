#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 14:39:04 2018

@author: abauville
"""

# Discussion figure about where sand stands in this theory



import numpy as np
from numpy import pi,sin,cos
import matplotlib.pyplot as plt
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '../CriticalTaper')
sys.path.insert(0, '../')
#import CritTaper_dataMaker
from CritTaper_utils import Taper
import CritTaper_Style
import Figz_Utils
from numpy import array as arr
from CritTaper_WedgeVisu import plotWedge


## Create window, Style, etc...
Style = CritTaper_Style.Style()


deg = 180.0/pi

## Figz and axes
# ===========================================
#fig = Figz_Utils.Figure(1,mode="draft",height=11.0)
fig = Figz_Utils.Figure(201,height=15.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00,rightMarginPad = 4.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
graphxPad = graphAxes['info']['xPad']
graphAxes['12'].axis('off')

domainAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.5,topMarginPad=graphH+1.0,rightMarginPad = 4.0+graphW+graphxPad)


## Create taper and get data
# ===========================================

chi_list = arr([1e-4,.15,.3,.45,.6,.75,.9,.99])
chi_list = arr([.5,.15,.3,.45,.6,.75,.9,.99])


phiRef_list = 30.0*1.0/deg*np.ones(chi_list.shape)
LambdaRef_list =0.6*np.ones(chi_list.shape)

Icolor = arr([2,0,0,
              0,0,0,
              0,0,0])

    

print(np.tan(phiRef_list)*(1.0-LambdaRef_list))

tpr_list = []
nTpr = len(phiRef_list)

rho_w = 1000.0
rho = 2500.0

#chi_b = 0.2




for iTpr in range(nTpr):
    print("iTpr = %i" % iTpr)
    phiRef   = phiRef_list[iTpr]#45.0*pi/180.0
    LambdaRef= LambdaRef_list[iTpr]
    chi = chi_list[iTpr]
    LambdaWeak = (1.0-chi) * LambdaRef   + chi
    ## ============= RefTaper =================    
    tpr = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaWeak,
                rho_w=rho_w, rho=rho)
    tpr.computeAlphaVsBeta(n=2010)
    
    tpr_list.append(tpr)

#betaMinRef = np.min(tpr.beta_all)
#betaMaxRef = np.max(tpr.beta_all)



# =============================================================================
#                             Plot alpha vs beta
x0 = -45.0
x1 = 110.0
y0 = -45.0
y1 = 45.0
edgeColor = ["r","r"]
edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
Colors = [[.5,.75,.25],[.25,.5,.75],[.75,.25,.5]]
linestyle_list = ['-','--',':',
                  '-','--',':',
                  '-','--',':']


plt.sca(graphAxes['11'])
for iTpr in range(nTpr):
    tpr = tpr_list[iTpr]    

    chi = chi_list[iTpr]
    
    mueff = np.tan(phiRef_list[iTpr])*(1.0-chi)#*(1.0-LambdaRef_list[iTpr])
    
    beta_Coulomb2 = pi/4.0-(np.arctan(mueff))/2.0
    alpha_Coulomb2  = tpr.findAlpha(beta_Coulomb2,"average")
#    alpha_Coulomb2  = tpr.findAlpha(beta_Coulomb2,"lower")
    
    betaMax = np.max(tpr.beta_all)
    IbetaMax = np.argmax(tpr.beta_all)
    betaMin = np.min(tpr.beta_all)
    IbetaMin = np.argmin(tpr.beta_all)
    
    alpha_IbetaMax = tpr.alpha_all[IbetaMax]
    alpha_IbetaMin = tpr.alpha_all[IbetaMin]
    
    beta_Coulomb = 0.5*(betaMax+betaMin)
    alpha_Coulomb = 0.5*(alpha_IbetaMax+alpha_IbetaMin)
    
    
    
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=0.08,facecolor=Colors[Icolor[iTpr]])
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Colors[Icolor[iTpr]],linestyle=linestyle_list[iTpr],linewidth=0.5)
    plt.plot(beta_Coulomb*deg,alpha_Coulomb*deg,'+',color=Colors[Icolor[iTpr]])
    
    
#    plt.plot(beta_Coulomb2*deg,alpha_Coulomb2*deg,'x',color=Colors[Icolor[iTpr]])
#    
#    plt.fill(tpr.beta_all*deg+(tpr.alpha_all)*deg, (tpr.alpha_all)*deg,alpha=0.08,facecolor=Colors[Icolor[iTpr]])
#    plt.fill(tpr.beta_all*deg+(tpr.alpha_all)*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Colors[Icolor[iTpr]],linestyle='-',linewidth=0.5)
#    plt.plot(beta_Coulomb*deg+alpha_Coulomb*deg,alpha_Coulomb*deg,'+',color=Colors[Icolor[iTpr]])

#    plt.fill(tpr.beta_all*deg, tpr.beta_all*deg+(tpr.alpha_all)*deg,alpha=0.08,facecolor=Colors[Icolor[iTpr]])
#    plt.fill(tpr.beta_all*deg, tpr.beta_all*deg+(tpr.alpha_all)*deg,facecolor="None",edgecolor=Colors[Icolor[iTpr]],linestyle='-',linewidth=0.5)
#    plt.plot(beta_Coulomb*deg, beta_Coulomb*deg+alpha_Coulomb*deg,'+',color=Colors[Icolor[iTpr]])

#                             Plot alpha vs beta
# =============================================================================



# =============================================================================
#                             Plot domains
tprRef = tpr_list[0]
    
for iTpr in range(1,nTpr-1):
    tpr = tpr_list[iTpr]
    
    betaMax = np.max(tpr.beta_all)
    betaMin = np.min(tpr.beta_all)
    n = 100
    betas = np.linspace(betaMin,betaMax,n)
    alpha_up = np.zeros(n)
    Dalpha = np.zeros(n)
    alpha_ref = np.zeros(n)
    iB = 0
    for beta in betas:
        alpha_up[iB] = tpr.findAlpha(beta,"upper")
#        alpha_ref[iB] = tprRef.findAlpha(beta,"average")
        alpha_ref[iB] = tprRef.findAlpha(beta,"lower")
        Dalpha[iB] = (alpha_up[iB]-alpha_ref[iB])#/(alpha_ref[iB])
        iB+=1
    plt.sca(graphAxes['11'])   
    I = np.argmin(np.abs(Dalpha))
    plt.plot(betas[I]*deg,alpha_ref[I]*deg,'ok')
    
    domainIII_left = np.min(tpr.beta_all)
    domainIII_right = betas[I]
    
    domainII_left = domainIII_right
    I = np.argmax(alpha_up)
    domainII_right = betas[I]
    
    plt.plot(betas[I]*deg,alpha_up[I]*deg,'ok')
    
    
    domainI_left = domainII_right
    
    chi = chi_list[iTpr]
    mueff = np.tan(phiRef_list[iTpr])*(1.0-chi)#*(1.0-LambdaRef_list[iTpr])
    beta_Coulomb2 = pi/4.0-(np.arctan(mueff))/2.0
    domainI_right = beta_Coulomb2
    
    plt.sca(domainAxes['11'])
#    plt.plot([domainIII_left,domainIII_right],arr([2.0,2.0])+0.1*iTpr,'b')
#    plt.plot([domainII_left,domainII_right],arr([1.0,1.0])+0.1*iTpr,'y')
#    plt.plot([domainI_left,domainI_right],arr([.0,.0])+0.1*iTpr,'r')
#    
#    plt.plot(arr([domainIII_left,domainIII_right])*deg,arr([iTpr,iTpr]),'b')
#    plt.plot(arr([domainII_left,domainII_right])*deg,arr([iTpr,iTpr]),'y')
#    plt.plot(arr([domainI_left,domainI_right])*deg,arr([iTpr,iTpr]),'r')
    
    plt.plot(arr([domainIII_left,domainIII_right])*deg,arr([chi_list[iTpr],chi_list[iTpr]]),'b')
    plt.plot(arr([domainII_left,domainII_right])*deg,arr([chi_list[iTpr],chi_list[iTpr]]),'y')
    plt.plot(arr([domainI_left,domainI_right])*deg,arr([chi_list[iTpr],chi_list[iTpr]]),'r')
    plt.xlim(-30.0,45.0)
#    plt.plot(betas*deg,Dalpha*deg,'-')

#                             Plot domains
# =============================================================================


#
# =============================================================================
#                                   Axis stuff
            
plt.sca(graphAxes['11'])

#plt.text(x0-(x1-x0)*0.1,y1-(y1-y0)*-0.01,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
#plt.text(x1-(x1-x0)*0.125,y0-(y1-y0)*0.0825,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
Letters = "ABCD"
i = 0

ax = graphAxes['11']        
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
#
#xTickList = np.arange(x0+5.0,x1,20.0)
#yTickList = np.arange(y0+5.0,y1,20.0)
#ax.axes.get_xaxis().set_ticks(xTickList)
#ax.axes.get_yaxis().set_ticks(yTickList)

#
#xTickLabels = []
#for iTick in range(0,len(xTickList)-1):
#    xTickLabels.append("%.f" % xTickList[iTick])
#    
#ax.axes.get_xaxis().set_ticklabels(xTickLabels)
#xTickLabels.append('')


ax.grid(b=True, which='both', color='0.65', linestyle=':')
plt.sca(ax)
#    plt.xlabel("$\\beta$ [°]")

ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s" % Letters[i],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(Style.fontdict['size'])
i+=1
    


#                                   Axis stuff
# =============================================================================


