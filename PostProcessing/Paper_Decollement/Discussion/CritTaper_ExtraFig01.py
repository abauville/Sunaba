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
fig = Figz_Utils.Figure(201,height=11.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00,rightMarginPad = 4.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
graphAxes['12'].axis('off')



## Create taper and get data
# ===========================================
phiRef_list = arr([.5,10.0,20.0,30.0,40.0]) * pi/180.0
LambdaRef_list = arr([.4,.4,.4,.4,.4])
#phiRef_list = arr([40.0,40.0,40.0,40.0,40.0]) * pi/180.0
#LambdaRef_list = arr([.01,.25,.5,.75,.99])
phiRef_list = arr([7.5,15.0,30.0,45.0,
                   7.5,15.0,30.0,45.0]) * pi/180.0
LambdaRef_list = arr([.00001,.51,.773,.869,
                      .00001,.51,.773,.869])

chi_list = [.2,.2,.2,.2,
            1e-4,1e-4,1e-4,1e-4]

Icolor = arr([0,0,0,0,
              1,1,1,1])
     
phiRef_list = arr([15.0,30.0,45.0,
                   15.0,30.0,45.0,
                   15.0,30.0,45.0]) * pi/180.0
LambdaRef_list = arr([.00001,.536,.7321,
                      .00001,.536,.7321,
                      .00001,.536,.7321,])

chi_list = [.99,.99,.99,
            .5,.5,.5,
            1e-4,1e-4,1e-4]

Icolor = arr([0,0,0,
              1,1,1,
              2,2,2])

                  

#phiRef_list = arr([45.0]) * pi/180.0
#LambdaRef_list = arr([.869])

print(np.tan(phiRef_list)*(1.0-LambdaRef_list))

tpr_list = []
nTpr = len(phiRef_list)

rho_w = 1000.0
rho = 2500.0






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


