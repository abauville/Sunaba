#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""

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

#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

deg = 180.0/pi


#fig = Figz_Utils.Figure(2,mode="draft",height=11.0)
fig = Figz_Utils.Figure(2,height=10.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00,rightMarginPad = 4.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
graphAxes['12'].axis('off')

#drawAxes = Figz_Utils.makeAxes(fig,2,2,aspectRatio=0.47)

drawAxes  = Figz_Utils.makeAxes(fig,1,1,leftMarginPad=1.00+graphW,bottomMarginPad=fig.usableHeight-graphH)
drawW = drawAxes['info']['plotsWidth']
drawH = drawAxes['info']['plotsHeight']
drawAspectRatio = drawH/drawW# 0.295


#drawW = drawAxes['info']['plotsWidth']
#drawH = drawAxes['info']['plotsHeight']
#drawAspectRatio = drawH/drawW# 0.295

#
#drawAxes['12'].axis('off')
#drawAxes['21'].axis('off')
#drawAxes['31'].axis('off')


## Create taper and get data
rho_w = 1000.0
rho = 2500.0
phiRef   = 30.0*pi/180.0

chi = 1e-6

LambdaRef_list = np.arange(0.0,1.001,0.1)
LambdaRef_list[00] += 1e-6
LambdaRef_list[-1] -= 1e-6
LambdaRef_selectList = [0.9,0.0]

tpr_list = []
nTpr = len(LambdaRef_list)
LambdaRef_selectList_I = arr([0,0])
tpr_selectList = [0,0]
i = 0
for LambdaRef in LambdaRef_selectList:
    LambdaRef_selectList_I[i] = np.argmin(abs(LambdaRef_list-LambdaRef))
    i+=1
    
iTpr = 0
for LambdaRef in LambdaRef_list:
    LambdaWeak = (1.0-chi) * LambdaRef   + chi
    
    ## ============= RefTaper =================    
    tpr = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaWeak,
                rho_w=rho_w, rho=rho)
    tpr.computeAlphaVsBeta(n=2010)
    
    betaMinRef = np.min(tpr.beta_all)
    betaMaxRef = np.max(tpr.beta_all)
    tpr_list.append(tpr)
    if iTpr in LambdaRef_selectList_I:
        I = np.argmin(abs(LambdaRef_selectList_I-iTpr))
        tpr_selectList[I] = tpr
    iTpr+=1


# =============================================================================
# =============================================================================
#                             Plot alpha vs beta
iTpr = 0    
alpha_av_list = np.zeros(len(LambdaRef_selectList))
for tpr in tpr_list:
    x0 = -45.0
    x1 = 115.0
    y0 = -45.0
    y1 = 45.0
    edgeColor = ["r","r"]
    edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
    faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
    plt.sca(graphAxes['11'])
#    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=0.08,facecolor=[.5,.75,.25])
    
    if iTpr in LambdaRef_selectList_I:
        edgeColor = [.5,.75,.25]
#        edgeColor = [.8,.8,.9]
        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=edgeColor,linestyle='-',linewidth=1.0)
    else:
        edgeColor = [.8,.8,.9]
        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=edgeColor,linestyle='-',linewidth=.5)
    
    beta = 0.0
    alpha_av  = tpr.findAlpha(beta,"average")
    
    if iTpr in LambdaRef_selectList_I:
#        plt.plot([beta*deg,beta*deg],[y0,y1],':k',linewidth=0.5)
        plt.plot([beta*deg],[alpha_av*deg],'ok',markerFaceColor='None')
        I = np.argmin(abs(LambdaRef_selectList_I-iTpr))
        alpha_av_list[I] = alpha_av

    
    if not iTpr%2:
#        betaMin = np.min(tpr.beta_all)
#        IBetaMin = np.argmin(tpr.beta_all)
#        plt.text(betaMin*deg-2.0,tpr.alpha_all[IBetaMin]*deg,"$%.1f$" % LambdaRef_list[iTpr],horizontalAlignment = 'right',verticalAlignment='center')
        
        betaMax = np.max(tpr.beta_all)
        IBetaMax = np.argmax(tpr.beta_all)
        if iTpr==nTpr-1:
            plt.text(betaMax*deg+2.0,tpr.alpha_all[IBetaMax]*deg,"$\\lambda=%.1f$" % LambdaRef_list[iTpr],horizontalAlignment = 'left',verticalAlignment='center')
        else:
            plt.text(betaMax*deg+2.0,tpr.alpha_all[IBetaMax]*deg,"$%.1f$" % LambdaRef_list[iTpr],horizontalAlignment = 'left',verticalAlignment='center')
    
    iTpr+=1
    
plt.axis([x0,x1,y0,y1])
    
                
plt.sca(graphAxes['11'])

plt.text(x0-(x1-x0)*0.095,y1-(y1-y0)*0.025,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
plt.text(x1-(x1-x0)*0.175,y0-(y1-y0)*0.08,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
Letters = "ABCD"
i = 0
    
ax = graphAxes['11']        
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

xTickList = np.arange(x0+15.0,x1,30.0)
yTickList = np.arange(y0+15.0,y1,30.0)
ax.axes.get_xaxis().set_ticks(xTickList)
ax.axes.get_yaxis().set_ticks(yTickList)


xTickLabels = []
for iTick in range(0,len(xTickList)-1):
    xTickLabels.append("%.f" % xTickList[iTick])
    
ax.axes.get_xaxis().set_ticklabels(xTickLabels)
xTickLabels.append('')


ax.grid(b=True, which='both', color='0.65', linestyle=':')
plt.sca(ax)
#    plt.xlabel("$\\beta$ [°]")

ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s" % Letters[i],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(Style.fontdict['size'])
i+=1
    
#                             Plot alpha vs beta
# =============================================================================
# =============================================================================






# =============================================================================
# =============================================================================
#                           Plot Wedge illustration


plt.sca(drawAxes['11'])
ax = drawAxes['11']
ax.axis('off')
iTpr = 0
axis = arr([-.1,1.01,-.01,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio])
x0 = axis[0]; x1 = axis[1];
y0 = axis[2]; y1 = axis[3];
plt.axis(axis)
Letters = 'CB'
for tpr in tpr_selectList:

    #plotWedge(tpr,'lower',sy0=0.025)
    ##plt.sca(drawAxes['22'])
    ##plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
    pad = 0.025
    alpha_av = alpha_av_list[iTpr]
    maxH = sin(alpha_av)*(2.0-cos(alpha_av))
    origin=arr([0.0,iTpr*(2.0*pad+sin(alpha_av_list[iTpr-1])*(2.0-cos(alpha_av_list[iTpr-1])))])
    plotWedge(tpr,'lower',plotFaults=True,
              origin=origin,
              fx0_list_a = arr([0.5]),
              fy0_list_a = arr([0.003, 0.33, 0.66])*maxH,
              fx0_list_b = arr([0.2, 0.4, .6, .8 ]),
              fy0_list_b = arr([0.01]),
              sy0 = 0.05)

    ax.text(origin[0],origin[1]+0.05,"%s" % Letters[iTpr],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    iTpr+=1
##plt.sca(drawAxes['32'])
##plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
#plotWedge(tpr,'upper',origin=arr([0.0,2.0*pad+sin(alpha_low)*(2.0-cos(alpha_low))+sin(alpha_av)*(2.0-cos(alpha_av))]),sy0=0.05)

#plt.sca(drawAxes['21'])
#plt.axis([-.1,1.1,-.1,1.1])
#plt.plot([0.0, 1.0], [0.0,0.0])
#plt.plot([0.0, 1.0], [0.0, sin(alpha_)*(2.0-cos(alpha_up))])



#                           Plot Wedge illustration
# =============================================================================
# =============================================================================
    
    
    
    
    
    
    
#    
#    
#
## =============================================================================
## =============================================================================
##                             Plot alpha vs lambda
#behindSmallAxes = Figz_Utils.makeAxes(fig,1,1,leftMarginPad=1.00+graphW*0.55,bottomMarginPad=fig.usableHeight-graphH*.48,aspectRatio=1.50,setAspectRatioBasedOn='y')
#smallAxes       = Figz_Utils.makeAxes(fig,1,1,leftMarginPad=1.00+graphW*0.75,bottomMarginPad=fig.usableHeight-graphH*.4,aspectRatio=1.50,setAspectRatioBasedOn='y')
#    
#n = 10
#LambdaRef_list = np.linspace(1e-6,1.0-1e-6,n)
#tpr_list = []
#nTpr = len(LambdaRef_list)
#alphas = np.zeros(n)
#iL = 0
#for LambdaRef in LambdaRef_list:
#    LambdaWeak = (1.0-chi) * LambdaRef   + chi
#    
#    ## ============= RefTaper =================    
#    tpr = Taper(phi=phiRef, phi_b=phiRef,
#                Lambda=LambdaRef, Lambda_b=LambdaWeak,
#                rho_w=rho_w, rho=rho)
#    tpr.computeAlphaVsBeta(n=2010)
#    
#    betaMinRef = np.min(tpr.beta_all)
#    betaMaxRef = np.max(tpr.beta_all)
#    tpr_list.append(tpr)
#    alphas[iL] = tpr.findAlpha(beta,"average")
#    iL+=1
#
#ax = smallAxes['11']
#plt.sca(ax)
#plt.axis([0.0,1.0,0.0,25.0])
#plt.plot(LambdaRef_list,alphas*deg,'-k')
## Only show ticks on the left and bottom spines
## Hide the right and top spines
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#
## Hide the spines of the behind graph
## Thes axes are used to hide the minor grid of the graph axes
#ax = behindSmallAxes['11']
#ax.spines['left'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#ax.yaxis.set_ticks([])
#ax.xaxis.set_ticks([])
#
##                             Plot alpha vs lambda
## =============================================================================
## =============================================================================
