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


deg = 180.0/pi


#fig = Figz_Utils.Figure(1,mode="draft",height=11.0)
fig = Figz_Utils.Figure(1,height=11.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00,rightMarginPad = 4.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
graphAxes['12'].axis('off')

#drawAxes = Figz_Utils.makeAxes(fig,2,2,aspectRatio=0.47)

drawAxes  = Figz_Utils.makeAxes(fig,1,1,leftMarginPad=1.00+graphW,bottomMarginPad=fig.usableHeight-graphH)
drawW = drawAxes['info']['plotsWidth']
drawH = drawAxes['info']['plotsHeight']
drawAspectRatio = drawH/drawW# 0.295
#
#drawAxes['12'].axis('off')
#drawAxes['21'].axis('off')
#drawAxes['31'].axis('off')


## Create taper and get data
rho_w = 1000.0
rho = 2500.0
phiRef   = 30.0*pi/180.0
LambdaRef=0.6
chi = 0.5

LambdaWeak = (1.0-chi) * LambdaRef   + chi

## ============= RefTaper =================    
tpr = Taper(phi=phiRef, phi_b=phiRef,
            Lambda=LambdaRef, Lambda_b=LambdaWeak,
            rho_w=rho_w, rho=rho)
tpr.computeAlphaVsBeta(n=2010)

betaMinRef = np.min(tpr.beta_all)
betaMaxRef = np.max(tpr.beta_all)



# =============================================================================
# =============================================================================
#                             Plot alpha vs beta
x0 = -15.0
x1 = 90.0
y0 = -25.0
y1 = 25.0
edgeColor = ["r","r"]
edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
plt.sca(graphAxes['11'])
plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=0.08,facecolor=[.5,.75,.25])
plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=[.5,.75,.25],linestyle='-')

beta = 0.0
alpha_up  = tpr.findAlpha(beta,"upper")
alpha_low = tpr.findAlpha(beta,"lower")
alpha_av  = tpr.findAlpha(beta,"average")

plt.plot([beta*deg,beta*deg],[y0,y1],':k',linewidth=0.5)
plt.plot([beta*deg,beta*deg,beta*deg],[alpha_up*deg,alpha_low*deg,alpha_av*deg],'ok',markerFaceColor='None')
plt.axis([x0,x1,y0,y1])

            
plt.sca(graphAxes['11'])

plt.text(x0-(x1-x0)*0.075,y1-(y1-y0)*0.005,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
plt.text(x1-(x1-x0)*0.15,y0-(y1-y0)*0.065,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
Letters = "ABCD"
i = 0

ax = graphAxes['11']        
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

xTickList = np.arange(x0+5.0,x1,10.0)
yTickList = np.arange(y0+5.0,y1,10.0)
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
drawAxes['11'].axis('off')
plt.axis(arr([-.1,1.01,-.01,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'lower',sy0=0.025)
#plt.sca(drawAxes['22'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
pad = 0.025
plotWedge(tpr,'average',origin=arr([0.0,pad+sin(alpha_low)*(2.0-cos(alpha_low))]),plotFaults=False)
#plt.sca(drawAxes['32'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'upper',origin=arr([0.0,2.0*pad+sin(alpha_low)*(2.0-cos(alpha_low))+sin(alpha_av)*(2.0-cos(alpha_av))]),sy0=0.05)

#plt.sca(drawAxes['21'])
#plt.axis([-.1,1.1,-.1,1.1])
#plt.plot([0.0, 1.0], [0.0,0.0])
#plt.plot([0.0, 1.0], [0.0, sin(alpha_)*(2.0-cos(alpha_up))])





























#                           Plot Wedge illustration
# =============================================================================
# =============================================================================