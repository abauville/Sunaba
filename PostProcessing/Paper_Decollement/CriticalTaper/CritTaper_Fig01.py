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
plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=Style.alphaBW,facecolor=Style.colorBW)
plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Style.colorBW,linestyle='-')

beta = 0.0
alpha_up  = tpr.findAlpha(beta,"upper")
alpha_low = tpr.findAlpha(beta,"lower")
alpha_av  = tpr.findAlpha(beta,"average")

plt.plot([beta*deg,beta*deg],[y0,y1],':k',linewidth=0.5)
plt.plot([beta*deg,beta*deg,beta*deg],[alpha_up*deg,alpha_low*deg,alpha_av*deg],'ok',markerFaceColor='None')
plt.axis([x0,x1,y0,y1])


# Annotations
plt.text(beta*deg+3,alpha_up *deg,'B',fontdict=Style.fontdict,verticalAlignment='center',size=12)
plt.text(beta*deg+3,alpha_av *deg,'C',fontdict=Style.fontdict,verticalAlignment='center',size=12)
plt.text(beta*deg+3,alpha_low*deg,'D',fontdict=Style.fontdict,verticalAlignment='center',size=12)

plt.text(55,17.5,'extensionally critical',rotation=-49,horizontalAlignment='center')
plt.text(35,1,'stable',rotation=00,horizontalAlignment='center')
plt.text(23,-1,'compressively critical',rotation=-49,horizontalAlignment='center')



            
plt.sca(graphAxes['11'])
plt.text(x0-(x1-x0)*0.1,y1-(y1-y0)*+0.05,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
plt.text(x1-(x1-x0)*0.125,y0-(y1-y0)*0.0825,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
#plt.yticks([-20,-10,0,10,20],[-20,-10,0,10,''])
Letters = "ABCD"
i = 0

ax = graphAxes['11']        
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

xTickList = np.arange(x0+5.0,x1,20.0)
yTickList = np.arange(y0+5.0,y1,10.0)
ax.axes.get_xaxis().set_ticks(xTickList)
ax.axes.get_yaxis().set_ticks(yTickList)


xTickLabels = []
for iTick in range(0,len(xTickList)):
    xTickLabels.append("%.f" % xTickList[iTick])
xTickLabels.append('')

yTickLabels = []
for iTick in range(0,len(yTickList)-1):
    yTickLabels.append("%.f" % yTickList[iTick])
yTickLabels.append('')    

ax.axes.get_xaxis().set_ticklabels(xTickLabels)
ax.axes.get_yaxis().set_ticklabels(yTickLabels)



#ax.grid(b=True, which='both', color='0.65', linestyle=':')
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

sy0_list = [0.033,0.1,0.075]
pad = 0.025
yOrigin_list = [0.0,pad+sin(alpha_low)*(2.0-cos(alpha_low)),2.0*pad+sin(alpha_low)*(2.0-cos(alpha_low))+sin(alpha_av)*(2.0-cos(alpha_av))]
plt.sca(drawAxes['11'])
drawAxes['11'].axis('off')
plt.axis(arr([-.1,1.01,-.01,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'lower',origin=arr([0.0,yOrigin_list[0]]),sy0=sy0_list[0],
          plotFaultsArrow=False,
          fx0_list_a = arr([0.29,0.4, .55]),
          fy0_list_a = arr([0.00]),
          fx0_list_b = arr([0.29,0.4, .55]),
          fy0_list_b = arr([0.00]),
          faultPos = 0.5,
          colorFaults='k')
plotWedge(tpr,'lower',origin=arr([0.0,yOrigin_list[0]]),sy0=sy0_list[0],
          plotWedge=False,
          plotStress=False,
          fx0_list_a = arr([.75 ]),
          fy0_list_a = arr([0.00]),
          fx0_list_b = arr([.75 ]),
          fy0_list_b = arr([0.00]),
          faultPos = 0.5,
          colorFaults='k')

#plt.sca(drawAxes['22'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))

plotWedge(tpr,'average',origin=arr([0.0,yOrigin_list[1]]),sy0=sy0_list[1],
          plotFaults=False,plotStress=False)




#plt.sca(drawAxes['32'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'upper',origin=arr([0.0,yOrigin_list[2]]),sy0=sy0_list[2],
          plotFaultsArrow=False,
          fx0_list_a = arr([0.15,0.325, .50]),
          fy0_list_a = arr([0.00]),
          fx0_list_b = arr([0.15,0.325, .50]),
          fy0_list_b = arr([0.00]),
          colorFaults='k')
plotWedge(tpr,'upper',origin=arr([0.0,yOrigin_list[2]]),sy0=sy0_list[2],
          plotFaultsArrow=True,plotWedge=False,plotStress=False,
          fx0_list_a = arr([0.675]),
          fy0_list_a = arr([0.00]),
          fx0_list_b = arr([0.675]),
          fy0_list_b = arr([0.00]),
          colorFaults='k')
plt.text(0.91,yOrigin_list[2]+sy0_list[2]-.03,'$\mathbf{\sigma_1}$')

Letters = 'DCB'
for i in range(3):
    plt.text(0.0,yOrigin_list[i]+.02,Letters[i],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)

#plt.sca(drawAxes['21'])
#plt.axis([-.1,1.1,-.1,1.1])
#plt.plot([0.0, 1.0], [0.0,0.0])
#plt.plot([0.0, 1.0], [0.0, sin(alpha_)*(2.0-cos(alpha_up))])





























#                           Plot Wedge illustration
# =============================================================================
# =============================================================================