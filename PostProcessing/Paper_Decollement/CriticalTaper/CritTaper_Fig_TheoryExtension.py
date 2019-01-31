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
from numpy import pi, tan, sin, cos, arctan, arcsin
from CritTaper_WedgeVisu import plotWedge


## Create window, Style, etc...
Style = CritTaper_Style.Style()


deg = 180.0/pi


#fig = Figz_Utils.Figure(1,mode="draft",height=11.0)
fig = Figz_Utils.Figure(100,height=11.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=.6,rightMarginPad = 0.0)
graphW = graphAxes['info']['plotsWidth']
graphH = graphAxes['info']['plotsHeight']
graphAxes['12'].axis('off')

#drawAxes = Figz_Utils.makeAxes(fig,2,2,aspectRatio=0.47)



## Create taper and get data
rho_w = 1000.0
rho = 2500.0
phiRef   = 30.0*pi/180.0
LambdaRef=0.75
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
x0 = -45.0
x1 = 135.0
y0 = -16.0
y1 = 16.0
edgeColor = ["r","r"]
edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
plt.sca(graphAxes['11'])
plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=Style.alphaBW,facecolor=Style.colorBW)
#plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Style.colorBW,linestyle='-')

I1 = np.argmin(tpr.beta_all)
I0 = np.argmin(tpr.alpha_all)
I1u = np.argmax(tpr.beta_all)
I0u = np.argmax(tpr.alpha_all)
plt.plot(tpr.beta_all[I0:I1]*deg,tpr.alpha_all[I0:I1]*deg,color=[.1,.1,.9],linestyle='-')
plt.plot(tpr.beta_all[I0:I1]*deg,tpr.alpha_all[I0:I1]*deg,color=Style.colorBW,linestyle='-',linewidth=1.0)
plt.plot(tpr.beta_all[I0u:I1u]*deg,tpr.alpha_all[I0u:I1u]*deg,color=Style.colorBW,linestyle='-',linewidth=1.0)
#plt.plot(tpr.beta_all[:I0]*deg,tpr.alpha_all[:I0]*deg,color=Style.colorBW,linestyle='--',linewidth=1.0)
plt.plot(tpr.beta_all[I1:I0u]*deg,tpr.alpha_all[I1:I0u]*deg,color=Style.colorBW,linestyle='--',linewidth=1.0)
plt.plot(tpr.beta_all[I1u::]*deg,tpr.alpha_all[I1u::]*deg,color=Style.colorBW,linestyle='--',linewidth=1.0)

beta = 5.0*np.pi/180.0
alpha_up  = tpr.findAlpha(beta,"upper")
alpha_low = tpr.findAlpha(beta,"lower")
alpha_av  = tpr.findAlpha(beta,"average")

#plt.plot([beta*deg,beta*deg],[y0,y1],':k',linewidth=0.5)
#plt.plot([x0,x1],np.max(tpr.alpha_all)*arr([1,1])*deg,':k',linewidth=0.5)
#plt.plot([beta*deg,beta*deg,beta*deg],[alpha_up*deg,alpha_low*deg,alpha_av*deg],'ok',markerFaceColor='None')
#plt.axis([x0,x1,y0,y1])


# Annotations
#plt.text(beta*deg+3,alpha_up *deg-1,'B',fontdict=Style.fontdict,verticalAlignment='center',size=12)
#plt.text(beta*deg+3,alpha_av *deg,'C',fontdict=Style.fontdict,verticalAlignment='center',size=12)
#plt.text(beta*deg+3,alpha_low*deg,'D',fontdict=Style.fontdict,verticalAlignment='center',size=12)



#plt.text(40,14,'repose angle',rotation=.0,horizontalAlignment='center',size=9)
#plt.text(67,9,'over-critical',rotation=.0,horizontalAlignment='center',size=9)
#plt.text(10,-11,'under-critical',rotation=.0,horizontalAlignment='center',size=9)
#plt.text(61.5,3.5,'extensionally\ncritical',rotation=-58,horizontalAlignment='center',size=9)
#plt.text(35,1,'stable',rotation=00,horizontalAlignment='center',size=9)
#plt.text(20.0,-2.75,'compressively\ncritical',rotation=-55,horizontalAlignment='center',size=9)



            
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
yTickList = np.arange(y0+6.0,y1,10.0)
ax.axes.get_xaxis().set_ticks(xTickList)
ax.axes.get_yaxis().set_ticks(yTickList)


xTickLabels = []
for iTick in range(0,len(xTickList)):
    xTickLabels.append("%.f" % xTickList[iTick])
xTickLabels.append('')

yTickLabels = []
for iTick in range(0,len(yTickList)):
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
#                             Analytical solution


## Compute
n = 1000.0
alpha = np.linspace(-45.0,45.0,n) / deg
rho_w = 1000.0
rho = 2500.0

Lambda = LambdaRef
Lambda_b = LambdaWeak
phi=phiRef
phi_b=phiRef

mu = tan(phi)
mu_b = tan(phi_b)
mu_b_p = mu_b*((1.0-Lambda_b)/(1.0-Lambda))
phi_b_p = arctan(mu_b_p)

alpha_p = arctan( (1.0-rho_w/rho)/(1.0-Lambda)*tan(alpha) )

psi_b = 0.5*arcsin(sin(phi_b_p)/sin(phi)) - 0.5*phi_b_p
psi_0 = 0.5*arcsin(sin(alpha_p)/sin(phi)) - 0.5*alpha_p

taperAngle = psi_b-psi_0

beta =  taperAngle-alpha
beta2= pi-(taperAngle-alpha)-0.5

beta_center  = 0.5* (np.max(beta ) + np.min(beta ))
alpha_center = 0.5* (np.max(alpha) + np.min(alpha))


#plot
plt.plot(beta*deg,alpha*deg,'-k')
I = np.argmin(abs(alpha-0))
beta_at_alpha0 = beta[I]
beta_center = (np.pi/4.0-phi_b_p/2.0)
beta_shift = beta_center-beta_at_alpha0
plt.plot(beta*deg,alpha*deg,'-k')
plt.plot((beta+2.0*beta_shift)*deg,alpha*deg,'--k')

chi = 1.0-(1.0-Lambda_b)/(1.0-Lambda)


plt.plot(beta_center*deg,.0,'ok',markerFaceColor='None')
plt.text(beta_center*deg+2.0,.0-2.0,'$p_{center}$')
hl = 5.0
plt.arrow(beta_at_alpha0*deg,0.0,beta_at_alpha0+2.0*beta_shift*deg-hl,0.0,head_width=1.0,head_length=hl)
plt.text((beta_at_alpha0+beta_center)/2.0*deg,.7,'$\Delta \\beta$',horizontalAlignment='center')
#plt.plot(beta_center*deg,.0,'ok')

#                             Analytical solution
# =============================================================================
# =============================================================================

