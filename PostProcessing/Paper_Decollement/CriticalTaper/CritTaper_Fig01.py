#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
from numpy import pi, sin,cos,tan
import matplotlib.pyplot as plt
#import CritTaper_dataMaker
from CritTaper_utils import Taper
import CritTaper_Style
import Figz_Utils
from numpy import array as arr

## Create window, Style, etc...
Style = CritTaper_Style.Style()

#nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
#alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

deg = 180.0/pi


fig = Figz_Utils.Figure(4,mode="draft",height=20.0)
graphAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.00)
graphAxes['12'].axis('off')

#drawAxes = Figz_Utils.makeAxes(fig,2,2,aspectRatio=0.47)
drawAspectRatio = 1.00# 0.295
drawAxes = Figz_Utils.makeAxes(fig,1,2,aspectRatio=drawAspectRatio)
drawAxes['11'].axis('off')
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


def intersection(segment1_x, segment1_y, segment2_x, segment2_y):
    p = arr([segment1_x[0],segment1_y[0]])
    r = arr([segment1_x[1]-segment1_x[0],segment1_y[1]-segment1_y[0]])
    
    q = arr([segment2_x[0],segment2_y[0]])
    s = arr([segment2_x[1]-segment2_x[0],segment2_y[1]-segment2_y[0]])

    t = np.cross( (q-p) , (s/np.cross(r,s)) ) 
    

    return [(p+t*r),t]


def plotWedge(taper,enveloppe="within",beta=0.0,
              origin=arr([0.0,0.0]),
              fx0_list = arr([.2, .4, .6, .8]),
              fy0_list = arr([0.0]),plotFaults=True):

    if enveloppe=='lower':
        alpha = taper.findAlpha(beta,"lower")
        psi = taper.psi_bmin
    elif enveloppe=='upper':
        alpha = taper.findAlpha(beta,"upper")
        psi = taper.psi_bmax
    elif enveloppe=='average':
        alpha = taper.findAlpha(beta,"upper")
        alpha += taper.findAlpha(beta,"lower")
        alpha /= 2.0
    else:
        raise ValueError("Enveloppe should be 'upper' or 'lower'")
        
    # Plot wedge outline
    # ===================================
    
    
    surf_x = origin[0] + arr([0.0,1.0])
    surf_y = origin[1] + arr([0.0, sin(alpha)*(2.0-cos(alpha))])
    back_x = origin[0] + arr([1.0,1.0])
    back_y = origin[1] + arr([0.0,sin(alpha)*(2.0-cos(alpha))])
    
    plt.fill(np.concatenate((surf_x, arr([origin[0]+1.0]))),np.concatenate((surf_y, arr([origin[1]]))),color=[.8,.8,.95])
    plt.plot([origin[0], origin[0]+1.0], [origin[1],origin[1]],'-k')
    plt.plot(surf_x,surf_y,'-k')
    
    if plotFaults:
        # Plot stress orientation
        # ===================================
        
        sx0 = origin[0] + 0.9
        sy0 = origin[1] + 0.1
        sl = 0.2
        
        aHL = sl/4.0 # arrow head length
        aHW = sl/12.0
        aBL = sl/2.0-aHL
        aBW = aHW/3.0
        a_x = np.array([-aHL, -aHL, -aHL-aBL, -aHL-aBL, -aHL, -aHL, 0.0]) # arrow points x
        a_y = np.array([+aHW, +aBW,   +aBW  ,   -aBW  , -aBW, -aHW, 0.0]) # arrow points x
        
        rot_a_x = np.cos(psi)*a_x - np.sin(psi)*a_y
        rot_a_y = np.sin(psi)*a_x + np.cos(psi)*a_y
        
        plt.fill(sx0+rot_a_x,sy0+rot_a_y)
        plt.fill(sx0-rot_a_x,sy0-rot_a_y)
        
        
        # Plot faults
        # ===================================
        
    #    fl = 1.0
        
        ca = 30.0*pi/180.0 # Coulomb angle
        fa_list = psi + np.array([+ca, -ca])
        #fa_list = psi + np.array([+ca])
        
        fx0_list = origin[0] + fx0_list
        fy0_list = origin[1] + fy0_list
        
        for fx0 in fx0_list:
            for fy0 in fy0_list:
                for fa in fa_list:
                
                    # check intersection
                    fx = arr([fx0,fx0+cos(fa)])
                    fy = arr([fy0,fy0+sin(fa)])
                    
                    it, t = intersection(fx,fy, surf_x,surf_y)
                    fx = arr([fx0,it[0]])
                    fy = arr([fy0,it[1]])
                    
                    it, t = intersection(fx,fy, back_x,back_y)
                    if (t>0.0 and t<1.0):
                        fx = arr([fx0,it[0]])
                        fy = arr([fy0,it[1]])
                
                    plt.plot( fx , fy , '-r',linewidth=0.75)


plt.sca(drawAxes['12'])
plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'lower')
#plt.sca(drawAxes['22'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
pad = 0.1
plotWedge(tpr,'average',origin=arr([0.0,pad+sin(alpha_low)*(2.0-cos(alpha_low))]),plotFaults=False)
#plt.sca(drawAxes['32'])
#plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
plotWedge(tpr,'upper',origin=arr([0.0,2.0*pad+sin(alpha_low)*(2.0-cos(alpha_low))+sin(alpha_av)*(2.0-cos(alpha_av))]))

#plt.sca(drawAxes['21'])
#plt.axis([-.1,1.1,-.1,1.1])
#plt.plot([0.0, 1.0], [0.0,0.0])
#plt.plot([0.0, 1.0], [0.0, sin(alpha_)*(2.0-cos(alpha_up))])





























#                           Plot Wedge illustration
# =============================================================================
# =============================================================================