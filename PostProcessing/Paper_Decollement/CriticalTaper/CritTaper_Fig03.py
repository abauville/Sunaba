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


fig = Figz_Utils.Figure(3,mode="draft",height=20.0)
#fig         = Figz_Utils.Figure(3,height=20.0)
graphAxes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.0,rightMarginPad = 1.0)
graphW  = graphAxes['info']['plotsWidth']
graphH  = graphAxes['info']['plotsHeight']
yPad    = graphAxes['info']['yPad']
drawAxes    = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.3,rightMarginPad = 1.0,topMarginPad=graphH+0.0*yPad)
drawW = drawAxes['info']['plotsWidth']
drawH = drawAxes['info']['plotsHeight']
drawAspectRatio = drawH/drawW



# =============================================================================
#                       Create taper and get data

rho_w = 1000.0
rho = 2500.0
phiRef   = 30.0*pi/180.0

chi = 1e-7

LambdaRef = 0.7
chi_list = [1e-6,0.3,0.6]
tpr_list = []
nTpr = len(chi_list)
iTpr = 0
for chi in chi_list:
    LambdaWeak = (1.0-chi) * LambdaRef   + chi
    
    ## ============= RefTaper =================    
    tpr = Taper(phi=phiRef, phi_b=phiRef,
                Lambda=LambdaRef, Lambda_b=LambdaWeak,
                rho_w=rho_w, rho=rho)
    tpr.computeAlphaVsBeta(n=2010)
    
    tpr_list.append(tpr)

        
#                       Create taper and get data
# =============================================================================
    




#                       ============================






# =============================================================================
#                          Plot tapers alpha vs beta
plt.sca(graphAxes['11'])
x0 = -18.0
x1 = 30.0
y0 = 0.0#-20.0
y1 = 17.0
transparency = 0.08
Color  = [[.0,.0,.0],[.25,.5,.5],[.25,.25,.75]]
iTpr = 0
for tpr in tpr_list:

    #                             Plot alpha vs beta
    
    
#    edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
#    faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
    plt.sca(graphAxes['11'])
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=transparency,facecolor=Color[iTpr])
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Color[iTpr],linestyle='-')
    
    beta = 0.0
    alpha_up  = tpr.findAlpha(beta,"upper")
    alpha_low = tpr.findAlpha(beta,"lower")
    alpha_av  = tpr.findAlpha(beta,"average")
    
#    plt.plot([beta,x1],[alpha_up*deg,alpha_up*deg],':',color=Color[iTpr],linewidth=0.5)
#    plt.plot([beta,x1],[alpha_low*deg,alpha_low*deg],':',color=Color[iTpr],linewidth=0.5)

    plt.plot([beta*deg,beta*deg],[alpha_up*deg,alpha_low*deg],'o',color=Color[iTpr],markerFaceColor='None')
    plt.axis([x0,x1,y0,y1])
    
    iTpr+=1
    
# end iTpr

plt.plot([beta*deg,beta*deg],[y0,y1],':k',linewidth=0.5)

plt.text(x0-(x1-x0)*0.075,y1-(y1-y0)*0.005,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
plt.text(x1-(x1-x0)*0.125,y0-(y1-y0)*0.065,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
Letters = "ABCD"

ax = graphAxes['11']        
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

xTickList = np.arange(x0+8.0,x1,10.0)
yTickList = np.arange(y0,y1,5.0)
ax.axes.get_xaxis().set_ticks(xTickList)
ax.axes.get_yaxis().set_ticks(yTickList)


xTickLabels = []
for iTick in range(0,len(xTickList)):
    xTickLabels.append("%.f" % xTickList[iTick])
    
ax.axes.get_xaxis().set_ticklabels(xTickLabels)
xTickLabels.append('')

#ax.grid(b=True, which='both', color='0.65', linestyle=':')
plt.sca(ax)
#    plt.xlabel("$\\beta$ [°]")

ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s" % Letters[0],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(Style.fontdict['size'])

     
#                          Plot tapers alpha vs beta
# =============================================================================
    




#                       ============================






# =============================================================================
#                              Plot alpha vs time

plt.sca(graphAxes['12'])
#y0 = 0.0
x0 = 0.0
x1 = 1.0

alphaRef = tpr_list[0].findAlpha(beta,"average")
alpha_list = [tpr_list[0].findAlpha(beta,"average"),
              tpr_list[1].findAlpha(beta,"lower"),
              tpr_list[2].findAlpha(beta,"upper")]
iTpr = 1
enveloppe_list = ['average','lower','upper']
for tpr in tpr_list[slice(1,3)]:
    alpha_up  = tpr.findAlpha(beta,"upper")
    alpha_low = tpr.findAlpha(beta,"lower")
    

    
#    alpha_av  = tpr.findAlpha(beta,"average")
    plt.sca(graphAxes['1%i' % (iTpr+1)])
    plt.axis([x0,x1,y0,y1])
    plt.plot([x0,x1],[alphaRef*deg ,alphaRef*deg ],'-',color=Color[0],lineWidth=0.5)
    plt.fill([x0,x1,x1,x0],[alpha_up*deg ,alpha_up*deg, alpha_low*deg, alpha_low*deg  ],lineStyle='None',color=Color[iTpr],alpha=transparency)
    plt.plot([x0,x1],[alpha_up*deg ,alpha_up*deg ],'-',color=Color[iTpr],lineWidth=0.5)
    plt.plot([x0,x1],[alpha_low*deg,alpha_low*deg],'-',color=Color[iTpr],lineWidth=0.5)

    alpha = alpha_list[iTpr]
    plt.plot(arr([.0,    .2  ,  .4 ,  1.0]),
         arr([.0,alphaRef,alpha,alpha])*deg,color=Color[iTpr])
    
    iTpr+=1
# end iTpr
    

# plot history for moderate chi


#
## plot history for high chi
#alpha = tpr_list[2].findAlpha(beta,"upper")
#plt.plot(arr([.0,    .3  ,  .5 ,  1.0]),
#         arr([.0,alphaRef,alpha,alpha])*deg,color=Color[2])







ax_list = [graphAxes['12']  , graphAxes['13'] ]
i = 1
for ax in ax_list:
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    xTickList = []#np.arange(x0+5.0,x1,10.0)
#    yTickList = []#np.arange(y0+5.0,y1,10.0)
    ax.axes.get_xaxis().set_ticks(xTickList)
    ax.axes.get_yaxis().set_ticks(yTickList)
    
    
    #xTickLabels = []
    #for iTick in range(0,len(xTickList)-1):
    #    xTickLabels.append("%.f" % xTickList[iTick])
    #    
    ax.axes.get_xaxis().set_ticklabels([])
    ax.axes.get_yaxis().set_ticklabels([])
    #xTickLabels.append('')
    
    #ax.grid(b=True, which='both', color='0.65', linestyle=':')
    plt.sca(ax)
    #    plt.xlabel("$\\beta$ [°]")
    
    ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s" % Letters[i],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(Style.fontdict['size'])
    i+=1
#                              Plot alpha vs time
# =============================================================================
    




#                       ============================






# =============================================================================
#                              Plot wedge illustrations
plt.sca(drawAxes['11'])
ax = drawAxes['11']
ax.axis('off')
from CritTaper_WedgeVisu import plotArrow

axis = arr([-.05,2.15,-.05,2.15])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio])
x0 = axis[0]; x1 = axis[1];
y0 = axis[2]; y1 = axis[3];
plt.axis(axis)
Letters = 'CB'

maxH_list = sin(alpha_list)*(2.0-cos(alpha_list))

#origin_list = [arr([0.0,.5*y0+.5*y1     - 0.5*maxH_list[0] ]),
#               arr([1.1,.5*y0+.5*y1    - 0.25*(maxH_list[1]+maxH_list[2]) ]),
#               arr([1.1,.5*y0+.5*y1    + 0.25*(maxH_list[1]+maxH_list[2]) ])]

origin_list = [arr([0.55,0.33 ]),
               arr([0.00,0.0 ]),
               arr([1.10,0.0 ])]

iTpr = 0
for tpr in tpr_list:

    #plotWedge(tpr,'lower',sy0=0.025)
    ##plt.sca(drawAxes['22'])
    ##plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
    pad = 0.025
#    alpha = alpha_av_list[iTpr]
    maxH = maxH_list[iTpr]#sin(alpha_av)*(2.0-cos(alpha_av))
    origin=origin_list[iTpr]
    if iTpr == 0:
        plotWedge(tpr,enveloppe_list[iTpr],plotFaults=True,
                  origin=origin,
                  fx0_list_a = arr([0.5]),
                  fy0_list_a = arr([0.003, 0.33, 0.66])*maxH,
                  fx0_list_b = arr([0.2, 0.4, .6, .8 ]),
                  fy0_list_b = arr([0.01]),
                  sy0 = 0.05)
    else:
        plotWedge(tpr,enveloppe_list[iTpr],plotFaults=True,
                  origin=origin,
                  fx0_list_a = arr([0.25, .5, .75]),
                  fy0_list_a = arr([0.0])*maxH,
                  fx0_list_b = arr([0.25, .5, .75]),
                  fy0_list_b = arr([0.0])*maxH,
                  sy0 = 0.05)
            

#    ax.text(origin[0],origin[1]+0.05,"%s" % Letters[iTpr],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    iTpr+=1
    
    
    
# plot arrows
x = [origin_list[0][0]+0.5-0.05 , 0.5]
y = [origin_list[0][1]-0.02 , maxH_list[1]]
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=0.01,headLength=0.1,headWidth=0.1/3.0,color=Color[1])

x = [origin_list[0][0]+0.5+0.05 , 1.5]
y = [origin_list[0][1]-0.02 , maxH_list[1]]
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=0.01,headLength=0.1,headWidth=0.1/3.0,color=Color[2])
#                              Plot wedge illustrations
# =============================================================================
