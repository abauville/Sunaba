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
from CritTaper_WedgeVisu import plotArrow

## Create window, Style, etc...
Style = CritTaper_Style.Style()



deg = 180.0/pi


#fig = Figz_Utils.Figure(3,mode="draft",height=16.0)
fig         = Figz_Utils.Figure(3,height=16.0)
graphAxes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.0)
#graphAxes['12'].axis('off')
#graphAxes['13'].axis('off')
graphW      = graphAxes['info']['plotsWidth']
graphH      = graphAxes['info']['plotsHeight']
yPad        = graphAxes['info']['yPad']
graphLmPad  = graphAxes['info']['leftMarginPad']
graphAxes['11'].grid(b=True, which='both', color='0.65', linestyle=':')

#evoAxes     = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,rightMarginPad = 1.0,leftMarginPad=graphW+graphLmPad+1.0)
#evoH        = evoAxes['info']['plotsHeight']
#evoAxes['11'].axis('off')
evoAxes = {}
evoAxes['info'] = graphAxes['info']
evoAxes['11'] = graphAxes['12']
evoAxes['12'] = graphAxes['13']


drawAxes    = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.35,topMarginPad=graphH+1.0)
drawW = drawAxes['info']['plotsWidth']
drawH = drawAxes['info']['plotsHeight']
drawAspectRatio = drawH/drawW
drawAxes['11'].axis('off')
drawAxes['11'].grid(b=True, which='both', color='0.65', linestyle=':')

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
transparency = 0.2
Color  = arr([[.0,.0,.0],[.25,.5,.5],[.25,.25,.75]])
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
plt.text(x0-(x1-x0)*0.125,y1-(y1-y0)*0.050,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
plt.text(x1-(x1-x0)*0.15,y0-(y1-y0)*0.095,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
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
    
yTickLabels = []
for iTick in range(0,len(yTickList)-1):
    yTickLabels.append("%.f" % yTickList[iTick])
    
ax.axes.get_xaxis().set_ticklabels(xTickLabels)
ax.axes.get_yaxis().set_ticklabels(yTickLabels)
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

#plt.sca(evoAxes['12'])
#y0 = 0.0
x0 = 0.0
x1 = 1.0

alphaRef = tpr_list[0].findAlpha(beta,"average")
alpha_list = [tpr_list[0].findAlpha(beta,"average"),
              tpr_list[1].findAlpha(beta,"lower"),
              tpr_list[2].findAlpha(beta,"upper")]
iTpr = 1
enveloppe_list = ['average','lower','upper']
evo_x = arr([.0,    .3  ,  .6 ,  1.0])
letters =     [ 'a' , 'b' , 'c' , 'a' , 'b' , "c'"]
xMod    = arr([ .00 , .0  , .025 ])*(x1-x0)
yMod    = arr([-.03 , .05 , .05  ])*(y1-y0)
for tpr in tpr_list[slice(1,3)]:
    alpha_up  = tpr.findAlpha(beta,"upper")
    alpha_low = tpr.findAlpha(beta,"lower")
    
    alpha = alpha_list[iTpr]
    evo_y = arr([.0,alphaRef,alpha,alpha])*deg
    
#    alpha_av  = tpr.findAlpha(beta,"average")
    plt.sca(evoAxes['1%i' % (iTpr)])
    plt.axis([x0,x1,y0,y1])
    plt.plot([x0,x1],[alphaRef*deg ,alphaRef*deg ],'-',color=Color[0],lineWidth=0.5)
    plt.fill([x0,x1,x1,x0],[alpha_up*deg ,alpha_up*deg, alpha_low*deg, alpha_low*deg  ],lineStyle='None',color=Color[iTpr],alpha=transparency)
    plt.plot([x0,x1],[alpha_up*deg ,alpha_up*deg ],'-',color=Color[iTpr],lineWidth=0.5)
    plt.plot([x0,x1],[alpha_low*deg,alpha_low*deg],'-',color=Color[iTpr],lineWidth=0.5)

    
    plt.plot(evo_x,evo_y,color=Color[iTpr])
    
    for i in range(len(evo_x)-1):
        I = (iTpr-1)*3+i
        plt.text(evo_x[i]+xMod[i],evo_y[i]+yMod[i],letters[I],verticalAlignment='center',horizontalAlignment='center',weight='bold')
        plt.plot(evo_x[i],evo_y[i],'o',color=Color[iTpr])
    
    
    # Add some text
    # ========================
    if iTpr==1:
        plt.text(0.15,alpha_up*deg+0.02*(y1-y0),'extensionally critical taper',verticalAlignment='baseline',fontsize = 8)
        plt.text(0.4,alphaRef*deg-0.02*(y1-y0),'reference critical taper',verticalAlignment='top',fontsize = 8)
        plt.text(0.15,alpha_low*deg-0.02*(y1-y0),'compressionally critical taper',verticalAlignment='top',fontsize = 8)
        
        
        # Delta alpha
        plotArrow([0.72,0.72],arr([alphaRef*deg,alpha_up*deg])+arr([+0.05,-0.05]),0.0,length=2*(alpha_up+alphaRef)/2.0*deg,style='single',headWidth=0.015,headLength=0.9,bodyWidth = 0.004)
        plt.text(0.75,(alpha_up+alphaRef)/2.0*deg,'$\\Delta \\alpha>0$',verticalAlignment='center')
        
        
    elif iTpr==2:
        y = 15
        plt.plot(evo_x[0:2]+arr([.02,-.02]),[y,y],'-',color=[.7,.7,.7])
        plt.plot(evo_x[1:3]+arr([.02,-.02]),[y,y],'-',color=[.7,.7,.7])
        plt.plot(evo_x[2:4]+arr([.02,-.02]),[y,y],'-',color=[.7,.7,.7])
        plt.text(np.mean(evo_x[0:2]),y+0.02*(y1-y0),'build-up',verticalAlignment='baseline',horizontalAlignment='center',fontsize = 8)
        plt.text(np.mean(evo_x[1:3]),y+0.02*(y1-y0),'weakening',verticalAlignment='baseline',horizontalAlignment='center',fontsize = 8)
        plt.text(np.mean(evo_x[2:4]),y+0.02*(y1-y0),'steady-state',verticalAlignment='baseline',horizontalAlignment='center',fontsize = 8)
        
        # Delta alpha
        plotArrow([0.72,0.72],arr([alphaRef*deg,alpha_up*deg])+arr([-0.05, 0.05]),0.0,length=2*(alpha_up+alphaRef)/2.0*deg,style='single',headWidth=0.015,headLength=0.9,bodyWidth = 0.004)
        plt.text(0.75,(alpha_up+alphaRef)/2.0*deg,'$\\Delta \\alpha<0$',verticalAlignment='center')
        
    iTpr+=1
# end iTpr
    

# plot history for moderate chi


#
## plot history for high chi
#alpha = tpr_list[2].findAlpha(beta,"upper")
#plt.plot(arr([.0,    .3  ,  .5 ,  1.0]),
#         arr([.0,alphaRef,alpha,alpha])*deg,color=Color[2])







ax_list = [evoAxes['11']  , evoAxes['12'] ]
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
    
    # x label
    plt.text(x1-(x1-x0)*0.5,y0-(y1-y0)*0.095,"time",rotation=00,fontdict=Style.fontdict,size=12,horizontalAlignment='center')
    

    
    
    
#                              Plot alpha vs time
# =============================================================================
    




#                       ============================






# =============================================================================
#                              Plot wedge illustrations
plt.sca(drawAxes['11'])
ax = drawAxes['11']



axis = arr([-.01,2.11,-.05,2.11])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio])
x0 = axis[0]; x1 = axis[1];
y0 = axis[2]; y1 = axis[3];
plt.axis(axis)
plt.fill([x0,x0,x1,x1],[y0,y1,y1,y0],color=[.95,.96,.95,1.0])
ax.text(x0+0.007*(x1-x0),y1-0.045*(y1-y0),"%s" % Letters[3],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='top',size=12)
Letters = 'CB'

maxH_list = sin(alpha_list)*(2.0-cos(alpha_list))

#origin_list = [arr([0.0,.5*y0+.5*y1     - 0.5*maxH_list[0] ]),
#               arr([1.1,.5*y0+.5*y1    - 0.25*(maxH_list[1]+maxH_list[2]) ]),
#               arr([1.1,.5*y0+.5*y1    + 0.25*(maxH_list[1]+maxH_list[2]) ])]

origin_list = [arr([0.55,0.28 ]),
               arr([0.00,0.0 ]),
               arr([1.10,0.0 ])]

# Styling
colorWedge = arr([[.9,.9,.95,0.5],
                 np.concatenate((Color[1,:],[transparency])),
                 np.concatenate((Color[2,:],[transparency]))])

iTpr = 0
letters = ['a','b','c',"c'"]
for tpr in tpr_list:

    #plotWedge(tpr,'lower',sy0=0.025)
    ##plt.sca(drawAxes['22'])
    ##plt.axis(arr([-.1,1.1,-.1,1.1])*arr([1.0,1.0,drawAspectRatio,drawAspectRatio]))
    pad = 0.025
#    alpha = alpha_av_list[iTpr]
    maxH = maxH_list[iTpr]#sin(alpha_av)*(2.0-cos(alpha_av))
    origin=origin_list[iTpr]
    sx0 = [.9,.9,.92]
    sy0=[.05,.05,.07]
    if iTpr == 0:
        plotWedge(tpr,enveloppe_list[iTpr],plotFaults=True,
                  origin=origin,
                  fx0_list_a = arr([0.5]),
                  fy0_list_a = arr([0.003])*maxH,
                  fx0_list_b = arr([0.2, 0.4, .6, .8 ]),
                  fy0_list_b = arr([0.01]),
                  sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                  colorWedge=colorWedge[iTpr])
    elif iTpr == 1:
#        fl = arr([.32, .5, .787])
        fl = arr([.47, .73])
        plotWedge(tpr,enveloppe_list[iTpr],plotFaults=True,
                  origin=origin,
                  fx0_list_a = fl,
                  fy0_list_a = arr([0.0])*maxH,
                  fx0_list_b = fl,
                  fy0_list_b = arr([0.0])*maxH,
                  sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                  colorWedge=colorWedge[iTpr],
                  colorBase=Color[iTpr],lineWidthBase=2.0)
    elif iTpr == 2:
#        fl = arr([.35,.45,.585,.75])
        fl = arr([.35,.53, .75])
        plotWedge(tpr,enveloppe_list[iTpr],plotFaults=True,
                  origin=origin,
                  fx0_list_a = fl,
                  fy0_list_a = arr([0.0])*maxH,
                  fx0_list_b = fl,
                  fy0_list_b = arr([0.0])*maxH,
                  sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                  colorWedge=colorWedge[iTpr],
                  colorBase=Color[iTpr],lineWidthBase=2.0)
            
    plt.text(origin_list[iTpr][0],origin_list[iTpr][1]+.03,letters[iTpr+1],weight='bold',size=12)
#    ax.text(origin[0],origin[1]+0.05,"%s" % Letters[iTpr],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    iTpr+=1
 
    
# Plot a initial horizontal state
x0 = 0.65
w = .8
y1 = .70
h = .125
plt.fill(x0+arr([.0,w,w,.0]),y1-arr([.0,.0,h,h]),color=colorWedge[0])
plt.plot(x0+arr([.0,w]),y1-arr([.0,.0]),color='k',lineWidth=1.0)
plt.plot(x0+arr([w,.0]),y1-arr([ h, h]),color='k',lineWidth=1.0)

plt.text(x0-.06,y1-h,letters[0],weight='bold',size=12)

# arrow properties
bodyWidth   = 0.005
headLength  = 0.05
headWidth   =headLength/3.0

# plot arrows Top row
x = x0+w/2.0-.15 + arr([.0,.0])
y = y1-h-.01     - arr([.0,.125])
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=bodyWidth,headLength=headLength,headWidth=headWidth,color=Color[1])
plt.text(1.05,np.mean(y),'build-up',verticalAlignment='center',horizontalAlignment='center',fontsize = 12)

x = x0+w/2.0+.15 + arr([.0,.0])
y = y1-h-.01     - arr([.0,.125])
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=bodyWidth,headLength=headLength,headWidth=headWidth,color=Color[2])    
    
# plot arrows Bottom row
x = [origin_list[0][0]+0.5-0.08 , 0.6]
y = [origin_list[0][1]-0.017 , maxH_list[1]]
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=bodyWidth,headLength=headLength,headWidth=headWidth,color=Color[1])
plt.text(1.05,np.mean(y),'weakening',verticalAlignment='center',horizontalAlignment='center',fontsize = 12)

x = [origin_list[0][0]+0.5+0.08 , 1.5]
y = [origin_list[0][1]-0.017 , maxH_list[1]]
plotArrow(x,y,0.0/deg,length=1.0,bodyWidth=bodyWidth,headLength=headLength,headWidth=headWidth,color=Color[2])


#plt.text(np.mean(evo_x[1:3]),y+0.02*(y1-y0),'weakening',verticalAlignment='baseline',horizontalAlignment='center',fontsize = 8)

#                              Plot wedge illustrations
# =============================================================================



## Ti do
#plt.text(0.1,(y0+y1)/2.0,"Add a double arrow \nto indicate $\\Delta \\alpha>0$",color='r',size=15)
# Add "D" for the lower panel
# Try ading the letters to the graphAxes











#plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/CritTaper/Fig03",dpi=300)