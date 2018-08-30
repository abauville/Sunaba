#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import CritTaper_Style
import Figz_Utils

Style = CritTaper_Style.Style()

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_Ref_all - alphas_WB_up_all

deg = 180.0/pi

#fig = Figz_Utils.Figure(4,mode="draft",height=10.0)
fig = Figz_Utils.Figure(4,height=10.0)
myAxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.00)


edgeColor = ["r","r"]
edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
linestyle = ["-","-","-"]
i = 0
iCount = 0

colors = np.random.rand(nLambda,4)
colors[:,-1] = 1.0




chiShortList = [0.1, 0.5, .999]

LambdaShortList = np.array([0.01, 0.45, 0.9])
iTapers = []
for Lambda in LambdaShortList:
    iTapers.append(np.argmin(np.abs(LambdaRef_list-Lambda)))


axList = [myAxes['11'], myAxes['12'], myAxes['13']]
AxCount = 0
AxCount2=0
x0 = 30.0-75.0
x1 = 30.0+75.0
y0 = -50.0
y1 = 50.0
for iTaper in iTapers:

    for iSub in range(len(chiShortList)):
        plt.sca(axList[AxCount])
        
        I = np.argmin(abs(chi_list-chiShortList[iSub]))
        i=0
        tpr = Taper_WB[iTaper*nChi+I]

        color = edgeColorWeak[iSub]

        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=0.08,facecolor=edgeColorWeak[iSub])
        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=edgeColorWeak[iSub],linestyle=linestyle[iSub])
           
        plt.axis([x0,x1,y0,y1])

        plt.plot([30.0,30.0],[y0,y1],':k',linewidth=0.5)
        plt.plot([x0,x1],[0.0,0.0],':k',linewidth=0.5)        

        Lambda = LambdaRef_list[iTaper]
        chi = chiShortList[iSub]
    
    tpr = Taper_Ref[iTaper]
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor="r",linestyle=linestyle[iSub])
    plt.sca(axList[AxCount])

    AxCount+=1

            
            
plt.sca(myAxes['11'])
#plt.ylabel("$\\alpha$ [°]")
plt.text(x0-(x1-x0)*0.12,y1-(y1-y0)*0.025,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
#plt.text(x0-(x1-x0)*0.1,y0+(y1-y0)*0.05,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict)
#plt.text(x0-(x1-x0)*0.1,y0-(y1-y0)*0.1,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict)
plt.sca(myAxes['13'])
plt.text(x1-(x1-x0)*0.18,y0-(y1-y0)*0.1,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
Letters = "ABCD"
i = 0
for ax in axList:
        
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    xTickList = np.arange(x0+15.0,x1-14.999,30.0)
    yTickList = np.arange(y0+20.0,y1-14.999,30.0)
    ax.axes.get_xaxis().set_ticks(xTickList)
    ax.axes.get_yaxis().set_ticks(yTickList)
    
    if ax == myAxes['13']:
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
    
myAxes['12'].axes.get_yaxis().set_ticklabels([])
myAxes['13'].axes.get_yaxis().set_ticklabels([])
