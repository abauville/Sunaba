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
from numpy import array as arr

Style = CritTaper_Style.Style()

nChi, nBeta, nLambda, LambdaRef_list, chi_list, betas_all, alphas_Ref_all, alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, Taper_Ref, Taper_WB, Taper_WF = CritTaper_dataMaker.getCritTaperFigData(Compute=False)
alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

deg = 180.0/pi

#fig = Figz_Utils.Figure(4,mode="draft",height=9.25)
fig = Figz_Utils.Figure(3,height=9.25)
myAxes = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.00)
#AxesLegend = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.33,
#                                 leftMarginPad = 13.25, rightMarginPad = 0.5)


edgeColor = ["r","r"]
edgeColorWeak = [Style.colorBW*1.5,Style.colorBW,Style.colorBW*.75]

#faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
linestyle = ["-","-","-"]
i = 0
iCount = 0





chiShortList = [0.2, 0.5, .999]
chiExtendedShortList = np.concatenate([ [.0] , chiShortList])
LambdaShortList = np.array([0.4, 0.6, 0.8])
iTapers = []
for Lambda in LambdaShortList:
    iTapers.append(np.argmin(np.abs(LambdaRef_list-Lambda)))


axList = [myAxes['11'], myAxes['12'], myAxes['13']]
AxCount = 0
AxCount2=0
x0 = 30.0-65.0
x1 = 30.0+65.0
y0 = -30.1
y1 = 30.1
iCount = 0
for iTaper in iTapers:

    for iSub in range(len(chiShortList)):
        plt.sca(axList[AxCount])
        
        I = np.argmin(abs(chi_list-chiShortList[iSub]))
        i=0
        tpr = Taper_WB[iTaper*nChi+I]

        color = edgeColorWeak[iSub]

        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=Style.alphaBW,facecolor=edgeColorWeak[iSub])
        plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=edgeColorWeak[iSub],linestyle=linestyle[iSub])
           
        plt.axis([x0,x1,y0,y1])

        if AxCount == 2:
            plt.plot([30.0,30.0],[y0,12],':k',linewidth=0.5)
            plt.plot([30.0,30.0],[20,y1],':k',linewidth=0.5)
        else:
            plt.plot([30.0,30.0],[y0,y1],':k',linewidth=0.5)
        plt.plot([x0,x1],[0.0,0.0],':k',linewidth=0.5)  
#        plt.plot([x0,x1],[30.0,30.0],':k',linewidth=0.5)  
#        plt.plot([x0,x1],[-30.0,-30.0],':k',linewidth=0.5)  

        Lambda = LambdaRef_list[iTaper]
        chi = chiShortList[iSub]
    tpr = Taper_Ref[iTaper]
    plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Style.colorRef,linestyle=linestyle[iSub])
    
    plt.sca(axList[AxCount])
#    plt.text(x0+.7*(x1-x0),y0+0.85*(y1-y0),"$\\lambda = %i$%%" % int(Lambda*100),fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12,weight='normal')
    AxCount+=1

            
            
plt.sca(myAxes['11'])
#plt.ylabel("$\\alpha$ [°]")
plt.text(x0-(x1-x0)*0.12,y1-(y1-y0)*0.25,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
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

    xTickList = np.arange(x0+5.0,x1-4.999,30.0)
    yTickList = np.arange(y0+.1,y1-.0009,30.0)
    ax.axes.get_xaxis().set_ticks(xTickList)
    ax.axes.get_yaxis().set_ticks(yTickList)
    
    if ax == myAxes['13']:
        xTickLabels = []
        for iTick in range(0,len(xTickList)-1):
            xTickLabels.append("%.f" % xTickList[iTick])
            
        ax.axes.get_xaxis().set_ticklabels(xTickLabels)
        xTickLabels.append('')
    
    
#    ax.grid(b=True, which='both', color='0.65', linestyle=':')
    plt.sca(ax)
#    plt.xlabel("$\\beta$ [°]")
    
    ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s. $\mathbf{\lambda=%i}$%%" % (Letters[i],LambdaShortList[i]*100),fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(Style.fontdict['size'])
    i+=1
    

    
myAxes['12'].axes.get_yaxis().set_ticklabels([])
myAxes['13'].axes.get_yaxis().set_ticklabels([])


## Legend
## ====================================
#plt.sca(AxesLegend['11'])
#plt.axis([0,1,0,1])
#plt.axis('off')
#lx0 = .05
#ly0 = .95
#lw = .2
#lh = .2
#lxPad = [.0, .0, .5, .5]
#lyPad = [.0, -.3, .0, -.3]
#
#plt.plot([lx0,lx0+lw],[ly0,ly0],'r')
#for i in range(3):
#    plt.fill(lxPad[i+1]+lx0+arr([0,lw,lw,0]),lyPad[i+1]+ly0+arr([0,0,-lw,-lw]),color=edgeColorWeak[i],alpha=0.08,linewidth=0.0)
#    plt.fill(lxPad[i+1]+lx0+arr([0,lw,lw,0]),lyPad[i+1]+ly0+arr([0,0,-lw,-lw]),edgeColor=edgeColorWeak[i],faceColor='None')
#

note_y0 = [15,20,15,20]
note_x0 = [-34,  2.5, 20.0, 44.0]
note_bar_x = arr([[note_x0[0]+12,note_x0[0]+22],
                  [note_x0[1]+5,note_x0[1]+5],
                  [note_x0[2]+5,note_x0[2]+5],
                  [note_x0[3]+5,note_x0[3]+5]])
    
note_bar_y = arr([[note_y0[0]-1,note_y0[0]-3.5],
                  [note_y0[1]-1,note_y0[1]-8.5],
                  [note_y0[2]-1,note_y0[2]-3.5],
                  [note_y0[3]-1,note_y0[3]-8.5]]) 
plt.text(note_x0[0],note_y0[0],'$\chi=$%i%%' % (chiExtendedShortList[0]*100.0))
for i in range(1,4):
    plt.text(note_x0[i],note_y0[i],'%i%%' % (chiExtendedShortList[i]*100.0))


plt.plot(note_bar_x.T,note_bar_y.T,'-k',linewidth=.5)
