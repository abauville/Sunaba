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
from numpy import sin, tan

def plotNeutralWedge(tpr,beta=.0,origin=[.0,.0],enveloppe='lower',colorWedge=arr([.8,.8,.9])):
    plotWedge(tpr,enveloppe,beta=beta,plotFaults=True,plotStress=False,
              origin=origin,
              fx0_list_a = arr([0.5]),
              fy0_list_a = arr([0.003]),
              fx0_list_b = arr([]),
              fy0_list_b = arr([]),
              faultPos = 0.40,
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')


    plotWedge(tpr,enveloppe,beta=beta,plotFaults=True,plotWedge=False,plotFaultsArrow=False,plotStress=False,
              origin=origin,
              fx0_list_a = arr([0.5]),
              fy0_list_a = arr([0.15, 0.3, .45]),
              fx0_list_b = arr([0.3, 0.5, .7 ]),
              fy0_list_b = arr([0.01]),
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')
    
    plotWedge(tpr,enveloppe,beta=beta,plotFaults=True,plotWedge=False,plotStress=False,
              origin=origin,
              fx0_list_a = arr([]),
              fy0_list_a = arr([]),
              fx0_list_b = arr([.9 ]),
              fy0_list_b = arr([0.01]),
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')



def plotStableWedge(tpr,beta,origin=[.0,.0],enveloppe='lower',colorWedge=arr([.8,.8,.9])):
    plotWedge(tpr,enveloppe,beta=beta,plotFaults=False,plotStress=False,plotFaultsArrow=False,
              origin=origin,
              fx0_list_a = arr([0.4]),
              fy0_list_a = arr([0.003]),
              fx0_list_b = arr([.4,.8]),
              fy0_list_b = arr([0.003]),
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')



def plotExtensionalWedge(tpr,beta,origin=[.0,.0],enveloppe='lower',colorWedge=[.7,.7,.8,.5]):
    plotWedge(tpr,enveloppe,beta=beta,plotFaults=True,plotStress=False,plotFaultsArrow=False,
              origin=origin,
              fx0_list_a = arr([0.4]),
              fy0_list_a = arr([0.003]),
              fx0_list_b = arr([.4,.8]),
              fy0_list_b = arr([0.003]),
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')
    plotWedge(tpr,enveloppe,beta=beta,plotFaults=True,plotStress=False,plotWedge=False,
              origin=origin,
              fx0_list_a = arr([.8]),
              fy0_list_a = arr([0.003]),
              fx0_list_b = arr([.6]),
              fy0_list_b = arr([0.003]),
#              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
              faultPos = 0.5,
              f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
              colorWedge=colorWedge,
              colorFaults='k')



## Create window, Style, etc...
Style = CritTaper_Style.Style()



deg = 180.0/pi


#fig = Figz_Utils.Figure(4,mode="draft",height=16.0)
#fig         = Figz_Utils.Figure(4,height=29.7)
#graphAxes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.0)
#DalphaAxes  = Figz_Utils.makeAxes(fig,1,3,aspectRatio=.5,topMarginPad=graphAxes['info']['plotsHeight']+0.5)
#


#drawAxes['1'] = drawAx1['11']

createTapers = True

#graphAxes['12'].axis('off')
#graphAxes['13'].axis('off')
graphW      = graphAxes['info']['plotsWidth']
graphH      = graphAxes['info']['plotsHeight']
yPad        = graphAxes['info']['yPad']
graphLmPad  = graphAxes['info']['leftMarginPad']

if createTapers:
    tpr_dict = {}

for iFinalState in range(2):
#for iFinalState in range(1):
    fig         = Figz_Utils.Figure(4+iFinalState,height=29.7)
    aspectRatio = .8
    graphAxes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=aspectRatio)
    DalphaAxes  = Figz_Utils.makeAxes(fig,1,3,aspectRatio=.5,topMarginPad=graphAxes['info']['plotsHeight']+0.5)
    
    
    plt.sca(graphAxes['12'])
    plt.axis('off')
    plt.sca(graphAxes['13'])
    plt.axis('off')
    
    plt.sca(DalphaAxes['12'])
    plt.axis('off')
    plt.sca(DalphaAxes['13'])
    plt.axis('off')

#    plt.sca(graphAxes['1%i' % (iFinalState+1)])
    plt.sca(graphAxes['11'])
    if createTapers:
        # =============================================================================
        #                       Create taper and get data
        
        rho_w = 1000.0
        rho = 2500.0
        phiRef   = 30.0*pi/180.0
        
        chi = 1e-7
        
        LambdaRef = 0.6
        #chi_list = [.05,0.7,0.7]
#        beta_list = np.linspace(35.0,-5.0,9)/deg
        if iFinalState==0:
            beta_list = arr([0.0,.0,15.0])/deg
            beta_list[0] = 0.0
        else:
            beta_list = arr([0.0,.0,40.0])/deg
            beta_list[0] = 0.0
    
        chi_list = 0.8*np.ones(beta_list.shape)
        chi_list[0] = 1e-6
        #beta_list = [0.0,15.0*1.0/deg,0.0]
        tpr_list = []
        nTpr = len(chi_list)
        iTpr = 0
        for chi in chi_list:
            LambdaWeak = (1.0-chi) * LambdaRef   + chi
            
            if iFinalState==0:
                Lambda=LambdaRef
            else:
                Lambda=LambdaWeak-1e-6
            ## ============= RefTaper =================    
            tpr = Taper(phi=phiRef, phi_b=phiRef,
                        Lambda=Lambda, Lambda_b=LambdaWeak,
                        rho_w=rho_w, rho=rho)
            tpr.computeAlphaVsBeta(n=2010)
            
            tpr_list.append(tpr)
        
                
        #                       Create taper and get data
        # =============================================================================
    
        tpr_dict['%i' % iFinalState] =tpr_list
    #                       ============================
    
    
    # =============================================================================
    #                          Plot tapers alpha vs beta
#    plt.sca(graphAxes['11'])
    x0 = -30.0
    x1 = 50.0
    y0 = -7#-20.0
    y1 = 24.0
    transparency = 0.2
#    Color  = arr([[.25,.5,.5],[.85,.15,.25],[.85,.15,.25]])
#    Color  = arr([[.85,.15,.25],[.25,.5,.5],[.25,.5,.5]])
    Color  = arr([Style.colorRef,Style.colorBW,Style.colorFW,Style.colorBW])
#    Color_w_transparency = arr([[.25,.5,.5,transparency],
#                                [.85,.15,.25,transparency],
#                                [.85,.15,.25,transparency]])
    Color_w_transparency = arr([np.concatenate([Style.colorRef,[Style.alphaBW]]),
                                np.concatenate([Style.colorBW,[Style.alphaBW]]),
                                np.concatenate([Style.colorFW,[Style.alphaBW]]),
                                np.concatenate([Style.colorBW,[Style.alphaBW]])])
    iTpr = 0
    linewidth = 1.0
    
    tpr_list = tpr_dict['%i' % iFinalState]
    for tpr in tpr_list:
    
        #                             Plot alpha vs beta
        

        if iTpr<2:
            if iFinalState==1:
                J = 1
            else:
                J=0
            plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=Style.alphaBW,facecolor=Color[iTpr+J],lineWidth=0)
#            plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Color[iTpr+J],linestyle='-',linewidth=linewidth/2.0)
            n = int(tpr.alpha_all.shape[0]/2)
            I = np.argmin(tpr.beta_all)
            if iTpr==1:
                plt.plot(tpr.beta_all[I:I+n]*deg, (tpr.alpha_all[I:I+n])*deg,color=Color[iTpr+J],linestyle='-',linewidth=linewidth)
            if iTpr==0:
                plt.plot(tpr.beta_all[:I]*deg, (tpr.alpha_all[:I])*deg,color=Color[iTpr],linestyle='-',linewidth=linewidth)
        
        beta = beta_list[iTpr]
        alpha_up  = tpr.findAlpha(beta,"upper")
        alpha_low = tpr.findAlpha(beta,"lower")
        alphaRef_up  = tpr_list[0].findAlpha(beta,"upper")
        alphaRef_low = tpr_list[0].findAlpha(beta,"lower")
        
        if iTpr>0:
            Dalpha = alpha_up-alphaRef_low
            
            if Dalpha<0.0:
                markersize=1.0
                ratio = (x1-x0)/(y1-y0)
                if iFinalState == 0:
                    I = 3
                else:
                    I = 2
                
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/1.0*ratio,alpha_up*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/1.0*1.0/aspectRatio,color=Color[I],lineWidth=0)
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0*1.0/aspectRatio,color=Color[0],lineWidth=0)
#                if iTpr==4 or iTpr == 7:
                plotArrow([beta*deg,beta*deg],arr([alphaRef_low*deg,alpha_up*deg])+arr([+0.05,-0.05]),0.0,length=2*(alpha_up+alphaRef_low)/2.0*deg,style='single',headWidth=1.2,headLength=1.4,bodyWidth = 0.15,color='k')
#                else:
#                    plotArrow([beta*deg,beta*deg],arr([alphaRef_low*deg,alpha_up*deg])+arr([+0.05,-0.05]),0.0,length=2*(alpha_up+alphaRef_low)/2.0*deg,style='single',headWidth=.5,headLength=1.,bodyWidth = 0.02,color=[.7,.7,.75])
            else:
                if iFinalState == 0:
                    I = 3
                else:
                    I = 2
                markersize=1.0
                ratio = (x1-x0)/(y1-y0)
                if iFinalState==0:
                    plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/1.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/1.0,color=Color[I],lineWidth=0)
                else:
                    plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/1.0*ratio,alpha_up*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/1.0*1.0/aspectRatio,color=Color[I],lineWidth=0)
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0*1.0/aspectRatio,color=Color[0],lineWidth=0)
               
        
        if iTpr==1:
            if iFinalState==0:
                betaTemp = np.linspace(-5.0,12.0,30)/deg
            else:
                betaTemp = np.linspace(29.0,31.0,30)/deg
                
            alpha_UpFinTemp = np.zeros(betaTemp.shape)
            alpha_IniTemp = np.zeros(betaTemp.shape)
            iB = 0
            for beta in betaTemp:
                alpha_UpFinTemp[iB] = tpr_list[1].findAlpha(beta,"upper")
                alpha_IniTemp[iB] = tpr_list[0].findAlpha(beta,"average")
                iB+=1
            I = np.argmin(np.abs(alpha_UpFinTemp-alpha_IniTemp))
            alpha_p0 = alpha_IniTemp[I]
            beta_p0 = betaTemp[I]
            
            plt.plot(beta_p0*deg,alpha_p0*deg,'ok',markerFaceColor='None')
            if iFinalState==0:
                plt.text(beta_p0*deg+2.5,alpha_p0*deg+.0,'$p_{cross}$',horizontalAlignment='left',size=12)
            else:
                plt.text(beta_p0*deg+0.5,alpha_p0*deg-2.25,'$p_{cross}$',horizontalAlignment='right',size=12)
            
#            if iFinalState == 0:
#                
#                I = np.argmax(tpr_list[1].alpha_all)
#                beta_p1 = tpr_list[1].beta_all[I]
#                plt.plot(tpr_list[1].beta_all[I]*deg,tpr_list[1].alpha_all[I]*deg,'ok',markerFaceColor='None')
#                plt.text(tpr_list[1].beta_all[I]*deg,tpr_list[1].alpha_all[I]*deg+1.5,'$p_{max}$',horizontalAlignment='center')
#                
            
         

            
        plt.axis([x0,x1,y0,y1])
    
        iTpr+=1
       
    #end iTpr
    plt.text(x0-(x1-x0)*0.125,y1-(y1-y0)*0.050,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
#
#    if iFinalState==0:
#        plt.text(x0-(x1-x0)*0.125,y1-(y1-y0)*0.050,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
#    else:
#        plt.text(x1-(x1-x0)*0.15,y0-(y1-y0)*0.1,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
    Letters = "ABCD"
    
    ax = plt.gca()
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
#    if iFinalState==0:
    xTickList = np.arange(x0+0.0,x1,15.0)
    yTickList = np.arange(y0+2,y1,5.0)
    ax.axes.get_xaxis().set_ticks(xTickList)
    ax.axes.get_yaxis().set_ticks(yTickList)
#    else:
#        xTickList = np.arange(x0+8.0,x1,10.0)
#        yTickList = np.arange(y0,y1,5.0)
#        ax.axes.get_xaxis().set_ticks(xTickList)
#        ax.axes.get_yaxis().set_ticks(yTickList)

    
#    if iFinalState==0:
#        xMod = 0
#    else:
#        xMod = 1
    xTickLabels = []
    for iTick in range(0,len(xTickList)-1):
        xTickLabels.append("%.f" % xTickList[iTick])
    
    
    yTickLabels = []
    for iTick in range(0,len(yTickList)-1):
        yTickLabels.append("%.f" % yTickList[iTick])
    
        
    ax.axes.get_xaxis().set_ticklabels([])
    
    ax.axes.get_yaxis().set_ticklabels(yTickLabels)
#    if iFinalState==0:
#        ax.axes.get_yaxis().set_ticklabels(yTickLabels)
#    else:
#        ax.axes.get_yaxis().set_ticklabels([])
    xTickLabels.append('')
    
    text = ['A', 'B']
    ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),'A',fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(Style.fontdict['size'])
    
         
    #                          Plot tapers alpha vs beta
    # =============================================================================
        
    



#                       ============================
    
    
    


    

    
    # =============================================================================
    #                              Plot wedge illustrations
    

#    tpr_list
        
    drawAxes = {}
    
#    width = 0.15*3
#    height = 0.075*4
    xPad = 0.05
    yPad = 0.05
    
    xCross = 0.4
    width = ( 0.5+xCross/2.0-xPad/2.0 ) * 1.3
    
    yMax_list = []
    tpr_newList = []
    beta_newList = []

    

    beta = beta_list[1]
    beta2 = beta_list[2]
    alpha_Ini   = tpr_list[0].findAlpha(beta,"average")
    alpha_Final = tpr_list[1].findAlpha(beta,"upper")
    
    alpha_Ini2   = tpr_list[0].findAlpha(beta2,"average")
    alpha_Final2 = tpr_list[1].findAlpha(beta2,"upper")

#    yMax_list.append(sin(alpha_Ini))
    beta_newList.append(beta)
    beta_newList.append(beta)
    height = tan(alpha_Ini+beta)*width * 1.05
    yMax_list.append(height/width)
    
    if iFinalState==0:
        yPlot=0.8
        text = 'C'
        suffix = '$_{bw}$'
        yShift = -.075
        
        yPlot2=0.35
        yShift2 = -.175
        text2 = 'D'
    else:
        yPlot=-.025
        text = 'C'
        suffix = '$_{fw}$'
        yShift = -0.075
        
        yPlot2=0.35
        yShift2 = -.175
        text2 = 'D'
    

    plt.text(beta*deg,alpha_Ini*deg+1,text + '$_i$',horizontalAlignment='center') 
    plt.text(beta*deg,alpha_Final*deg-1,text + suffix,verticalAlignment='top',horizontalAlignment='center')     
    
    if iFinalState==0:
        plt.text(beta2*deg+5.5,alpha_Ini2*deg+00,text2 + '$_{i,bw}$',horizontalAlignment='center') 
    else:
        plt.text(beta2*deg,alpha_Ini2*deg+6,text2 + '$_i$',horizontalAlignment='center') 
        plt.text(beta2*deg-6,alpha_Final2*deg-2,text2 + suffix,verticalAlignment='top',horizontalAlignment='center')     
    
    heights = [tan(alpha_Ini+beta)*width * 1.05,
               tan(alpha_Ini2+beta2)*width * 1.05,
               tan(alpha_Final+beta)*width * 1.05,
               tan(alpha_Final2+beta2)*width * 1.05]
#    xPlots = arr([-.5, -.5,  .4,  .4,  .4,  .4])-.5
#    yPlots = arr([ .0,  .6, -.1,  .1,  .4,  .8])-1.0
    xPlots = arr([-.5, -.5, .4,  .4])-.5
    yPlots = arr([ .8,  .0, .8,  .0])-.5
    
#    xPlots = arr([-.5, .4, -.5,  .4])-.5
#    yPlots = arr([ .6,  .6, .0,  .0])-1.0
    
    yPlot=0.8
#    if iFinalState == 0:
    drawAxes['11'                   ] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ xPlots[0],  yPlots[0],  width,heights[0]])['11']
    plt.axis([.0,1.0,.0,heights[0]/width]); plt.axis('off')
    drawAxes['21'                   ] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ xPlots[1],  yPlots[1],  width,heights[1]])['11']
    plt.axis([.0,1.0,.0,heights[1]/width]); plt.axis('off')
        
    drawAxes['12'] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ xPlots[2],  yPlots[2],  width,heights[2]])['11']
    plt.axis([.0,1.0,.0,heights[2]/width]); plt.axis('off')
    drawAxes['22'] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ xPlots[3],  yPlots[3],  width,heights[3]])['11']
    plt.axis([.0,1.0,.0,heights[3]/width]); plt.axis('off')
    
    
    Itpr = [0,1,0,0]
    ItprColor = [0,1,0,1]
    ItprColorFirstGroup = [0,2,0,1]
    
    beta_newList = [beta,beta2]
    # Draw the initial states 1 and 2
    if iFinalState == 0:
        # plot initial states
        plt.sca(drawAxes['11'])            
        plotNeutralWedge(tpr_list[0],beta,enveloppe='upper',colorWedge=Style.colorRef_a)
        
        plt.sca(drawAxes['21'])            
        plotNeutralWedge(tpr_list[0],beta2,enveloppe='upper',colorWedge=Style.colorRef_a)
    else:
        # plot initial states
        plt.sca(drawAxes['11'])            
        plotNeutralWedge(tpr_list[0],beta,enveloppe='upper',colorWedge=Style.colorRef_a)
        
        plt.sca(drawAxes['21'])            
        plotNeutralWedge(tpr_list[0],beta2,enveloppe='upper',colorWedge=Style.colorRef_a)


#    plt.sca(drawAxes['1%i' % (iFinalState+2)])   
    plt.sca(drawAxes['12'])   
    if iFinalState == 0:
        origin=arr([.0,.0])
        tpr = tpr_list[1]
        plotExtensionalWedge(tpr_list[1],beta,enveloppe='upper',colorWedge=Style.colorBW_a)
    else:
        plotNeutralWedge(tpr_list[1],beta,enveloppe='upper',colorWedge=Style.colorFW_a)
        
#    plt.sca(drawAxes['2%i' % (iFinalState+2)])   
    plt.sca(drawAxes['22'])   
    if iFinalState == 0:
        origin=arr([.0,.0])
        tpr = tpr_list[1]
        plotStableWedge(tpr_list[0],beta2,enveloppe='upper',colorWedge=Style.colorBW_a)
    else:
        plotNeutralWedge(tpr_list[1],beta2,enveloppe='upper',colorWedge=Style.colorFW_a)
            
        

    
    #                              Plot wedge illustrations     
    # =============================================================================



    #                       ============================



    
    # =============================================================================
    #                              Plot Dalpha
    
    plt.sca(DalphaAxes['11'])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axes.get_xaxis().set_ticks(xTickList)
#    xTickLabels[-2] = ''
    ax.axes.get_xaxis().set_ticklabels(xTickLabels)

    y1 = 1
    y0 = -1.0
    plt.xlim([x0,x1])
    plt.ylim([y0,y1])
#    plt.ylim([y0,y1])
#    plt.text(x0-(x1-x0)*0.125,y1-(y1-y0)*0.050,"$\\bf{\\bar{\\Delta \\alpha}}$ []",rotation=90,fontdict=Style.fontdict,size=12)
    plt.text(x0-(x1-x0)*0.15,y1-(y1-y0)*0.50,"$\\bf{\\bar{\\Delta \\alpha}}$ []",rotation=90,fontdict=Style.fontdict,size=12,horizontalAlignment='center')
    plt.text(x1-(x1-x0)*0.15,y0-(y1-y0)*0.18,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
    
    
    
    betaDense=np.linspace(np.min(tpr.beta_all)*deg,x1,100)/deg
    alpha_up = np.zeros(100)
    alphaRef_low = np.zeros(100)
    iB = 0
    plt.plot([x0,x1],[0,0],':k',linewidth=.5)
    plt.plot(arr([beta_p0,beta_p0])*deg,[y0,y1],':k',linewidth=.5)
#    plt.plot(arr([beta_p1,beta_p1])*deg,[y0,y1],':k',linewidth=.5)
    for beta in betaDense:
        alpha_up[iB]  = tpr.findAlpha(beta,"upper")
        alphaRef_low[iB] = tpr_list[0].findAlpha(beta,"lower")
        iB+=1
    Dalpha = (alpha_up-alphaRef_low)/(alphaRef_low+betaDense)
    plt.plot(betaDense*deg,Dalpha,'-k')

    ax.text(x0+0.025*(x1-x0),y0+0.05*(y1-y0),'B',fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)

#    maxDalpha = np.max(Dalpha)
#    I = np.argmin(np.abs(Dalpha-maxDalpha))
#    plt.plot(betaDense[I]*deg,Dalpha[I],'ok')

    #                              Plot Dalpha
    # =============================================================================



    

