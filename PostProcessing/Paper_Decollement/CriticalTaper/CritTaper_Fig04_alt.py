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

## Create window, Style, etc...
Style = CritTaper_Style.Style()



deg = 180.0/pi


#fig = Figz_Utils.Figure(4,mode="draft",height=16.0)
fig         = Figz_Utils.Figure(4,height=29.7)
graphAxes   = Figz_Utils.makeAxes(fig,1,3,aspectRatio=1.0)

plt.sca(graphAxes['13'])
plt.axis('off')

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
    plt.sca(graphAxes['1%i' % (iFinalState+1)])
    if createTapers:
        # =============================================================================
        #                       Create taper and get data
        
        rho_w = 1000.0
        rho = 2500.0
        phiRef   = 30.0*pi/180.0
        
        chi = 1e-7
        
        LambdaRef = 0.7
        #chi_list = [.05,0.7,0.7]
        beta_list = np.linspace(35.0,-5.0,9)/deg
        beta_list[0] = 0.0
        
    
        chi_list = 0.7*np.ones(beta_list.shape)
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
    x0 = -18.0
    x1 = 35.0
    y0 = -7#-20.0
    y1 = 18.0
    transparency = 0.2
    Color  = arr([[.25,.5,.5],[.85,.15,.25],[.85,.15,.25]])
    Color_w_transparency = arr([[.25,.5,.5,transparency],
                                [.85,.15,.25,transparency],
                                [.85,.15,.25,transparency]])
    iTpr = 0
    linewidth = 1.0
    
    tpr_list = tpr_dict['%i' % iFinalState]
    for tpr in tpr_list:
    
        #                             Plot alpha vs beta
        
        
    #    edgeColorWeak = [[.5,.75,.25],[.25,.5,.5],[.25,.25,.75]]
    #    faceColor = [np.array([202,231,202])/255,[0,0,0],[0,0,0]]
    
#        plt.sca(graphAxes['11'])
        if iTpr<2:
            plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,alpha=transparency,facecolor=Color[iTpr],linewidth=0)
            plt.fill(tpr.beta_all*deg, (tpr.alpha_all)*deg,facecolor="None",edgecolor=Color[iTpr],linestyle='-',linewidth=linewidth/2.0)
            n = int(tpr.alpha_all.shape[0]/2)
            I = np.argmin(tpr.beta_all)
            if iTpr==1:
                plt.plot(tpr.beta_all[I:I+n]*deg, (tpr.alpha_all[I:I+n])*deg,color=Color[iTpr],linestyle='-',linewidth=linewidth)
    #            plt.plot(tpr.beta_all[:I]*deg, (tpr.alpha_all[:I])*deg,color=Color[iTpr],linestyle='-',linewidth=linewidth/2.0)
            if iTpr==0:
                plt.plot(tpr.beta_all[:I]*deg, (tpr.alpha_all[:I])*deg,color=Color[iTpr],linestyle='-',linewidth=linewidth)
    #            plt.plot(tpr.beta_all[n::]*deg, (tpr.alpha_all[n::])*deg,color=Color[iTpr],linestyle='-',linewidth=linewidth/2.0)
        
        beta = beta_list[iTpr]
        alpha_up  = tpr.findAlpha(beta,"upper")
        alpha_low = tpr.findAlpha(beta,"lower")
        alphaRef_up  = tpr_list[0].findAlpha(beta,"upper")
        alphaRef_low = tpr_list[0].findAlpha(beta,"lower")
        
        if iTpr>0:
            Dalpha = alpha_up-alphaRef_low
            
        #    plt.plot([beta,x1],[alpha_up*deg,alpha_up*deg],':',color=Color[iTpr],linewidth=0.5)
        #    plt.plot([beta,x1],[alpha_low*deg,alpha_low*deg],':',color=Color[iTpr],linewidth=0.5)
            if Dalpha<0.0:
                markersize=1.0
                ratio = (x1-x0)/(y1-y0)
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alpha_up*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0,color=Color[2],lineStyle='None')
                markersize=.5
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0,color=Color[0],lineStyle='None')
                plotArrow([beta*deg,beta*deg],arr([alphaRef_low*deg,alpha_up*deg])+arr([+0.05,-0.05]),0.0,length=2*(alpha_up+alphaRef_low)/2.0*deg,style='single',headWidth=.5,headLength=1.,bodyWidth = 0.025)
            else:
                markersize=1.0
                ratio = (x1-x0)/(y1-y0)
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0,color=Color[2],lineStyle='None')
                markersize=.5
                plt.fill(beta*deg+arr([-1.0,1.0,1.0,-1.0])*markersize/2.0*ratio,alphaRef_low*deg+arr([1.0,1.0,-1.0,-1.0])*markersize/2.0,color=Color[0],lineStyle='None')
        #        plt.plot([beta*deg],[alphaRef_low*deg],'s',color=Color[iTpr],markeredgewidth=0.0,markersize=8.0)
        #        plt.plot([beta*deg],[alphaRef_low*deg],'s',color=Color[0],markeredgewidth=0.0,markersize=4.0)
        
        if iTpr==1:
            if iFinalState==0:
                betaTemp = np.linspace(5.0,12.0,30)/deg
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
            plt.text(beta_p0*deg+1.75,alpha_p0*deg+.0,'$p_{cross}$',horizontalAlignment='left')
            
            if iFinalState == 0:
                I = np.argmax(tpr_list[1].alpha_all)
                plt.plot(tpr_list[1].beta_all[I]*deg,tpr_list[1].alpha_all[I]*deg,'ok',markerFaceColor='None')
                plt.text(tpr_list[1].beta_all[I]*deg,tpr_list[1].alpha_all[I]*deg+1.5,'$p_{max}$',horizontalAlignment='center')
                
            
            
            
        plt.axis([x0,x1,y0,y1])
    
        iTpr+=1
       
    if iFinalState==0:
        plt.text(x0-(x1-x0)*0.125,y1-(y1-y0)*0.050,"$\\bf \\alpha$ [°]",rotation=90,fontdict=Style.fontdict,size=12)
    else:
        plt.text(x1-(x1-x0)*0.15,y0-(y1-y0)*0.1,"$\\bf \\beta$ [°]",rotation=00,fontdict=Style.fontdict,size=12)
    Letters = "ABCD"
    
    ax = plt.gca()
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    if iFinalState==0:
        xTickList = np.arange(x0+8.0,x1,10.0)
        yTickList = np.arange(y0,y1,5.0)
        ax.axes.get_xaxis().set_ticks(xTickList)
        ax.axes.get_yaxis().set_ticks(yTickList)
    else:
        xTickList = np.arange(x0+8.0,x1,10.0)
        yTickList = np.arange(y0,y1,5.0)
        ax.axes.get_xaxis().set_ticks(xTickList)
        ax.axes.get_yaxis().set_ticks(yTickList)

    
    if iFinalState==0:
        xMod = 0
    else:
        xMod = 1
    xTickLabels = []
    for iTick in range(0,len(xTickList)-xMod):
        xTickLabels.append("%.f" % xTickList[iTick])
    
    
    yTickLabels = []
    for iTick in range(0,len(yTickList)):
        yTickLabels.append("%.f" % yTickList[iTick])
            
        
    ax.axes.get_xaxis().set_ticklabels(xTickLabels)
    if iFinalState==0:
        ax.axes.get_yaxis().set_ticklabels(yTickLabels)
    else:
        ax.axes.get_yaxis().set_ticklabels([])
    xTickLabels.append('')
    
    ax.text(x0+0.025*(x1-x0),y0+0.025*(y1-y0),"%s" % Letters[iFinalState],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
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
    width = 0.5+xCross/2.0-xPad/2.0
    
    yMax_list = []
    tpr_newList = []
    beta_newList = []

    

    beta = beta_list[-2]
    alpha_Ini   = tpr_list[0].findAlpha(beta,"average")
    alpha_Final = tpr_list[1].findAlpha(beta,"upper")

#    yMax_list.append(sin(alpha_Ini))
    beta_newList.append(beta)
    beta_newList.append(beta)
    height = tan(alpha_Ini+beta)*width * 1.05
    yMax_list.append(height/width)
    
    if iFinalState==0:
        yPlot=0.8
        text = 'C1'
        yShift = -.075
    else:
        yPlot=-.025
        text = 'C3'
        yShift = -0.075
    

    plt.text(beta*deg,alpha_Ini*deg+1,text + 'i',horizontalAlignment='center') 
    plt.text(beta*deg,alpha_Final*deg-1,text + 'f',verticalAlignment='top',horizontalAlignment='center')     
    
    drawAxes['1%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ 0.0,  yPlot,  width,height])['11']
    
#    drawAxes['1%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ 0.0,  -0.2,  width,height])['11']
#    drawAxes['1%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ (beta*deg-x0)/(x1-x0)+xPad,  (alpha_Ini*deg-y0)/(y1-y0)+yPad,  width,height])['11']
#    height = sin(alpha_Ini)*width
    height = tan(alpha_Final+beta)*width * 1.05
    yMax_list.append(height/width)
    drawAxes['2%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ 0.0+width+xPad-xCross, yPlot+yShift,  width,height])['11']
#    drawAxes['2%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ (beta*deg-x0)/(x1-x0)-xPad-width,  (alpha_Final*deg-y0)/(y1-y0)-yPad-height,  width,height])['11']
#    
    plt.sca(graphAxes['13'])
    plt.axis([.0,1.0,.0,1.0])
    
    if text=='C1':
        plt.text(0.0,yPlot+0.05,text+'i')
        plt.text(1.0,yPlot+0.025,text+'f',horizontalAlignment='right')
    if iFinalState==0:
        yPlot=0.35
        yShift = -.175
        text = 'C2'
        beta = beta_list[4]
        beta_newList.append(beta)
        beta_newList.append(beta)
        alpha_Ini   = tpr_list[0].findAlpha(beta,"average")
        alpha_Final = alpha_Ini
        
        print("alphaFinal+beta = %.5f" % (alpha_Final+beta))

        plt.sca(graphAxes['1%i' % (iFinalState+1)])
        plt.text(beta*deg,alpha_Ini*deg+0.75,text + 'i,f') 
        
        height = tan(alpha_Ini+beta)*width * 1.05
        yMax_list.append(height/width)
        drawAxes['3%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ 0.0,  yPlot,  width,height])['11']
#        drawAxes['3%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ (beta*deg-x0)/(x1-x0)+xPad,  (alpha_Ini*deg-y0)/(y1-y0)+yPad,  width,height])['11']
        height = tan(alpha_Final+beta)*width * 1.05
        yMax_list.append(height/width)
        drawAxes['4%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['13'],1,1,box=[ 0.0+width+xPad-xCross, yPlot+yShift,  width,height])['11']
#        drawAxes['4%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ 0.0+width+xPad,  -0.5,  width,height])['11']
#        drawAxes['4%i' % (iFinalState+1)] = Figz_Utils.makeSubAxes(graphAxes['1%i' % (iFinalState+1)],1,1,box=[ (beta*deg-x0)/(x1-x0)-xPad-width,  (alpha_Final*deg-y0)/(y1-y0)-yPad-height,  width,height])['11']
        
        
    plt.sca(graphAxes['13'])
    
    if text=='C3':
        plt.text(0.0,yPlot+0.05,text+'i')
        plt.text(1.0,yPlot-0.00,text+'f',horizontalAlignment='right')
    else:
        plt.text(0.0,yPlot+0.07,text+'i')
        plt.text(1.0,yPlot+0.125,text+'f',horizontalAlignment='right')
    
#    drawAxes['2%i' % (iFinalState+1)] = Figz_Utils.makeAxes(fig,1,1,box=[3,0,width,height])['11']
    



    for i in range(len(drawAxes)):
        plt.sca(drawAxes['%i%i' % (i+1,iFinalState+1)])
        ax = plt.gca()
        ax.patch.set_color('b')
        ax.patch.set_alpha(0.5)
        plt.axis([.0,1.0,.0,yMax_list[i]])
        plt.axis('off')
    #    ax.bac
    

    
    Itpr = [0,1,0,0]
    ItprColor = [0,1,0,1]
    
    for i in range(len(drawAxes)):
        plt.sca(drawAxes['%i%i' % (i+1,iFinalState+1)])
        origin = arr([0.00,0.00 ])
        
        if i==0 or i==2 or (iFinalState==1 and i==1):
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=True,plotStress=False,
                      origin=origin,
                      fx0_list_a = arr([0.5]),
                      fy0_list_a = arr([0.003]),
                      fx0_list_b = arr([]),
                      fy0_list_b = arr([]),
                      faultPos = 0.40,
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
            
            
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=True,plotWedge=False,plotFaultsArrow=False,plotStress=False,
                      origin=origin,
                      fx0_list_a = arr([0.5]),
                      fy0_list_a = arr([0.15, 0.3]),
                      fx0_list_b = arr([0.3, 0.5, .7 ]),
                      fy0_list_b = arr([0.01]),
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
            
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=True,plotWedge=False,plotStress=False,
                      origin=origin,
                      fx0_list_a = arr([]),
                      fy0_list_a = arr([]),
                      fx0_list_b = arr([.9 ]),
                      fy0_list_b = arr([0.01]),
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
        elif i==3:
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=False,plotStress=False,plotFaultsArrow=False,
                      origin=origin,
                      fx0_list_a = arr([0.4]),
                      fy0_list_a = arr([0.003]),
                      fx0_list_b = arr([.4,.8]),
                      fy0_list_b = arr([0.003]),
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
        else:
            truc=0
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=True,plotStress=False,plotFaultsArrow=False,
                      origin=origin,
                      fx0_list_a = arr([0.4]),
                      fy0_list_a = arr([0.003]),
                      fx0_list_b = arr([.4,.8]),
                      fy0_list_b = arr([0.003]),
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
            plotWedge(tpr_list[Itpr[i]],'upper',beta=beta_newList[i],plotFaults=True,plotStress=False,
                      origin=origin,
                      fx0_list_a = arr([.8]),
                      fy0_list_a = arr([0.003]),
                      fx0_list_b = arr([.6]),
                      fy0_list_b = arr([0.003]),
        #              sx0 = sx0[iTpr], sy0 = sy0[iTpr],
                      faultPos = 0.5,
                      f_arrowLength=.04, f_arrowHeadLength = .5,f_arrowSpacing=0.02,
                      colorWedge=Color_w_transparency[ItprColor[i]],
                      colorFaults='k')
#    
    
    
    print("text")
    print(arr(beta_newList)*deg)
    
    
#    # test
#    if iFinalState==0:
#        betas = np.linspace(0.0,35.0,100)/deg
#        alphasRef = np.zeros(betas.shape)
#        alphasFinal = np.zeros(betas.shape)
#        alphasUp = np.zeros(betas.shape)
#        alphasDown = np.zeros(betas.shape)
#        for iB in range(len(betas)):
#            beta = betas[iB]
#            alphasRef[iB] = tpr_list[0].findAlpha(beta,"average")
#            alphasFinal[iB] = tpr_list[1].findAlpha(beta,"average")
#            alphasUp[iB] = tpr_list[1].findAlpha(beta,"upper")
#            alphasDown[iB] = tpr_list[1].findAlpha(beta,"lower")
#            
#        Icross = np.argmin(np.abs(alphasRef-alphasFinal))
#        
#        
#        alphasNorm = (alphasRef-alphasDown)*(alphasUp-alphasDown)
#        psiNorm = (tpr_list[0].psi_bmax - tpr_list[1].psi_bmin)*(tpr_list[1].psi_bmax -tpr_list[1].psi_bmin)
#        Icross2 = np.argmin(np.abs(alphasNorm-psiNorm))
#        
#        Ip1 = np.argmax(alphasUp)
#        print("===")
#        print(tpr_list[1].psi_bmin*deg)
#        print(tpr_list[1].psi_bmax*deg)
#        print("===")
        

#plt.savefig("/Users/abauville/Output/Paper_Decollement/Figz/CritTaper/Fig03",dpi=300)