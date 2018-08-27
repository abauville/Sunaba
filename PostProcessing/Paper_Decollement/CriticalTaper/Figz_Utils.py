#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:06:54 2018

@author: abauville
"""

import matplotlib.pyplot as plt
import numpy as np

def makeFigure( nVertical=1,nHorizontal=1, 
                figNum=1,orientation="portrait", 
                figWidth="default", figHeight="default", 
                xMargin = 1.25, yMargin = 1.5,
                xPad = 0.5, yPad = 0.5,
                leftMarginPad = 1.0, rightMarginPad = 0.25,
                bottomMarginPad = 1.0, topMarginPad = 0.0,
                setAspectRatioBasedOn="x", aspectRatio="default",
                mode="production",
                dpi=220):
    # setAspectRatioBasedOn can take "x" or "y"
    # mode: "draft" or "production"
    cm = 1.0
    cm2inch = 0.393701
    if orientation == "portrait":
        figW = 21.0 * cm
        figH = 29.7 * cm
    elif orientation == "landscape":
        figW = 29.7 * cm
        figH = 21.0 * cm 
    else:
        raise ValueError("Unknown orientation %s" % orientation)
    
    
    if figWidth!="default":
        if not np.isreal(figWidth):
            raise ValueError("figWidth should be a number")
        figW = figWidth
    if figHeight!="default":
        if not np.isreal(figHeight):
            raise ValueError("figHeight should be a number")
        figH = figHeight
        
    fig = plt.figure(figNum)
    fig.set_size_inches(figW*cm2inch,figH*cm2inch, forward=True)
    fig.set_dpi(dpi)
    plt.clf()
        
    myAxes = {}
    myAxes['00'] = plt.axes([0.0,0.0,1.0,1.0])
    if mode=="draft":
#    xMargin = 1.25 * cm
#    yMargin = 1.5 * cm
        plt.fill([0.0, xMargin, xMargin, 0.0], [0.0, 0.0, figH, figH],color=[.9,.9,.9])
        plt.fill([figW, figW-xMargin, figW-xMargin, figW], [0.0, 0.0, figH, figH],color=[.9,.9,.9])
        
        plt.fill([0.0, 0.0, figW, figW], [0.0, yMargin, yMargin, 0.0],color=[.9,.9,.9])
        plt.fill([0.0, 0.0, figW, figW], [figH, figH-yMargin, figH-yMargin, figH],color=[.9,.9,.9])
    elif mode=="production":
        daijoubu=1
    else: 
        raise ValueError("Unknwon mode, should be 'draft' or production")
    plt.axis([0.0,figW,0.0,figH])
    plt.axis("off")
    
    #plotsH = 5.0 * cm

#    xPad = 0.5 * cm
#    yPad = 0.5 * cm
    leftMargin = xMargin + leftMarginPad
    rightMargin = xMargin + rightMarginPad
    bottomMargin = yMargin + bottomMarginPad
    topMargin = yMargin + topMarginPad
    plotsW = (figW-leftMargin-rightMargin-(nHorizontal-1)*xPad)/nHorizontal
    plotsH = (figH-bottomMargin-topMargin-(nVertical-1)*yPad)/nVertical
    if aspectRatio != "default":
        if not np.isreal(aspectRatio):
            raise ValueError("aspectRatio should be a number")
            
        if setAspectRatioBasedOn=="x":
            plotsH = plotsW
        elif setAspectRatioBasedOn=="y":
            plotsW = plotsH
        else:
            raise ValueError("setAspectRatioBasedOn must be 'x' or 'y'.")
            
    scalePosition = np.array([1.0/figW, 1.0/figH, 1.0/figW, 1.0/figH])
    
    for i in range (nVertical):
        for j in range (nHorizontal):
            myAxes['%i%i' % (i+1,j+1)] = plt.axes(np.array([leftMargin+j*(plotsW+xPad), bottomMargin+i*(plotsH*yPad),plotsW, plotsH])*scalePosition)

    return [fig, myAxes]