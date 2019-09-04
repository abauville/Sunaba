#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:06:54 2018

@author: abauville
"""

import matplotlib.pyplot as plt
import numpy as np

#class customFig():
#    width = 0.0
#    height = 0.0

class Figure():
    
    def __init__( self, figNum=1,orientation="portrait", 
                  width="default", height="default", 
                  leftMargin = 1.25, rightMargin = 1.25, 
                  bottomMargin = 1.25, topMargin = 1.25,
                  mode="production", dpi=220):
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
        
        
        if width!="default":
            if not np.isreal(width):
                raise ValueError("figWidth should be a number")
            figW = width
        if height!="default":
            if not np.isreal(height):
                raise ValueError("figHeight should be a number")
            figH = height
        
        if mode=="crop":
            figW = figW-leftMargin-rightMargin
            figH = figH-topMargin-bottomMargin
            leftMargin = 0.0
            rightMargin= 0.0
            topMargin = 0.0
            bottomMargin = 0.0
            
        fig = plt.figure(figNum)
        fig.set_size_inches(figW*cm2inch,figH*cm2inch, forward=True)
        fig.set_dpi(dpi)
        plt.clf()
    
    
        
    
        self.width = figW
        self.height = figH
        self.usableWidth = figW-leftMargin-rightMargin
        self.usableHeight = figH-topMargin-bottomMargin
        self.leftMargin = leftMargin
        self.rightMargin = rightMargin
        self.bottomMargin = bottomMargin
        self.topMargin = topMargin
        self.backgroundAxes = plt.axes([0.0,0.0,1.0,1.0],label='background')
        self.figure = fig
        
        if mode=="draft":
    #    xMargin = 1.25 * cm
    #    yMargin = 1.5 * cm
            plt.fill([0.0, leftMargin, leftMargin, 0.0], [0.0, 0.0, figH, figH],color=[.9,.9,.9])
            plt.fill([figW, figW-rightMargin, figW-rightMargin, figW], [0.0, 0.0, figH, figH],color=[.9,.9,.9])
            
            plt.fill([0.0, 0.0, figW, figW], [0.0, bottomMargin, bottomMargin, 0.0],color=[.9,.9,.9])
            plt.fill([0.0, 0.0, figW, figW], [figH, figH-topMargin, figH-topMargin, figH],color=[.9,.9,.9])
        elif mode=="production":
            daijoubu=1
        elif mode == "crop":
            daijoubu=1
        else: 
            raise ValueError("Unknwon mode, should be 'draft', 'production' or 'crop'.")
        plt.axis([0.0,figW,0.0,figH])
        plt.axis("off")






def makeAxes(   figure, nVertical=1,nHorizontal=1, 
                xPad = 0.5, yPad = 0.5,
                leftMarginPad = 1.0, rightMarginPad = 0.25,
                bottomMarginPad = 0.0, topMarginPad = 0.25,
                box = 'automatic',
                setAspectRatioBasedOn="x", aspectRatio="default"):  
    
    # box contains [x0,y0,width,height], where x0,y0 are the coordinates of the bottom left corner and coordinates start from the bottom left corner of the usable fig
    
    if not isinstance(figure, Figure):
        raise TypeError('the figure obj passed as first argument should be an instance of the class Figz_Utils.Figure')
    
    
    plt.figure(figure.figure.number)
    myAxes = {}
    
    figW = figure.width
    figH = figure.height
    
    if box!='automatic':
        if type(box)==int:
            raise AttributeError('box should be [x0,y0,width,height]')
        if type(box)==list:
            if len(box)!=4:
                raise ValueError('box should be [x0,y0,width,height]')
        if type(box)==np.ndarray:
            if len(box)!=4 or box.size!=4:
                raise ValueError('box should be [x0,y0,width,height]')

        x0 = box[0]
        y0 = box[1]
        totalWidth = box[2]            
        totalHeight = box[3]
        
        leftMarginPad   = x0
        rightMarginPad  = figure.width-figure.leftMargin-figure.rightMargin-leftMarginPad-totalWidth
#        bottomMarginPad = figure.height-figure.bottomMargin-figure.topMargin-y0
#        topMarginPad = figure.height-figure.bottomMargin-figure.topMargin-totalHeight-y0
        bottomMarginPad = figure.height-figure.bottomMargin-figure.topMargin-y0
        topMarginPad = y0-totalHeight
        
            
        
    #plotsH = 5.0 * cm

#    xPad = 0.5 * cm
#    yPad = 0.5 * cm
    leftMargin = figure.leftMargin + leftMarginPad
    rightMargin = figure.rightMargin + rightMarginPad
    bottomMargin = figure.bottomMargin + bottomMarginPad
    topMargin = figure.topMargin + topMarginPad
    plotsW = (figW-leftMargin-rightMargin-(nHorizontal-1)*xPad)/nHorizontal
    plotsH = (figH-bottomMargin-topMargin-(nVertical-1)*yPad)/nVertical
#    print(plotsW)
#    print(plotsH)
    if aspectRatio != "default":
        if not np.isreal(aspectRatio):
            raise ValueError("aspectRatio should be a number")

        if setAspectRatioBasedOn=="x":
            plotsH = aspectRatio*plotsW
        elif setAspectRatioBasedOn=="y":
            plotsW = aspectRatio*plotsH
        else:
            raise ValueError("setAspectRatioBasedOn must be 'x' or 'y'.")
            
    scalePosition = np.array([1.0/figW, 1.0/figH, 1.0/figW, 1.0/figH])
    
#    print(plotsW)
#    print(plotsH)
    
    for i in range (nVertical):
        for j in range (nHorizontal):
            myAxes['%i%i' % (i+1,j+1)] = plt.axes(np.array([leftMargin+j*(plotsW+xPad), figH-topMargin-(i+1)*(plotsH)-i*yPad,plotsW, plotsH])*scalePosition,label='%i' % int(np.random.rand(1)[0]*1000000000))

    myAxes['info'] = {}
    myAxes['info']['plotsWidth']    = plotsW
    myAxes['info']['plotsHeight']   = plotsH
    myAxes['info']['leftMargin']    = leftMargin
    myAxes['info']['rightMargin']   = rightMargin
    myAxes['info']['topMargin']     = topMargin
    myAxes['info']['bottomMargin']  = bottomMargin
    myAxes['info']['leftMarginPad']    = leftMarginPad
    myAxes['info']['rightMarginPad']   = rightMarginPad
    myAxes['info']['topMarginPad']     = topMarginPad
    myAxes['info']['bottomMarginPad']  = bottomMarginPad
    myAxes['info']['xPad']          = xPad
    myAxes['info']['yPad']          = yPad
    myAxes['info']['parentFigure']        = figure
    myAxes

    return myAxes

def makeSubAxes(parentAxes, nVertical=1,nHorizontal=1, 
                xPad = 0.2, yPad = 0.2,
                box = 'automatic'): 
    
    # box should be [x0,y0,x1,y1] x0,y0 are coordinate of the bottom left of this axes from the bottom left corner of the parent axes
    # coordinates are normalized to the size of the parent axes
    
    x0 = box[0]
    y0 = box[1]
    width=box[2]
    height=box[3]
    
    plotsW = (width-(nHorizontal-1)*xPad)/nHorizontal
    plotsH = (height-(nVertical-1)*yPad)/nVertical
    
    
#    scalePosition = np.array([1.0/figW, 1.0/figH, 1.0/figW, 1.0/figH])
    
    # Parent axes position
    px0 = parentAxes.get_position().get_points()[0,0]
    py0 = parentAxes.get_position().get_points()[0,1]
    px1 = parentAxes.get_position().get_points()[1,0]
    py1 = parentAxes.get_position().get_points()[1,1]   
    
    pw = px1-px0
    ph = py1-py0
    
    x0 = px0+box[0]*pw
    y0 = py0+box[1]*ph
    
    w = box[2]*pw
    h = box[3]*ph
    
    
    
    myAxes = {}
    for i in range (nVertical):
        for j in range (nHorizontal):
            myAxes['%i%i' % (i+1,j+1)] = plt.axes(np.array([x0,y0,w,h]))
    
    return myAxes


    
    
    
    
    
    
    
    