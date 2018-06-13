#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 17:59:27 2018

@author: abauville
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import exp, sqrt, pi
from numpy.fft import fft2, ifft2
def createRandomSurface(N,rL=1.0,h=1.0,clx=0.1,cly=0.1,angle=0.0):
    # 2D Gaussian random rough surface with Gaussian autocovariance function
    x = np.linspace(-1.0,1.0,N)
    y = np.linspace(-1.0,1.0,N)
    
    X, Y =np.meshgrid(x,y)
    X = X.T
    Y = Y.T

    Z = np.random.randn(N,N)
    #Z(:) = h*randraw( SurfFun, distribParams, N*N );
    
    # Case for 0 angle
    #F = exp(-(X.^2/(clx^2/2)+Y.^2/(cly^2/2))); # Gaussian filter
    #f = 2/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z).*fft2(F)); # correlated surface generation including convolution (faltning) and inverse Fourier transform
    
    
    
#    if angle == 0.0:
    F = exp(-(X**2/(0.5*clx**2)+Y**2/(0.5*cly**2))); # Gaussian filter
    f = 2.0/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z)*fft2(F)); # correlated surface generation including convolution (faltning) and inverse Fourier transform
#   else:
#       #tilted Gaussian
#       a = ((cosd(angle)^2) / (2*clx^2)) + ((sind(angle)^2) / (2*cly^2));
#       b = -((sind(2*angle)) / (4*clx^2)) + ((sind(2*angle)) / (4*cly^2));
#       c = ((sind(angle)^2) / (2*clx^2)) + ((cosd(angle)^2) / (2*cly^2));
#        
#       F = exp(-(a*(X).^2 + 2*b*(X).*(Y) + c*(Y).^2)); % tilted Gaussian filter
#       f = 2/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z).*fft2(F)); % correlated surface generation including convolution (faltning) and inverse Fourier transform

    
    
    f = f.astype("float64")
    
    #plt.pcolor(X,Y,f)
    
    return f



