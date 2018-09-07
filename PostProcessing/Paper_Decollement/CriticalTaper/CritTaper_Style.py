#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:27:07 2018

@author: abauville
"""
#import numpy as np

class Style():
    Setting = "Paper"
#    Setting = "Presentation"
    if Setting == "Paper":
        fontdict = {'family': 'Montserrat',
                    'weight': 'bold',
                    'size': int(11)
                    }
    
    elif Setting == "Presentation":
        fontdict = {'family': 'Montserrat',
                    'weight': 'bold',
                    'size': int(18)
                    }
    
    colormap = "seismic"