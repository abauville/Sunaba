#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 16:09:39 2018

@author: abauville
"""
import sys
import matplotlib.pyplot as plt
dpi = 800
OutputFolder = "/Users/abauville/Output/Paper_Decollement/Figz/Numerics/"

#
#plt.close('all')
#print("Fig01")
#import Num_Fig01
#plt.savefig(OutputFolder + "03_Numerics_Fig01",dpi=dpi)
#
#plt.close('all')
#print("Fig02")
#import Num_Fig02
#plt.savefig(OutputFolder + "03_Numerics_Fig02",dpi=dpi)

#plt.close('all')
#print("Fig03")
#import Num_Fig03
#plt.savefig(OutputFolder + "03_Numerics_Fig03",dpi=dpi)
#
#plt.close('all')
#print("Fig04")
#import Num_Fig04
#plt.savefig(OutputFolder + "03_Numerics_Fig04",dpi=dpi)

sys.path.insert(0, '../ResTest')
plt.close('all')
print("FigResTest")
import FigResTest
plt.savefig(OutputFolder + "04_ResTest_Fig01",dpi=dpi)