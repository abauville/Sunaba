#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:13:13 2018

@author: abauville
"""

import sys
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
import CritTaper_dataMaker
import Figz_Utils
import CritTaper_Style
from numpy import array as arr
from PaperDecollement_Utils import getColormap, get_XYandPattern


#(nChi, nBeta, nLambda, LambdaRef_list, 
# chi_list, betas_all, alphas_Ref_all, 
# alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False,nChi=51, nBeta=101, nLambda = 51, enveloppeRes = 6001)
#

alphas_diff_all = alphas_WB_up_all - alphas_Ref_all
taper_angles_all  = betas_all+alphas_Ref_all


fig = Figz_Utils.Figure(8,height=21.0)
myAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.00)

betas = np.linspace(0.0,30.0,100)*np.pi/180.0
taper_angles = np.linspace(-0.0,45.0,20)*np.pi/180.0
nB = len(betas)
nTA = len(taper_angles)

#
# Plot 1 chi_transition vs Taper angle
chi_transition_I_II = np.zeros((nLambda,nTA))
chi_transition_II_III = np.zeros((nLambda,nTA))
alphas_diff = np.zeros(nChi)
#alphas_diff_all2 = np.zeros((nLambda,nChi,nB))
for iL in range(0,nLambda,1):
#    for iB in range(nB):
    for iTA in range(nTA):
#        beta = betas[iB]
        taper_angle = taper_angles[iTA]
        for iW in range(nChi):
#            IB = np.argmin(abs(betas_all[iL,iW,:]-beta))
            IB = np.argmin(abs(taper_angles_all[iL,iW,:]-taper_angle))
            alphas_diff[iW] = alphas_diff_all[iL,iW,IB]#/taper_angles_all[iL,iW,IB]
#            alphas_diff_all2[iL,iW,iB] = alphas_diff_all[iL,iW,IB]
        # end iW
        chi_transition_I_II [iL, iTA] = chi_list[np.argmax(alphas_diff)]
        chi_transition_II_III[iL, iTA] = chi_list[np.argmin(np.abs(alphas_diff))]
    # end iB
#    plt.plot(betas*180.0/np.pi, chi_transition[iL, :],'-o')    
    plt.plot(taper_angles*180.0/np.pi, chi_transition_I_II[iL, :],'-')    
#    plt.plot(taper_angles*180.0/np.pi, chi_transition_II_III[iL, :],'o-')    
# end iL

#        taper_angles[iL,iW]  = betas_all[iL,iW,IB]+alphas_Ref_all[iL,iW,IB]
    
#    

##
## Plot 1
#chi_transition_I_II = np.zeros((nLambda,nTA))
#chi_transition_II_III = np.zeros((nLambda,nTA))
#alphas_diff = np.zeros(nChi)
##alphas_diff_all2 = np.zeros((nLambda,nChi,nB))
#for iL in range(0,nLambda,1):
##    for iB in range(nB):
#    for iTA in range(nTA):
##        beta = betas[iB]
#        taper_angle = taper_angles[iTA]
#        for iW in range(nChi):
##            IB = np.argmin(abs(betas_all[iL,iW,:]-beta))
#            IB = np.argmin(abs(taper_angles_all[iL,iW,:]-taper_angle))
#            alphas_diff[iW] = alphas_diff_all[iL,iW,IB]#/taper_angles_all[iL,iW,IB]
##            alphas_diff_all2[iL,iW,iB] = alphas_diff_all[iL,iW,IB]
#        # end iW
#        chi_transition_I_II [iL, iTA] = chi_list[np.argmax(alphas_diff)]
#        chi_transition_II_III[iL, iTA] = chi_list[np.argmin(np.abs(alphas_diff))]
#    # end iB
##    plt.plot(betas*180.0/np.pi, chi_transition[iL, :],'-o')    
##    plt.plot(taper_angles*180.0/np.pi, chi_transition_I_II[iL, :],'x')    
#    plt.plot(taper_angles*180.0/np.pi, chi_transition_II_III[iL, :],'o-')    
## end iL
#
##        taper_angles[iL,iW]  = betas_all[iL,iW,IB]+alphas_Ref_all[iL,iW,IB]
#    
##    






#chi_transition_I_II = np.zeros((nLambda,nB))
#chi_transition_II_III = np.zeros((nLambda,nB))
#
#taperAngle_transition_I_II = np.zeros((nLambda,nB))
#taperAngle_transition_II_III = np.zeros((nLambda,nB))
#
#alphas_diff = np.zeros(nChi)
#alphas = np.zeros(nChi)
#
##alphas_diff_all2 = np.zeros((nLambda,nChi,nB))
#for iL in range(0,nLambda,10):
#    for iB in range(nB):
#        beta = betas[iB]
#        
#        for iW in range(nChi):
#            IB = np.argmin(abs(betas_all[iL,iW,:]-beta))
#            alphas_diff[iW] = alphas_diff_all[iL,iW,IB]#/taper_angles_all[iL,iW,IB]
#            
#            alphas[iW] = alphas_Ref_all[iL,iW,IB]
##            alphas_diff_all2[iL,iW,iB] = alphas_diff_all[iL,iW,IB]
#        # end iW
#        IW = np.argmax(alphas_diff)
#        IB = np.argmin(abs(betas_all[iL,IW,:]-beta))
#        taperAngle_transition_I_II [iL, iB] = taper_angles_all[iL,IW,IB]
#        IW = np.argmin(np.abs(alphas_diff))
#        IB = np.argmin(abs(betas_all[iL,IW,:]-beta))
#        taperAngle_transition_II_III[iL, iB] = taper_angles_all[iL,IW,IB]
#    # end iB
##    plt.plot(betas*180.0/np.pi, chi_transition[iL, :],'-o')    
##    plt.plot(betas*180.0/np.pi, taperAngle_transition_I_II[iL, :]*180.0/np.pi,'-x')    
##    plt.plot(betas*180.0/np.pi, taperAngle_transition_II_III[iL, :]*180.0/np.pi,'-o')  
#
##    plt.plot(betas*180.0/np.pi, taperAngle_transition_II_III[iL, :]*180.0/np.pi,'-o')  
#    plt.plot(betas*180.0/np.pi, taperAngle_transition_II_III[iL, :]*180.0/np.pi,'-o')  
#    plt.axis([0.0,30.0,0.0,30.0])
## end iL
#
##        taper_angles[iL,iW]  = betas_all[iL,iW,IB]+alphas_Ref_all[iL,iW,IB]
#    
    