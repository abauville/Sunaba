#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 12:59:21 2018

@author: abauville
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 11:21:03 2018

@author: abauville
"""

## Illsutration of stress and fault orientations in accretionary prism
## Based on the critical taper theory
## Nomenclature follows Buiter, 2011, "A review of brittle compressional wedge models"


import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan
import numpy.matlib as matlib
from CritTaper_utils import Taper


def getCritTaperFigData(nChi=51, nBeta=51, nLambda = 51, Compute=False, beta_list = np.array([]),
                        enveloppeRes = 2001, alphaMin="default",alphaMax="default"):
    LambdaRef_list =np.linspace(0,1.0,nLambda)
    LambdaRef_list[ 0] += 1e-10
    LambdaRef_list[-1] -= 1e-10

    if isinstance(beta_list, np.float):
        fixedBetas = True
        nBeta = 1
        if np.abs(beta_list)>np.pi/2.0:
               raise Warning("beta_list should be given in radians")
    else:
        if len(beta_list) == 0:
            fixedBetas = False
        else:
            fixedBetas = True
            nBeta = len(beta_list)
            for beta in beta_list:
               if np.abs(beta)>np.pi/2.0:
                   raise Warning("beta_list should be given in radians")

        

    
    chi_list = np.linspace(0.0,1.0,nChi)
    chi_list[ 0] += 1e-10
    chi_list[-1] -= 1e-10
    
    beta = 0.0
    
    
    
    
    Taper_WB = []
    Taper_WF = [] # Fully weakened
    Taper_Ref = []
        
    
    alphas_Ref = np.zeros(nBeta)
    alphas_WF = np.zeros(nBeta)
    alphas_WB_up = np.zeros(nBeta)
    alphas_WB_low = np.zeros(nBeta)
    
    betas_all = np.zeros((nLambda,nChi,nBeta))
    chis_all = np.zeros((nLambda,nChi,nBeta))
    alphas_Ref_all = np.zeros((nLambda,nChi,nBeta))
    alphas_WF_all   = np.zeros((nLambda,nChi,nBeta))
    alphas_WB_up_all = np.zeros((nLambda,nChi,nBeta))
    alphas_WB_low_all = np.zeros((nLambda,nChi,nBeta))
    
    Lambdas_Ref_all = np.zeros((nLambda,nChi,nBeta))
    
    
    chi_small = chi_list.copy()
    chi_small = matlib.repmat(chi_small,nBeta,1)
    chi_small = chi_small.T
    
    
    
    if Compute:
        Counter = 0
        maxCounter = nLambda*nChi
        for iTaper in range(nLambda):
            chis_all[iTaper,:,:] = chi_small
            
            if Counter%1==0:
                print("Counter = %i/%i" % (Counter,maxCounter))
            
            rho_w = 1000.0
            rho = 2500.0
            phiRef   = 30.0*pi/180.0
            LambdaRef=LambdaRef_list[iTaper]
            
            
            
            ## ============= RefTaper =================    
            thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                             Lambda=LambdaRef, Lambda_b=LambdaRef,
                             rho_w=rho_w, rho=rho)
            thisTaper.computeAlphaVsBeta(n=enveloppeRes,alphaMin=alphaMin,alphaMax=alphaMax)
            
            betaMinRef = np.min(thisTaper.beta_all)
            betaMaxRef = np.max(thisTaper.beta_all)
            
            Taper_Ref.append(thisTaper)
            ## ========================================
            
            for iChi in range(nChi):
    
                
                chiFac = chi_list[iChi]
                LambdaWeak = (1.0-chiFac) * LambdaRef   + chiFac
            
                ## ============= Weak base =================    
                thisTaper = Taper(phi=phiRef, phi_b=phiRef,
                                      Lambda=LambdaRef, Lambda_b=LambdaWeak,
                                      rho_w=rho_w, rho=rho)    
                thisTaper.computeAlphaVsBeta(n=enveloppeRes,alphaMin=alphaMin,alphaMax=alphaMax)
                
                betaMinWB = np.min(thisTaper.beta_all)
                betaMaxWB = np.max(thisTaper.beta_all)              
                Taper_WB.append(thisTaper)
                ## =========================================   
                
                ## ============= Fully weakened =================    
                thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                                      Lambda=LambdaWeak, Lambda_b=LambdaWeak,
                                      rho_w=rho_w, rho=rho)    
                thisTaper.computeAlphaVsBeta(n=enveloppeRes,alphaMin=alphaMin,alphaMax=alphaMax)
#                
#                betaMinWB = np.min(thisTaper.beta_all)
#                betaMaxWB = np.max(thisTaper.beta_all)              
                Taper_WF.append(thisTaper)
                ## =========================================  
                
                if fixedBetas:
                    betas_all[iTaper,iChi,:] = beta_list
                else:
                    betaMin = np.max([betaMinRef,betaMinWB])
                    betaMax = np.min([betaMaxRef,betaMaxWB])
                    betas_all[iTaper,iChi,:] = np.linspace(betaMin,betaMax, nBeta)
                    ## Use the exact value of beta = 0.0 in this array (to produce smoother graph afterward)
                    I = np.argmin(abs(betas_all[iTaper,iChi,:] - 0.0))
                    betas_all[iTaper,iChi,I] = 0.0 
                

                
                
                for iB in range(nBeta):
                    beta = betas_all[iTaper,iChi,iB]
                    alphas_Ref[iB]  = Taper_Ref[iTaper].findAlpha(beta,"average")
                    alphas_WF [iB]  = Taper_WF [iTaper].findAlpha(beta,"average")
                    
                    alphas_WB_up[iB] = Taper_WB[iTaper*nChi+iChi].findAlpha(beta,"upper",tol=1e-3)
                    alphas_WB_low[iB] = Taper_WB[iTaper*nChi+iChi].findAlpha(beta,"lower",tol=1e-3)
                    
                    Lambdas_Ref_all[iTaper,iChi,iB] = LambdaRef
                # end iB
                    
    #                alphas_Diff[iB] = alphas_Ref[iB] - alphas_WB_up[iB]
                
                alphas_Ref_all[iTaper,iChi,:] = alphas_Ref
                alphas_WF_all[iTaper,iChi,:] = alphas_WF
                alphas_WB_up_all[iTaper,iChi,:] = alphas_WB_up
                alphas_WB_low_all[iTaper,iChi,:] = alphas_WB_low            
    #            alphas_Diff_all[iTaper,iW,:] = alphas_Diff
                
                Counter+=1
                # end iChi
        # end for iTaper
        if fixedBetas:
            fileName  = "/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2_fixedBetas.npz"
        else:
            fileName  = "/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz"
        
        np.savez(fileName,
                 nChi = nChi,
                 nBeta = nBeta,
                 nLambda = nLambda,
                 LambdaRef_list = LambdaRef_list,
                 chi_list = chi_list,
                 betas_all = betas_all,
                 alphas_Ref_all = alphas_Ref_all,
                 alphas_WF_all = alphas_WF_all,
                 alphas_WB_up_all = alphas_WB_up_all,
                 alphas_WB_low_all = alphas_WB_low_all,
                 Lambdas_Ref_all = Lambdas_Ref_all,
                 chis_all = chis_all,
                 Taper_Ref = Taper_Ref,
                 Taper_WB = Taper_WB,
                 Taper_WF = Taper_WF
                 
                 )
        
    else: #if Compute   
        
        if fixedBetas:
            loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2_fixedBetas.npz")
        else:
            loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz")
        nChi = loadedData["nChi"][()]
        nBeta = loadedData["nBeta"][()]
        nLambda = loadedData["nLambda"][()]
        LambdaRef_list = loadedData["LambdaRef_list"][()]
        chi_list = loadedData["chi_list"][()]
        betas_all = loadedData["betas_all"][()]
        alphas_Ref_all = loadedData["alphas_Ref_all"][()]
        alphas_WF_all = loadedData["alphas_WF_all"][()]
        alphas_WB_up_all = loadedData["alphas_WB_up_all"][()]
        alphas_WB_low_all = loadedData["alphas_WB_low_all"][()]
        Lambdas_Ref_all = loadedData["Lambdas_Ref_all"][()]
        chis_all = loadedData["chis_all"][()]
        Taper_Ref = loadedData["Taper_Ref"][()]
        Taper_WB = loadedData["Taper_WB"][()]
        Taper_WF = loadedData["Taper_WF"][()]
        
        
        
    return [nChi,
            nBeta,
            nLambda,
            LambdaRef_list,
            chi_list,
            betas_all,
            alphas_Ref_all,
            alphas_WF_all,
            alphas_WB_up_all,
            alphas_WB_low_all,
            Lambdas_Ref_all,
            chis_all,
            Taper_Ref,
            Taper_WB,
            Taper_WF]
    
