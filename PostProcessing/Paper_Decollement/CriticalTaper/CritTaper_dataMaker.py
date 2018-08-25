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


def getCritTaperFigData(nChi=51, nBeta=51, nLambda = 51, Compute=False):
    LambdaRef_list =np.linspace(0,1.0,nLambda)
    LambdaRef_list[ 0] += 1e-10
    LambdaRef_list[-1] -= 1e-10

    
    
    chi_list = np.linspace(0.01,0.99,nChi)
    
    beta = 0.0
    
    enveloppeRes = 2001
    
    
    
    Taper_WB = []
    Taper_WF = [] # Fully chiened
    Taper_Ref = []
        
    
    alphas_Ref = np.zeros(nBeta)
    alphas_WB_up = np.zeros(nBeta)
    alphas_WB_low = np.zeros(nBeta)
    
    betas_all = np.zeros((nLambda,nChi,nBeta))
    chis_all = np.zeros((nLambda,nChi,nBeta))
    alphas_Ref_all = np.zeros((nLambda,nChi,nBeta))
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
            thisTaper.computeAlphaVsBeta(n=2010)
            
            betaMinRef = np.min(thisTaper.beta_all)
            betaMaxRef = np.max(thisTaper.beta_all)
            
            Taper_Ref.append(thisTaper)
            ## ========================================
            
            for iChi in range(nChi):
    
                
                chiFac = chi_list[iChi]
                Lambdachi = (1.0-chiFac) * LambdaRef   + chiFac
            
                thisTaper = Taper(phi=phiRef, phi_b=phiRef,
                                      Lambda=LambdaRef, Lambda_b=Lambdachi,
                                      rho_w=rho_w, rho=rho)
                
                
                
                
                thisTaper.computeAlphaVsBeta(n=enveloppeRes)
                
                betaMinChiB = np.min(thisTaper.beta_all)
                betaMaxWB = np.max(thisTaper.beta_all)
                                   
                Taper_WB.append(thisTaper)
                
                betaMin = np.max([betaMinRef,betaMinChiB])
                betaMax = np.min([betaMaxRef,betaMaxWB])
                
    
                betas_all[iTaper,iChi,:] = np.linspace(betaMin,betaMax, nBeta)
                
                ## Use the exact value of beta = 0.0 in this array (to produce smoother graph afterward)
                I = np.argmin(abs(betas_all[iTaper,iChi,:] - 0.0))
                betas_all[iTaper,iChi,I] = 0.0 
                
                ## Fully chiened
                thisTaper = Taper(phi=phiRef, phi_b=phiRef-1e-6,
                                  Lambda=Lambdachi, Lambda_b=Lambdachi,
                                  rho_w=rho_w, rho=rho)            
                thisTaper.computeAlphaVsBeta(n=enveloppeRes)
                Taper_WF.append(thisTaper)
                
                
                for iB in range(nBeta):
                    beta = betas_all[iTaper,iChi,iB]
                    alphas_Ref[iB]  = Taper_Ref[iTaper].findAlpha(beta,"average")
                    
                    alphas_WB_up[iB] = Taper_WB[iTaper*nChi+iChi].findAlpha(beta,"upper",tol=1e-3)
                    alphas_WB_low[iB] = Taper_WB[iTaper*nChi+iChi].findAlpha(beta,"lower",tol=1e-3)
                    
                    Lambdas_Ref_all[iTaper,iChi,iB] = LambdaRef
                # end iB
                    
    #                alphas_Diff[iB] = alphas_Ref[iB] - alphas_WB_up[iB]
                
                alphas_Ref_all[iTaper,iChi,:] = alphas_Ref
                alphas_WB_up_all[iTaper,iChi,:] = alphas_WB_up
                alphas_WB_low_all[iTaper,iChi,:] = alphas_WB_low            
    #            alphas_Diff_all[iTaper,iW,:] = alphas_Diff
                
                Counter+=1
                # end iChi
        # end for iTaper
        np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz",
                 nChi = nChi,
                 nBeta = nBeta,
                 nLambda = nLambda,
                 LambdaRef_list = LambdaRef_list,
                 chi_list = chi_list,
                 betas_all = betas_all,
                 alphas_Ref_all = alphas_Ref_all,
                 alphas_WB_up_all = alphas_WB_up_all,
                 alphas_WB_low_all = alphas_WB_low_all,
                 Lambdas_Ref_all = Lambdas_Ref_all,
                 chis_all = chis_all,
                 Taper_Ref = Taper_Ref,
                 Taper_WB = Taper_WB,
                 Taper_WF = Taper_WF
                 
                 )
        
    else: #if Compute   
        loadedData = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/CritTaper_Type1vs2.npz")
        nChi = loadedData["nChi"][()]
        nBeta = loadedData["nBeta"][()]
        nLambda = loadedData["nLambda"][()]
        LambdaRef_list = loadedData["LambdaRef_list"][()]
        chi_list = loadedData["chi_list"][()]
        betas_all = loadedData["betas_all"][()]
        alphas_Ref_all = loadedData["alphas_Ref_all"][()]
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
            alphas_WB_up_all,
            alphas_WB_low_all,
            Lambdas_Ref_all,
            chis_all,
            Taper_Ref,
            Taper_WB,
            Taper_WF]
    
