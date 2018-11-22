#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

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

(nChi, nBeta, nLambda, LambdaRef_list, 
 chi_list, betas_all, alphas_Ref_all, 
 alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
 Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=np.linspace(0.0,30.0,13.0)*np.pi/180.0, nChi=61, nLambda=61,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)



beta = 0.0*np.pi/180.0

#(nChi, nBeta, nLambda, LambdaRef_list, 
# chi_list, betas_all, alphas_Ref_all, 
# alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False)

alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
## Lambda vs chi @ beta=0
fig    = Figz_Utils.Figure(5,height=13.0,mode='production')
#fig    = Figz_Utils.Figure(7,height=13.0,mode='draft')
#AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
#AxesDum['12'].axis('off')
#Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=10.5,topMarginPad = 1.0)
xPad = 2.0
Axes   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=1.5,topMarginPad = 1.0,xPad = xPad)
AxesW = Axes['info']['plotsWidth']
AxesH = Axes['info']['plotsHeight']
AxesxPad = Axes['info']['xPad']
AxeslPad = Axes['info']['leftMarginPad']
AxesrPad = Axes['info']['rightMarginPad']
AxestPad = Axes['info']['topMarginPad']
AxesbPad = Axes['info']['bottomMarginPad']
cBaryPad    = 0.0
cBartPad = .75
cBarbPad = 0.0
cBarlPad = 0.4
cBarrPad = 0.0
cBarW = 0.4
#cBarAxes   = Figz_Utils.makeAxes(fig,1,aspectRatio=0.15,leftMarginPad=AxeslPad+cBarlPad,rightMarginPad=AxesrPad+1*AxesxPad+1*AxesW+cBarrPad,topMarginPad=AxesH+cBaryPad)
cBarAxes   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+AxesW+cBarlPad,
                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+AxesW+cBarlPad)-cBarW),
                                       topMarginPad=AxestPad+cBartPad,
                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)
cBarAxes2   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+2.0*AxesW+cBarlPad+xPad,
                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+2.0*AxesW+cBarlPad+xPad)-cBarW),
                                       topMarginPad=AxestPad+cBartPad,
                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)
#Axes['12'].axis('off')
plt.sca(Axes['11'])



deg = 180.00/np.pi


Lambdas = Lambdas_Ref_all[:,:,0]
chis = chis_all[:,:,0]

alphas_diff = np.zeros((nLambda,nChi))
alphas_Ref = np.zeros((nLambda,nChi))
alphas_width = np.zeros((nLambda,nChi))
#chis = np.zeros((nLambda,nChi))
taper_angles = np.zeros((nLambda,nChi))
alphas_WB_up = np.zeros((nLambda,nChi))
alphas_WB_low = np.zeros((nLambda,nChi))

for iL in range(nLambda):
    for iW in range(nChi):
        iB = np.argmin(abs(betas_all[iL,iW,:]-beta))
        alphas_diff[iL,iW]   = alphas_diff_all[iL,iW,iB]
        alphas_width[iL,iW]  = alphas_WB_up_all[iL,iW,iB] - alphas_WB_low_all[iL,iW,iB]
        alphas_WB_up[iL,iW]  = alphas_WB_up_all[iL,iW,iB]
        alphas_WB_low[iL,iW] = alphas_WB_low_all[iL,iW,iB]
        alphas_Ref[iL,iW]    = alphas_Ref_all[iL,iW,iB]
        taper_angles[iL,iW]  = betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB]
        
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,1000)
plt.sca(Axes['11'])
CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(-1.0001,1.0001,20),vmin=-1.00,vmax=1.00)
#plt.text(75,90,"Extensional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)
#plt.text(35,20,"Compressional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)

plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
plt.ylabel("$\\mathbf{\\chi}$ [%]",weight='bold',verticalAlignment='center')

ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([40.0,100.0,100.0,.0])






#   Define interpolated versions of arrays
# ============================================
denseFac = 3
chi_list_dense = np.linspace(chi_list[0],chi_list[-1],denseFac*nChi)
LambdaRef_list_dense = np.linspace(LambdaRef_list[0],LambdaRef_list[-1],denseFac*nLambda)

from scipy import interpolate
f = interpolate.RectBivariateSpline(LambdaRef_list,chi_list,alphas_diff)
alphas_diff_dense = f(LambdaRef_list_dense,chi_list_dense)

f = interpolate.RectBivariateSpline(LambdaRef_list,chi_list,taper_angles)
taper_angles_dense = f(LambdaRef_list_dense,chi_list_dense)

f = interpolate.RectBivariateSpline(LambdaRef_list,chi_list,alphas_WB_up)
alphas_WB_up_dense = f(LambdaRef_list_dense,chi_list_dense)

#f = interpolate.interp1d(LambdaRef_list, chis_vDiv_0)
#chis_vDiv_0_dense = f(LambdaRef_list_dense)




# Add indication of max delta alpha
plt.sca(Axes['11'])
#chis_alpha_diff_min = chi_list_dense[np.argmin(alphas_diff_dense/taper_angles_dense,axis=1)]
#chis_alpha_diff_max = chi_list_dense[np.argmax(alphas_diff_dense/taper_angles_dense,axis=1)]
#
#
#chis_alpha_WB_up_max = chi_list_dense[np.argmax(alphas_WB_up_dense,axis=1)]
#
#
#chis_vDiv_0 = chis_alpha_WB_up_max#chi_list_centered [np.argmin(vDiv,axis=1)]
#chis_vDiv_0[-1] = 0.0#chis_vDiv_0[-2]
#
#


chis_alpha_diff_min = chi_list_dense[np.argmin(alphas_diff_dense/taper_angles_dense,axis=1)]
chis_alpha_diff_max = chi_list_dense[np.argmax(alphas_diff_dense/taper_angles_dense,axis=1)]

bDalpha = alphas_diff_dense/taper_angles_dense
dum = (bDalpha[1:,:]-bDalpha[0:-1,:])/((LambdaRef_list_dense[1]-LambdaRef_list_dense[0])*100.0)
dAlpha_dLambda = (dum[:,1:] + dum[:,0:-1])/2.0
dum = (bDalpha[:,1:]-bDalpha[:,0:-1])/((chi_list_dense[1]-chi_list_dense[0])*100.0)
dAlpha_dChi = (dum[1:,:] + dum[0:-1,:])/2.0
vDiv = dAlpha_dLambda+dAlpha_dChi

chis_dense,Lambdas_dense = np.meshgrid(chi_list_dense,LambdaRef_list_dense)


dum = (Lambdas_dense[0:-1,:] + Lambdas_dense[1:,:])/2.0
Lambdas_centered = (dum[:,0:-1] + dum[:,1:])/2.0
dum = (chis_dense[0:-1,:] + chis_dense[1:,:])/2.0
chis_centered = (dum[:,0:-1] + dum[:,1:])/2.0

chi_list_centered = (chi_list_dense[0:-1] + chi_list_dense[1:])/2.0
LambdaRef_list_centered = (LambdaRef_list_dense[0:-1] + LambdaRef_list_dense[1:])/2.0
chis_vDiv_0 = chi_list_centered [np.argmin(vDiv**2,axis=1)]
chis_vDiv_0 = np.concatenate([chis_vDiv_0,[0.0]])
chis_vDiv_0[-13:-1] = 0.0

#chis_alpha_WB_up_max = chi_list_dense[np.argmax(alphas_WB_up_dense,axis=1)]
#chis_vDiv_0 = chis_alpha_WB_up_max#chi_list_centered [np.argmin(vDiv,axis=1)]
#chis_vDiv_0[-1] = 0.0#chis_vDiv_0[-2]


width = 4
chis_vDiv_0_old = chis_vDiv_0.copy()
for i in range (0,len(chis_vDiv_0)):
    if i<=width:
        thisWidth = i
    elif i>len(chis_vDiv_0)-width-1:
        thisWidth = len(chis_vDiv_0)-1-i
    else:
        thisWidth = width
    sumVal = 0
    for j in range(-thisWidth,thisWidth+1):
        sumVal += chis_vDiv_0_old[i+j]
        
    chis_vDiv_0[i] = sumVal/(2.0*thisWidth+1.0)



#f = interpolate.interp1d(LambdaRef_list, chis_vDiv_0)
#chis_vDiv_0_dense = f(LambdaRef_list_dense)



# Find types
# ============================================

Type = np.zeros([denseFac*nLambda,denseFac*nChi])
#chi_boundary = 
for iL in range(denseFac*nLambda):
    chi_boundary = chis_vDiv_0[iL]
    for iC in range(denseFac*nChi):
        chi = chi_list_dense[iC]
        if chi<chi_boundary:
            Type[iL,iC] = 1
        else:
            if alphas_diff_dense[iL][iC]<0.0:
                Type[iL,iC]=3
            else:
                Type[iL,iC]=2
    #end iC
#end iL
    
#   Find boundaries
# ============================================
bound_10 = []
bound_12 = []
bound_21 = []
bound_23 = []
bound_32 = []
bound_30 = []
for iC in range(denseFac*nChi-1):
    y = chi_list_dense[iC]
    for iL in range(denseFac*nLambda-1):
        x = LambdaRef_list_dense[iL]
        if Type[iL,iC]==1 and (iC==0):
            bound_10.append([x,y])
        elif Type[iL,iC]==2 and (Type[iL-1,iC]==1 or Type[iL,iC-1]==1):
            bound_21.append([x,y])
        elif Type[iL,iC]==3 and (Type[iL-1,iC]==2 or Type[iL,iC-1]==2):
            bound_32.append([x,y])
            
        if Type[iL,iC]==1 and (Type[iL+1,iC]==2 or Type[iL,iC+1]==2):
            bound_12.append([x,y])
        elif Type[iL,iC]==2 and (Type[iL,iC+1]==3):
            bound_23.append([x,y])        
        elif Type[iL,iC]==3 and (iC==denseFac*nChi-2):
            bound_30.append([x,y])
    #end iC
#end iL
    
bound_10 = arr(bound_10)
bound_12 = arr(bound_12)
bound_21 = arr(bound_21)
bound_23 = arr(bound_23)
bound_32 = arr(bound_32)
bound_30 = arr(bound_30)
    
    
    

#from matplotlib.colors import LinearSegmentedColormap
plt.sca(Axes['12'])
plt.cla()
#plt.contourf(Lambdas_centered*100.0,chis_centered*100.0,Type,vmin=0,vmax=7)
CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(0.000,1e3,2),vmin=-1.00,vmax=1.00)
plt.cla()
zeroContour = CS.allsegs[0][0]

deleteIndex = []
for i in range(zeroContour.shape[0]):
#    for j in zeroContour.shape[1]:
        if zeroContour[i,0] < 0.1 or zeroContour[i,1] < 2.0 or zeroContour[i,1] > 98.5:
            deleteIndex.append(i)
            
zeroContour = np.delete(zeroContour,deleteIndex,0)




bound_12 = arr([LambdaRef_list_dense,chis_vDiv_0]).T
bound_12 = np.concatenate([arr([[0.0,chis_vDiv_0[0]]]),bound_12])
bound_23 = np.flipud(zeroContour)/100.0
bound_23 = np.concatenate([arr([[1.0,0.0]]),bound_23,arr([[0.0,1.0]])])


domain1 = np.concatenate([arr([[0.0,0.0]]),bound_12])
domain2 = np.concatenate([bound_12,bound_23])
domain3 = np.concatenate([bound_23,arr([[1.0,1.0]])])


CMAP = [[.99,0.0,0.0],
        [0.0,1.0,0.0],
        [0.0,0.0,1.0]]

colors=[(0.75, 0.15, 0.15), (1, 0.75, 0.15), (0.15, 0.75, 0.15)]
#plt.set_cmap("Set2")

                


#   Attribute type
# ============================================
floatType = np.zeros([denseFac*nLambda,denseFac*nChi]) 
for iC in range(denseFac*nChi):
    y = chi_list_dense[iC]
    for iL in range(denseFac*nLambda):
        x = LambdaRef_list_dense[iL]
        if Type[iL,iC] == 1:
            Dist_low = np.min(np.sqrt( (x-bound_12[:,0])**2 + (y-bound_12[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_10[:,0])**2 + (y-bound_10[:,1])**2))

            
        elif Type[iL,iC] == 2:
            Dist_low = np.min(np.sqrt( (x-bound_23[:,0])**2 + (y-bound_23[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_21[:,0])**2 + (y-bound_21[:,1])**2))

            
        elif Type[iL,iC] == 3:
            Dist_low = np.min(np.sqrt( (x-bound_30[:,0])**2 + (y-bound_30[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_32[:,0])**2 + (y-bound_32[:,1])**2))
            
            
            
        if (Dist_low+Dist_up)!=0.0:
            floatType[iL,iC] = Type[iL,iC]+1 - Dist_low/(Dist_low+Dist_up)
        else:
            floatType[iL,iC] = Type[iL,iC]+1
            
            
            
#   Plot
# ============================================
plt.sca(Axes['12'])
plt.cla()
#plt.contourf(LambdaRef_list_dense*100.0,chi_list_dense*100.0,floatType.T,np.linspace(1.0,4.0,19),vmin=1.0,vmax=4.0)
plt.contourf(LambdaRef_list_dense*100.0,chi_list_dense*100.0,floatType.T,np.linspace(1.0,4.0,128),vmin=1.0,vmax=4.0)
#plt.contourf(LambdaRef_list_dense*100.0,chi_list_dense*100.0,floatType.T)

plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'-k',linewidth=1.5)
plt.sca(Axes['11'])
plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'--k',linewidth=0.5)
plt.sca(Axes['12'])




#   Colormap
# ============================================
CMAP, color_list = Style.getCmap_Type()
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")





#   Draw some text
# ============================================
plt.sca(Axes['11'])
plt.text(43,20,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(50,50,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(70,85,'III',family='Times New Roman',color='w',size=28,weight='bold')
plt.sca(Axes['12'])
plt.text(43,20,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(50,50,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(70,85,'III',family='Times New Roman',color='w',size=28,weight='bold')
#   axis work
# ============================================
plt.sca(Axes['12'])
ax = plt.gca()
ax.axes.get_yaxis().set_ticklabels([])
plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
plt.axis([40.0,100.0,100.0,.0])

#   Colorbar stuff
# ============================================
plt.sca(Axes['11'])
cbar = plt.colorbar(cax=cBarAxes['11'], ticks=[-1, 0, 1])

plt.sca(cBarAxes['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)

plt.sca(Axes['12'])
cbar = plt.colorbar(cax=cBarAxes2['11'], ticks=[1, 2, 3, 4])
plt.sca(cBarAxes2['11'])
plt.text(0.5,1.05,"$\\mathbf{Type}}$",horizontalAlignment='center',fontdict=Style.fontdict)

#cbar.set_ticks()
    
    
#   Save stuff
# ============================================  
np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0),
     Lambdas = LambdaRef_list_dense,
     chis    = chi_list_dense,
     floatType = floatType,
     alphas_diff = alphas_diff_dense,
     taper_angles = taper_angles_dense)