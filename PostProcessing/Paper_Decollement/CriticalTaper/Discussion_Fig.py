#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 13:42:52 2018

@author: abauville
"""
import sys
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '../../src/UserInput')
sys.path.insert(0, './CriticalTaper')
sys.path.insert(0, '../../Utils/')
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
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=np.linspace(0.0,30.0,13.0)*np.pi/180.0, nChi=61, nLambda=61,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)


# Limits of the colorbar
c0 = -1.0
c1 = 1.0


Lambda_hydro = 0.4
Lambda_ovRef_list = (LambdaRef_list-Lambda_hydro)/(1.0-Lambda_hydro)


beta = 0.0*np.pi/180.0
alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
## Lambda vs chi @ beta=0
fig    = Figz_Utils.Figure(5,height=10.0,mode='crop')
#fig    = Figz_Utils.Figure(7,height=13.0,mode='draft')
#AxesDum   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5)
#AxesDum['12'].axis('off')
#Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=10.5,topMarginPad = 1.0)
xPad = 1.0
Width = 4.75
xPad1 = 1.5
xPad2 = 0.5
Axes   = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin)-Width-1.5,topMarginPad = 1.2,xPad = xPad1)
Axes2  = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5+1.0*Width+xPad1,rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin)-2.0*Width-1.5-xPad1,topMarginPad = 1.2,xPad = xPad2)
Axes3  = Figz_Utils.makeAxes(fig,1,1,aspectRatio=1.0,leftMarginPad=1.5+2.0*Width+xPad1+xPad2,rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin)-3.0*Width-1.5-xPad1-xPad2,topMarginPad = 1.2,xPad = xPad2)

Axes['12'] = Axes2['11']
Axes['13'] = Axes3['11']

AxesW = Axes['info']['plotsWidth']
AxesH = Axes['info']['plotsHeight']
AxesxPad = Axes['info']['xPad']
AxeslPad = Axes['info']['leftMarginPad']
AxesrPad = Axes['info']['rightMarginPad']
AxestPad = Axes['info']['topMarginPad']
AxesbPad = Axes['info']['bottomMarginPad']
cBaryPad    = 0.0
cBartPad = .2
cBarbPad = 0.0
cBarlPad = 0.0
cBarrPad = 0.0
cBarW = 0.0
cBarH = 0.35
#cBarAxes   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+AxesW+cBarlPad,
#                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+AxesW+cBarlPad)-cBarW),
#                                       topMarginPad=AxestPad+cBartPad,
#                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)
#cBarAxes2   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+2.0*AxesW+cBarlPad+xPad,
#                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+2.0*AxesW+cBarlPad+xPad)-cBarW),
#                                       topMarginPad=AxestPad+cBartPad,
#                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)



cBarAxes   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+AxesW+xPad1,
                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+2.0*AxesW+cBarlPad)-cBarW)-xPad1,
                                       topMarginPad=AxestPad+cBartPad+AxesH,
                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad-cBarH-cBartPad)
cBarAxes2  = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+2.0*AxesW+xPad1+xPad2,
                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+3.0*AxesW+cBarlPad)-cBarW)-xPad1-xPad2,
                                       topMarginPad=AxestPad+cBartPad+AxesH,
                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad-cBarH-cBartPad)
#cBarAxes2   = Figz_Utils.makeAxes(fig,1,leftMarginPad=AxeslPad+2.0*AxesW+cBarlPad+xPad,
#                                       rightMarginPad=(fig.width-fig.rightMargin-fig.leftMargin-(AxeslPad+2.0*AxesW+cBarlPad+xPad)-cBarW),
#                                       topMarginPad=AxestPad+cBartPad,
#                                       bottomMarginPad=(fig.height-fig.topMargin-fig.bottomMargin-AxestPad-AxesH)+cBarbPad)

#Axes['12'].axis('off')
plt.sca(Axes['12'])



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
plt.sca(Axes['12'])
phiRef = 30.0*np.pi/180.0
Lambda_hydro = 0.4
Lambda_ovs = (Lambdas-Lambda_hydro)/(1-Lambda_hydro)
alphas_repose = np.arctan((1.0-Lambda_ovs)*np.tan(phiRef))
#CS = plt.contourf(Lambda_ovs*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(c0-0.0001,c1+.0001,20),vmin=c0,vmax=c1)
CS = plt.contourf(Lambda_ovs*100.0, chis*100.0, alphas_diff/alphas_repose,np.linspace(c0-0.0001,c1+.0001,20),vmin=c0,vmax=c1)
#plt.text(75,90,"Extensional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)
#plt.text(35,20,"Compressional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)

plt.xlabel("$\\mathbf{\\lambda^*}$ [%]",weight='bold',verticalAlignment='center')
plt.ylabel("$\\mathbf{\\chi}$ [%]",weight='bold',verticalAlignment='top')

plt.xticks([0,33,66,100],[0,33,66,100])
plt.yticks([0,33,66,100],[0,33,66,100])

ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.axis([0.0,100.0,100.0,.0])






#   Define interpolated versions of arrays
# ============================================
denseFac = 3
chi_list_dense = np.linspace(chi_list[0],chi_list[-1],denseFac*nChi)
Lambda_ovRef_list_dense = np.linspace(Lambda_ovRef_list[0],Lambda_ovRef_list[-1],denseFac*nLambda)

from scipy import interpolate
f = interpolate.RectBivariateSpline(Lambda_ovRef_list,chi_list,alphas_diff)
alphas_diff_dense = f(Lambda_ovRef_list_dense,chi_list_dense)

f = interpolate.RectBivariateSpline(Lambda_ovRef_list,chi_list,taper_angles)
taper_angles_dense = f(Lambda_ovRef_list_dense,chi_list_dense)

f = interpolate.RectBivariateSpline(Lambda_ovRef_list,chi_list,alphas_WB_up)
alphas_WB_up_dense = f(Lambda_ovRef_list_dense,chi_list_dense)

#f = interpolate.interp1d(Lambda_ovRef_list, chis_vGrad_0)
#chis_vGrad_0_dense = f(Lambda_ovRef_list_dense)




# Add indication of max delta alpha
plt.sca(Axes['12'])
#chis_alpha_diff_min = chi_list_dense[np.argmin(alphas_diff_dense/taper_angles_dense,axis=1)]
#chis_alpha_diff_max = chi_list_dense[np.argmax(alphas_diff_dense/taper_angles_dense,axis=1)]
#
#
#chis_alpha_WB_up_max = chi_list_dense[np.argmax(alphas_WB_up_dense,axis=1)]
#
#
#chis_vGrad_0 = chis_alpha_WB_up_max#chi_list_centered [np.argmin(vGrad,axis=1)]
#chis_vGrad_0[-1] = 0.0#chis_vGrad_0[-2]
#
#


chis_alpha_diff_min = chi_list_dense[np.argmin(alphas_diff_dense/taper_angles_dense,axis=1)]
chis_alpha_diff_max = chi_list_dense[np.argmax(alphas_diff_dense/taper_angles_dense,axis=1)]

bDalpha = alphas_diff_dense/taper_angles_dense
dum = (bDalpha[1:,:]-bDalpha[0:-1,:])/((Lambda_ovRef_list_dense[1]-Lambda_ovRef_list_dense[0])*100.0)
dAlpha_dLambda = (dum[:,1:] + dum[:,0:-1])/2.0
dum = (bDalpha[:,1:]-bDalpha[:,0:-1])/((chi_list_dense[1]-chi_list_dense[0])*100.0)
dAlpha_dChi = (dum[1:,:] + dum[0:-1,:])/2.0
vDiv = dAlpha_dLambda+dAlpha_dChi
vGrad = np.sqrt(dAlpha_dLambda**2+dAlpha_dChi**2)

chis_dense,Lambdas_dense = np.meshgrid(chi_list_dense,Lambda_ovRef_list_dense)


dum = (Lambdas_dense[0:-1,:] + Lambdas_dense[1:,:])/2.0
Lambdas_centered = (dum[:,0:-1] + dum[:,1:])/2.0
dum = (chis_dense[0:-1,:] + chis_dense[1:,:])/2.0
chis_centered = (dum[:,0:-1] + dum[:,1:])/2.0

chi_list_centered = (chi_list_dense[0:-1] + chi_list_dense[1:])/2.0
Lambda_ovRef_list_centered = (Lambda_ovRef_list_dense[0:-1] + Lambda_ovRef_list_dense[1:])/2.0
chis_vGrad_0 = chi_list_centered [np.argmin(vGrad**2,axis=1)]
chis_vGrad_0 = np.concatenate([chis_vGrad_0,[0.0]])
chis_vGrad_0[-3:-1] = 0.0

#chis_alpha_WB_up_max = chi_list_dense[np.argmax(alphas_WB_up_dense,axis=1)]
#chis_vGrad_0 = chis_alpha_WB_up_max#chi_list_centered [np.argmin(vGrad,axis=1)]
#chis_vGrad_0[-1] = 0.0#chis_vGrad_0[-2]


width = 5

for t in range(2):
    chis_vGrad_0_old = chis_vGrad_0.copy()
    for i in range (0,len(chis_vGrad_0)):
        if i<=width:
            thisWidth = i
        elif i>len(chis_vGrad_0)-width-1:
            thisWidth = len(chis_vGrad_0)-1-i
        else:
            thisWidth = width
        sumVal = 0
        for j in range(-thisWidth,thisWidth+1):
            sumVal += chis_vGrad_0_old[i+j]
            
        chis_vGrad_0[i] = sumVal/(2.0*thisWidth+1.0)
    


#f = interpolate.interp1d(Lambda_ovRef_list, chis_vGrad_0)
#chis_vGrad_0_dense = f(Lambda_ovRef_list_dense)



# Find types
# ============================================

Type = np.zeros([denseFac*nLambda,denseFac*nChi])
#chi_boundary = 
for iL in range(denseFac*nLambda):
    chi_boundary = chis_vGrad_0[iL]
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
        x = Lambda_ovRef_list_dense[iL]
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
        elif Type[iL,iC]==3 and (iC==denseFac*nChi-2) and iL==denseFac*nLambda-2:
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
plt.sca(Axes['13'])
plt.cla()
#plt.contourf(Lambdas_centered*100.0,chis_centered*100.0,Type,vmin=0,vmax=7)
CS = plt.contourf(Lambda_ovs*100.0, chis*100.0, alphas_diff,np.linspace(0.000,1e3,2),vmin=-1.00,vmax=1.00)
plt.cla()
zeroContour = CS.allsegs[0][0]

deleteIndex = []
for i in range(zeroContour.shape[0]):
#    for j in zeroContour.shape[1]:
        if zeroContour[i,0] < 0.1 or zeroContour[i,1] < 2.0 or zeroContour[i,1] > 98.5:
            deleteIndex.append(i)
            
zeroContour = np.delete(zeroContour,deleteIndex,0)




bound_12 = arr([Lambda_ovRef_list_dense,chis_vGrad_0]).T
bound_12 = np.concatenate([arr([[0.0,chis_vGrad_0[0]]]),bound_12])
bound_23 = np.flipud(zeroContour)/100.0
bound_23 = np.concatenate([arr([[1.0,0.0]]),bound_23])


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
        x = Lambda_ovRef_list_dense[iL]
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
plt.sca(Axes['13'])
plt.cla()
#plt.contourf(Lambda_ovRef_list_dense*100.0,chi_list_dense*100.0,floatType.T,np.linspace(1.0,4.0,19),vmin=1.0,vmax=4.0)
plt.contourf(Lambda_ovRef_list_dense*100.0,chi_list_dense*100.0,floatType.T,np.linspace(1.0,4.0,128),vmin=1.0,vmax=4.0)
#plt.contourf(Lambda_ovRef_list_dense*100.0,chi_list_dense*100.0,floatType.T)

plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'-k',linewidth=1.5)
plt.sca(Axes['12'])
plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'--k',linewidth=0.5)
plt.sca(Axes['13'])




#   Colormap
# ============================================
CMAP, color_list = Style.getCmap_Type()
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")





#   Draw some text
# ============================================
plt.sca(Axes['12'])
plt.text(5,20,'I',family='Times New Roman',color='w',size=22,weight='bold')
plt.text(15,58,'II',family='Times New Roman',color='w',size=22,weight='bold')
plt.text(30,85,'III',family='Times New Roman',color='w',size=22,weight='bold')
plt.sca(Axes['13'])
plt.text(5,20,'I',family='Times New Roman',color='w',size=22,weight='bold')
plt.text(15,58,'II',family='Times New Roman',color='w',size=22,weight='bold')
plt.text(30,85,'III',family='Times New Roman',color='w',size=22,weight='bold')
#   axis work
# ============================================
plt.sca(Axes['13'])
ax = plt.gca()
ax.axes.get_yaxis().set_ticklabels([])
plt.xlabel("$\\mathbf{\\lambda^*}$ [%]",weight='bold',verticalAlignment='center')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
plt.axis([00.0,100.0,100.0,.0])
plt.xticks([0,33,66,100])


#   Colorbar stuff
# ============================================
plt.sca(Axes['12'])
cbar = plt.colorbar(cax=cBarAxes['11'], ticks=[c0, 0, c1], orientation='horizontal')

plt.sca(cBarAxes['11'])
plt.text(0.5,-2.5,"$\\mathbf{\\Delta \\bar{\\alpha}} [°]$",horizontalAlignment='center',fontdict=Style.fontdict)

plt.sca(Axes['13'])
cbar = plt.colorbar(cax=cBarAxes2['11'], ticks=[1, 2, 3, 4], orientation='horizontal')
plt.sca(cBarAxes2['11'])
plt.text(0.5,-2.5,"$\\mathbf{M}}$",horizontalAlignment='center',fontdict=Style.fontdict)

#cbar.set_ticks()


plt.sca(Axes['11'])
plt.cla()
#   Load stuff
# ============================================  
loadedData  = np.load("/Users/abauville/Output/Paper_Decollement/Figz/Data/alpha_phi_b.npz")
alpha_up    = loadedData["alpha_up"][()]
alpha_low   = loadedData["alpha_low"][()]
phi_b       = loadedData["phi_b"][()]
Lambda_Ref  = loadedData["LambdaRef"][()]
phiRef      = loadedData["phiRef"][()]
    
## Plotting  
#alpha_repose = (1.0-Lambda_ov)*np.tan(phiRef)
Lambda_ov=0.0
alpha_repose = np.arctan((1.0-Lambda_ov)*np.tan(phiRef))
alpha_p_repose = alpha_repose/(1.0-Lambda_ov)      
#alpha_up = alphas_WB_up[0,:]       
#alpha_low= alphas_WB_low[0,:]      
            
alpha_up = alpha_up[:,0]
alpha_low = alpha_low[:,0]

import CritTaper_Style
Style = CritTaper_Style.Style()


plotLambda_ov = 0.0
#I = np.argmin(np.abs(plotLambda_ov-Lambda_ov))
#plotLambda_ov = Lambda_ov[I]
plotAlpha_repose = np.arctan((1.0-plotLambda_ov)*np.tan(phiRef))

chi = np.linspace(1,0,len(alpha_up))


plt.fill(np.concatenate([(1.0-chi),chi]) , np.concatenate([alpha_up,np.flipud(alpha_low)])/alpha_repose,color=Style.colorBW_a)
plt.plot((1.0-chi),alpha_up/alpha_repose ,'-',color=Style.colorBW,linewidth=1.0)
plt.plot((1.0-chi),alpha_low/alpha_repose,'--',color=Style.colorBW,linewidth=1.0)

plt.plot([.0,1.0],[alpha_up[-1]/alpha_repose,alpha_up[-1]/alpha_repose],'-',color=Style.colorRef,linewidth=0.75)
I = np.argmin(np.abs(alpha_up-alpha_repose))
plt.plot(arr([I,I])/len(alpha_up),[.0,1.0],'--k',linewidth=0.5)
I = np.argmin(np.abs(alpha_up[:-1]-alpha_up[-1]))
plt.plot(arr([I,I])/len(alpha_up),[.0,1.0],'--k',linewidth=0.5)


plt.text(.05,.85,'III',family='Times New Roman',color='k',size=18,weight='bold')
plt.text(.27,.85,'II' ,family='Times New Roman',color='k',size=18,weight='bold')
plt.text(.67,.85,'I'  ,family='Times New Roman',color='k',size=18,weight='bold')

plt.xlabel('$\\mathbf{1-\\chi=\\mu_b/\\mu}$')
plt.ylabel('$\mathbf{\\bar{\\alpha}=\\alpha/\\alpha_{max}}$')

plt.xticks([.0,.5,1.0])
plt.yticks([.0,.5,1.0])
plt.axis([.0,1.0,.0,1.0])

#Letters = 'ABCD'
#x0s = [.025,2.5,2.5]
#y0s = [.025,100.0-2.5,100-2.5]
#colors = ['k','w','w']
#for i in range(3):
#    plt.sca(Axes['1%i' % (i+1)])
#    ax = plt.gca()
#    ax.text(x0s[i],y0s[i],"%s" % (Letters[i]),color=colors[i],fontdict=Style.fontdict,horizontalAlignment='left',verticalAlignment='baseline',size=12)
#    
#   Save stuff
# ============================================  
np.savez("/Users/abauville/Output/Paper_Decollement/Figz/Data/floatType_beta%02d.npz" % round(beta*180.0/np.pi*10.0),
     Lambdas = Lambda_ovRef_list_dense,
     chis    = chi_list_dense,
     floatType = floatType,
     alphas_diff = alphas_diff_dense,
     taper_angles = taper_angles_dense)
#    
#    
#fig2    = Figz_Utils.Figure(500,height=13.0,mode='production')
#Axesb   = Figz_Utils.makeAxes(fig2,1,2,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=1.5,topMarginPad = 1.0,xPad = xPad)
#
#plt.sca(Axesb['11'])
#plt.contourf(Lambda_ovs*100.0,chis*100.0,alphas_diff/(1-Lambda_ovs)*deg,30)
##plt.contourf(Lambda_ovs*100.0,chis*100.0,alphas_diff*deg,30)
##plt.contourf(Lambda_ovs*100.0,chis*100.0,alphas_WB_up*deg,30)
#plt.axis([0.0,100.0,100.0,0.0])
#plt.colorbar()
#
#plt.sca(Axesb['12'])
##plt.contourf(Lambda_ovs*100.0,chis*100.0,alphas_Ref/(1-Lambda_ovs)*deg,30)
#plt.contourf(Lambda_ovs*100.0,chis*100.0,alphas_WB_up/(1-Lambda_ovs)*deg,30)
##plt.contourf(Lambda_ovs*100.0,chis*100.0,( alphas_diff/(1-Lambda_ovs) )/( alphas_Ref/(1-Lambda_ovs) - beta),30)
#plt.axis([0.0,100.0,100.0,0.0])
#plt.colorbar()