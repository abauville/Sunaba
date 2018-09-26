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

#(nChi, nBeta, nLambda, LambdaRef_list, 
# chi_list, betas_all, alphas_Ref_all, 
# alphas_WF_all, alphas_WB_up_all, alphas_WB_low_all, Lambdas_Ref_all, chis_all, 
# Taper_Ref, Taper_WB, Taper_WF) = CritTaper_dataMaker.getCritTaperFigData(Compute=False, beta_list=0.0, nChi=51, nLambda=51,enveloppeRes=6001,alphaMin=-1.0*np.pi/180.0)

alphas_diff_all = alphas_WB_up_all - alphas_Ref_all

Style = CritTaper_Style.Style()
plt.set_cmap(Style.colormap)
## Lambda vs chi @ beta=0
fig    = Figz_Utils.Figure(77,height=13.0,mode='production')
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
beta = 0.0
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
plt.text(75,90,"Extensional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)
plt.text(35,20,"Compressional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)

plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
plt.ylabel("$\\mathbf{\\chi}$ [%]",weight='bold',verticalAlignment='center')

ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])






#   Define interpolated versions of arrays
# ============================================
denseFac = 3
chi_list_dense = np.linspace(chi_list[0],chi_list[-1],denseFac*nChi)
LambdaRef_list_dense = np.linspace(LambdaRef_list[0],LambdaRef_list[-1],denseFac*nLambda)

from scipy import interpolate
f = interpolate.RectBivariateSpline(chi_list,LambdaRef_list,alphas_diff)
alphas_diff_dense = f(chi_list_dense,LambdaRef_list_dense)

f = interpolate.RectBivariateSpline(chi_list,LambdaRef_list,taper_angles)
taper_angles_dense = f(chi_list_dense,LambdaRef_list_dense)

f = interpolate.RectBivariateSpline(chi_list,LambdaRef_list,alphas_WB_up)
alphas_WB_up_dense = f(chi_list_dense,LambdaRef_list_dense)

#f = interpolate.interp1d(LambdaRef_list, chis_vGrad_min)
#chis_vGrad_min_dense = f(LambdaRef_list_dense)




# Add indication of max delta alpha
plt.sca(Axes['11'])
chis_alpha_diff_min = chi_list_dense[np.argmin(alphas_diff_dense/taper_angles_dense,axis=1)]
chis_alpha_diff_max = chi_list_dense[np.argmax(alphas_diff_dense/taper_angles_dense,axis=1)]


chis_alpha_WB_up_max = chi_list_dense[np.argmax(alphas_WB_up_dense,axis=1)]


chis_vGrad_min = chis_alpha_WB_up_max#chi_list_centered [np.argmin(vGrad,axis=1)]
chis_vGrad_min[-1] = 0.0#chis_vGrad_min[-2]
width = 4
chis_vGrad_min_old = chis_vGrad_min.copy()
for i in range (0,len(chis_vGrad_min)):
    if i<=width:
        thisWidth = i
    elif i>len(chis_vGrad_min)-width-1:
        thisWidth = len(chis_vGrad_min)-1-i
    else:
        thisWidth = width
    sumVal = 0
    for j in range(-thisWidth,thisWidth+1):
        sumVal += chis_vGrad_min_old[i+j]
        
    chis_vGrad_min[i] = sumVal/(2.0*thisWidth+1.0)



#f = interpolate.interp1d(LambdaRef_list, chis_vGrad_min)
#chis_vGrad_min_dense = f(LambdaRef_list_dense)



# Find types
# ============================================

Type = np.zeros([denseFac*nLambda,denseFac*nChi])
#chi_boundary = 
for iL in range(denseFac*nLambda):
    chi_boundary = chis_vGrad_min[iL]
    for iC in range(denseFac*nChi):
        chi = chi_list_dense[iC]
        if chi<chi_boundary:
            Type[iC][iL] = 1
        else:
            if alphas_diff_dense[iL][iC]<0.0:
                Type[iC][iL]=3
            else:
                Type[iC][iL]=2
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
for iC in range(denseFac*nChi):
    y = chi_list_dense[iC]
    for iL in range(denseFac*nLambda):
        x = LambdaRef_list_dense[iL]
        if Type[iC,iL]==1 and (iC==0):
            bound_10.append([x,y])
        elif Type[iC,iL]==2 and (Type[iC-1,iL]==1 or Type[iC,iL-1]==1):
            bound_21.append([x,y])
        elif Type[iC,iL]==3 and (Type[iC-1,iL]==2 or Type[iC,iL-1]==2):
            bound_32.append([x,y])
            
        if Type[iC,iL]==1 and (Type[iC+1,iL]==2 or Type[iC,iL+1]==2):
            bound_12.append([x,y])
        elif Type[iC,iL]==2 and (Type[iC,iL+1]==3):
            bound_23.append([x,y])        
        elif Type[iC,iL]==3 and (iC==denseFac*nChi-2):
            bound_30.append([x,y])
    #end iC
#end iL
    
bound_10 = arr(bound_10)
bound_12 = arr(bound_12)
bound_21 = arr(bound_21)
bound_23 = arr(bound_23)
bound_32 = arr(bound_32)
bound_30 = arr(bound_30)
    
    
    

from matplotlib.colors import LinearSegmentedColormap
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




bound_12 = arr([LambdaRef_list_dense,chis_vGrad_min]).T
bound_12 = np.concatenate([arr([[0.0,chis_vGrad_min[0]]]),bound_12])
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
floatType = np.zeros([denseFac*nChi,denseFac*nLambda]) 
for iC in range(denseFac*nChi):
    y = chi_list_dense[iC]
    for iL in range(denseFac*nLambda):
        x = LambdaRef_list_dense[iL]
        if Type[iC,iL] == 1:
            Dist_low = np.min(np.sqrt( (x-bound_12[:,0])**2 + (y-bound_12[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_10[:,0])**2 + (y-bound_10[:,1])**2))

            
        elif Type[iC,iL] == 2:
            Dist_low = np.min(np.sqrt( (x-bound_23[:,0])**2 + (y-bound_23[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_21[:,0])**2 + (y-bound_21[:,1])**2))

            
        elif Type[iC,iL] == 3:
            Dist_low = np.min(np.sqrt( (x-bound_30[:,0])**2 + (y-bound_30[:,1])**2))
            Dist_up = np.min(np.sqrt( (x-bound_32[:,0])**2 + (y-bound_32[:,1])**2))
            
            
            
        if (Dist_low+Dist_up)!=0.0:
            floatType[iC,iL] = -Type[iC,iL]+2 + Dist_low/(Dist_low+Dist_up)
        else:
            floatType[iC,iL] = -Type[iC,iL]+2
            
            
            
#   Plot
# ============================================
plt.sca(Axes['12'])
plt.cla()
plt.contourf(LambdaRef_list_dense*100.0,chi_list_dense*100.0,floatType,np.linspace(-1.000,2.0,130),vmin=-1.0,vmax=2.0)

plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'-k',linewidth=1.5)
plt.sca(Axes['11'])
plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'--k',linewidth=0.5)
plt.sca(Axes['12'])




#   Colormap
# ============================================
#colors = arr([[200, 30, 32],
#              [248,120, 64],
#              [248,180, 64],
#              [ 64,180,248],
#              [ 64,120,248],
#              [ 32, 30,200]]) / 255.0
colors = arr([[200, 30, 32],
              [248,160, 32],
              [248,200,  0],
              [  0,200,248],
              [ 32,160,248],
              [ 32, 30,200]]) / 255.0
colors = np.flipud(colors)
nSeg = colors.shape[0]-1

# algorithm to blend colors
def blendColorValue(a, b, t):
    return np.sqrt((1 - t) * a**2 + t * b**2)

C = 0
nSubSegs = np.array([40,1,40,1,40])
nTot = np.int(np.sum(nSubSegs)) + 1
blendedColors = np.zeros([nTot,3])
i = 0
for iSeg in range(nSeg):
    iSub = 0
    for t in np.linspace(0.0,1.0,nSubSegs[iSeg]+1):
        if (iSeg<nSeg-1 and iSub==nSubSegs[iSeg]):
            break
        else:
            blendedColors[i,:] = blendColorValue(colors[iSeg,:],colors[iSeg+1,:],t)
            i+=1
            iSub+=1
        

CMAP = LinearSegmentedColormap.from_list('custom',blendedColors,N=nTot)

plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")





#   Draw some text
# ============================================
plt.text(10,25,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(25,70,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(60,90,'III',family='Times New Roman',color='w',size=28,weight='bold')

#   axis work
# ============================================
ax = plt.gca()
ax.axes.get_yaxis().set_ticklabels([])
plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])

#   Colorbar stuff
# ============================================
plt.sca(Axes['11'])
cbar = plt.colorbar(cax=cBarAxes['11'], ticks=[-1, 0, 1])

plt.sca(cBarAxes['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)

plt.sca(Axes['12'])
cbar = plt.colorbar(cax=cBarAxes2['11'], ticks=[-1, 0, 1, 2])
plt.sca(cBarAxes2['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)

#cbar.set_ticks()