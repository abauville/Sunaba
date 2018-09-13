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
Axes   = Figz_Utils.makeAxes(fig,1,2,aspectRatio=1.0,leftMarginPad=1.5,rightMarginPad=1.5,topMarginPad = 1.0,xPad = 3.0)
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
#Axes['12'].axis('off')
plt.sca(Axes['11'])



deg = 180.00/np.pi


Lambdas = Lambdas_Ref_all[:,:,0]
chis = chis_all[:,:,0]
#alphas_diff = alphas_diff_all[iTaper,:,:]
#Lambdas = np.zeros((nLambda,nChi))
alphas_diff = np.zeros((nLambda,nChi))
alphas_width = np.zeros((nLambda,nChi))
#chis = np.zeros((nLambda,nChi))
taper_angles = np.zeros((nLambda,nChi))
alphas_WB_up = np.zeros((nLambda,nChi))
beta = 0.0
for iL in range(nLambda):
    for iW in range(nChi):
        iB = np.argmin(abs(betas_all[iL,iW,:]-beta))
        alphas_diff[iL,iW] = alphas_diff_all[iL,iW,iB]
        alphas_width[iL,iW] = alphas_WB_up_all[iL,iW,iB] - alphas_WB_low_all[iL,iW,iB]
        alphas_WB_up[iL,iW] = alphas_WB_up_all[iL,iW,iB]
        taper_angles[iL,iW] = betas_all[iL,iW,iB]+alphas_Ref_all[iL,iW,iB]
        
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,1000)

CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(-1.0001,1.0001,20),vmin=-1.00,vmax=1.00)
#CS = plt.contourf(Lambdas*100.0, chis*100.0, alphas_diff*deg,20,vmin=-20.0,vmax=20.0)
#plt.pcolormesh(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,shading='Gouraud',vmin=-1.0,vmax=1.0)

# = np.zeros((nLambda,nChi))

#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff-taper_angles, levels = [-1000.0, 0.0, 1000.0])
#ax.clabel(CS, CS.levels, inline=True, fontsize=16)

plt.text(75,90,"Extensional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)
plt.text(35,20,"Compressional",fontdict=Style.fontdict,horizontalAlignment="center",verticalAlignment="center",color="w",size=13)

plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
plt.ylabel("$\\mathbf{\\chi}$ [%]",weight='bold',verticalAlignment='center')

#plt.plot(alphas_width*180.0/pi)
#


ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])






#cbar = plt.colorbar()
#cbar.set_ticks([-1.0,0.0,1.0])
plt.sca(Axes['11'])
#plt.axis([-10.0,70.0,0.0,1.0])
cbar = plt.colorbar(cax=cBarAxes['11'], ticks=[-1, 0, 1])
#cbar = plt.colorbar(cax=cBarAxes['11'])

plt.sca(cBarAxes['11'])
plt.text(0.5,1.05,"$\\mathbf{\\bar{\\Delta \\alpha}}$",horizontalAlignment='center',fontdict=Style.fontdict)



#cbar.ax.set_xticklabels([-1.0,"",1.0])



## Add indication of max delta alpha
plt.sca(Axes['11'])
chis_alpha_diff_min = chi_list[np.argmin(alphas_diff/taper_angles,axis=1)]
chis_alpha_diff_max = chi_list[np.argmax(alphas_diff/taper_angles,axis=1)]
#Lambdas_alpha_diff_0 = LambdaRef_list[np.argmax(np.abs(alphas_diff[0:-1,:]),axis=0)]
#
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_max*100.0,'--k')
#plt.plot(Lambdas_alpha_diff_0*100.0,chi_list*100.0)
#plt.cla()
#plt.contour(alphas_diff)

#plt.sca(Axes['12'])
bDalpha = alphas_diff/taper_angles
dum = (bDalpha[1:,:]-bDalpha[0:-1,:])/((chi_list[1]-chi_list[0])*100.0)
dAlpha_dChi = (dum[:,1:] + dum[:,0:-1])/2.0
dum = (bDalpha[:,1:]-bDalpha[:,0:-1])/((LambdaRef_list[1]-LambdaRef_list[0])*100.0)
dAlpha_dLambda = (dum[1:,:] + dum[0:-1,:])/2.0

#plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, dAlpha_dChi,vmin=-.1,vmax=.1)
#plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, dAlpha_dLambda,1000,vmin=-.1,vmax=.1)
#CS = plt.contourf(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, -np.sqrt(dAlpha_dLambda**2+dAlpha_dChi**2),1000)
vGrad = np.sqrt(dAlpha_dLambda**2+dAlpha_dChi**2)



#CS = plt.contour(Lambdas*100.0, chis*100.0, alphas_diff/taper_angles,np.linspace(.00,0,2),vmin=-1.00,vmax=1.00)
#plt.colorbar()
r = 5
#plt.set_cmap('plasma')
dum = (Lambdas[0:-1,:] + Lambdas[1:,:])/2.0
Lambdas_centered = (dum[:,0:-1] + dum[:,1:])/2.0

dum = (chis[0:-1,:] + chis[1:,:])/2.0
chis_centered = (dum[:,0:-1] + dum[:,1:])/2.0

#CS = plt.contour(Lambdas[0:-1,0:-1]*100.0, chis[0:-1,0:-1]*100.0, vGrad,100)



#plt.quiver(Lambdas[0:-1:r,0:-1:r]*100.0, chis[0:-1:r,0:-1:r]*100.0, dAlpha_dLambda[::r,::r],dAlpha_dChi[::r,::r],scale=.5)
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_max*100.0,'--w')
#plt.plot(Lambdas_alpha_diff_0[:-1]*100.0,chi_list[:-1]*100.0,'--w')


chi_list_centered = (chi_list[0:-1] + chi_list[1:])/2.0
LambdaRef_list_centered = (LambdaRef_list[0:-1] + LambdaRef_list[1:])/2.0
chis_vGrad_min = chi_list_centered [np.argmin(vGrad,axis=1)]



#plt.sca(Axes['11'])
#plt.plot(LambdaRef_list_centered*100.0,chis_vGrad_min*100.0,'--k')


dy = ((chi_list[1]-chi_list[0])*100.0)
dx = ((LambdaRef_list[1]-LambdaRef_list[0])*100.0)
Z = vGrad
K = 0.000
end = nLambda-1
I = slice(1,end-1)
Ip = slice(2,end)
Im = slice(0,end-2)
for it in range(100):
    Z[I,I] += K * (  (Z[Ip,I]-2.0*Z[I,I]+Z[Im,I])/dy/dy  +  (Z[I,Ip]-2.0*Z[I,I]+Z[I,Im])/dx/dx  )


#for i in range(nLambda-1):
#    Z[i,:] -= np.min(Z[i,:])
#    Z[i,:] /= np.max(Z[i,:])

chis_vGrad_min = chi_list_centered [np.argmin(vGrad,axis=1)]
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
#chis_vGrad_min[3:] = (chis_vGrad_min[:-3] + chis_vGrad_min[1:-2] + chis_vGrad_min[1:])/2.0

plt.sca(Axes['12'])
plt.cla()
CS = plt.contour(Lambdas_centered*100.0, chis_centered*100.0, vGrad,[1.1e-2, 1e6])
plt.sca(Axes['11'])
plt.plot(LambdaRef_list_centered*100.0,chis_vGrad_min*100.0,'--k',linewidth=0.5)
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')
#
#Lambdas_alpha_diff_max = LambdaRef_list[np.argmax(alphas_diff,axis=0)]
#Lambdas_alpha_diff_min = LambdaRef_list[np.argmin(alphas_diff,axis=0)]
#plt.plot(Lambdas_alpha_diff_max*100.0,chi_list*100.0,'--k')
#plt.plot(Lambdas_alpha_diff_min*100.0,chi_list*100.0,'--k')
#plt.plot(LambdaRef_list*100.0,chis_alpha_diff_min*100.0,'--k')



## Find types
Type = np.zeros([nLambda-1,nChi-1])
for iL in range(nLambda-1):
    chi_boundary = chis_vGrad_min[iL]
    for iC in range(nChi-1):
        chi = chi_list_centered[iC]
        if chi<chi_boundary:
            Type[iC][iL] = 1
        else:
            if alphas_diff[iL][iC]<0.0:
                Type[iC][iL]=3
            else:
                Type[iC][iL]=2
    #end iC
#end iL
    
## Find boundaries
#Type = np.zeros([nLambda-1,nChi-1])
bound_10 = []
bound_12 = []
bound_21 = []
bound_23 = []
bound_32 = []
bound_30 = []
for iC in range(nChi-1):
    y = chi_list_centered[iC]
    for iL in range(nLambda-1):
        x = LambdaRef_list_centered[iL]
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
        elif Type[iC,iL]==3 and (iC==nChi-2):
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
        
        
#plt.plot(zeroContour[:,0],zeroContour[:,1],'-k')
#

#plt.plot(LambdaRef_list_centered*100.0,chis_vGrad_min*100.0,'-g',linewidth=0.5)


bound_12 = arr([LambdaRef_list_centered,chis_vGrad_min]).T
bound_12 = np.concatenate([arr([[0.0,chis_vGrad_min[0]]]),bound_12])
bound_23 = np.flipud(zeroContour)/100.0
bound_23 = np.concatenate([arr([[1.0,0.0]]),bound_23,arr([[0.0,1.0]])])


domain1 = np.concatenate([arr([[0.0,0.0]]),bound_12])
domain2 = np.concatenate([bound_12,bound_23])
domain3 = np.concatenate([bound_23,arr([[1.0,1.0]])])


#plt.fill(domain1[:,0]*100.0,domain1[:,1]*100.0,color=arr([219, 59, 38])/255.0)
#plt.fill(domain2[:,0]*100.0,domain2[:,1]*100.0,color=arr([239,189, 64])/255.0)
#plt.fill(domain3[:,0]*100.0,domain3[:,1]*100.0,color=arr([ 80,159,248])/255.0)
#plt.plot(bound_23[:,0]*100.0,bound_23[:,1]*100.0,'-k')

CMAP = [[.99,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,1.0]]

colors=[(0.75, 0.15, 0.15), (1, 0.75, 0.15), (0.15, 0.75, 0.15)]
#plt.set_cmap("Set2")

                


ax = plt.gca()
#ax.tick_params(axis='x',top=True,bottom=False,labeltop=True,labelbottom=False)
ax.xaxis.tick_top()
#ax.invert_yaxis()
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])

plt.text(10,25,'I',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(20,70,'II',family='Times New Roman',color='w',size=28,weight='bold')
plt.text(60,90,'III',family='Times New Roman',color='w',size=28,weight='bold')


nBP = bound_12.shape[0] # number of boundary points
#for iBP in nBP:
Dist_12 = np.zeros([nChi-1,nLambda-1])
Dist_23 = np.zeros([nChi-1,nLambda-1])
for iC in range(nChi-1):
    y = chi_list_centered[iC]
    for iL in range(nLambda-1):
        x = LambdaRef_list_centered[iL]
        if Type[iC,iL] == 2:
            Dist_12[iC,iL] = np.min( (x-bound_12[:,0])**2 + (y-bound_12[:,1])**2)
            Dist_23[iC,iL] = np.min( (x-bound_23[:,0])**2 + (y-bound_23[:,1])**2)

plt.sca(Axes['12'])

#plt.contourf(LambdaRef_list_centered*100.0,chi_list_centered*100.0,Dist_12,np.linspace(0.01,1.0,10))









plt.sca(Axes['12'])
plt.cla()
# Get distance to boundary
Dist_10 = np.zeros([nChi-1,nLambda-1]) 
Dist_12 = np.zeros([nChi-1,nLambda-1]) 
Dist_21 = np.zeros([nChi-1,nLambda-1]) 
Dist_23 = np.zeros([nChi-1,nLambda-1]) 
Dist_32 = np.zeros([nChi-1,nLambda-1]) 
Dist_30 = np.zeros([nChi-1,nLambda-1]) 
floatType = np.zeros([nChi-1,nLambda-1]) 
for iC in range(nChi-1):
    y = chi_list_centered[iC]
    for iL in range(nLambda-1):
        x = LambdaRef_list_centered[iL]
        if Type[iC,iL] == 1:
            Dist_low = np.min( (x-bound_10[:,0])**2 + (y-bound_10[:,1])**2)
            Dist_up = np.min( (x-bound_12[:,0])**2 + (y-bound_12[:,1])**2)
            
        elif Type[iC,iL] == 2:
            Dist_low = np.min( (x-bound_21[:,0])**2 + (y-bound_21[:,1])**2)
            Dist_up = np.min( (x-bound_23[:,0])**2 + (y-bound_23[:,1])**2)
            
        elif Type[iC,iL] == 3:
            Dist_low = np.min( (x-bound_32[:,0])**2 + (y-bound_32[:,1])**2)
            Dist_up = np.min( (x-bound_30[:,0])**2 + (y-bound_30[:,1])**2)
            
            
        if (Dist_low+Dist_up)!=0.0:
            floatType[iC,iL] = Type[iC,iL] + Dist_low/(Dist_low+Dist_up)
        else:
            floatType[iC,iL] = Type[iC,iL]
#plt.contourf(LambdaRef_list_centered*100.0,chi_list_centered*100.0,Dist_12,np.linspace(0.01,0.1,10))
plt.contourf(LambdaRef_list_centered*100.0,chi_list_centered*100.0,floatType,np.linspace(1.000,4.0,700),vmin=1.0,vmax=4.0)
#plt.plot(chi_list_centered*100.0,floatType[:,20])
#plt.plot(bound_23[:,0]*100.0,bound_23[:,1]*100.0,'-k')
##plt.fill(domain1[:,0]*100.0,domain1[:,1]*100.0,color=arr([219, 59, 38])/255.0)
#plt.plot(domain2[:,0]*100.0,domain2[:,1]*100.0,'--k')
#plt.fill(domain2[:,0]*100.0,domain2[:,1]*100.0,color=arr([239,189, 64])/255.0)
##plt.fill(domain3[:,0]*100.0,domain3[:,1]*100.0,color=arr([ 80,159,248])/255.0)
cbar = plt.colorbar()
cbar.set_ticks([1.0,2.0,3.0,4.0])
ax.xaxis.set_label_position('top')
plt.axis([.0,100.0,100.0,.0])
#plt.set_cmap("Set2")
#colors = arr([[219, 59, 38],
#              [239,189, 64],
#              [ 80,159,248],
#              [ 20, 50,150]]) / 255.0


w = 0.000 # width of the border
segPos = [0.0, 1.0/3.0-w, 1.0/3.0+w, 2.0/3.0-w, 2.0/3.0+w, 1.0]

colors = arr([[219, 59, 38],
              [239,189, 64],
              [255,255,255],
              [ 80,159,248],
              [ 20, 50,150]]) / 255.0

keys = ['red','green','blue']
segmentdata = {}
# flat version
for i in range (3):
    segmentdata[keys[i]] = [ (segPos[0],  colors[0][i], colors[0][i]),               
                             
                             (segPos[1],  colors[0][i], colors[1][i]),
                             (segPos[2],  colors[1][i], colors[2][i]),
                               
                             (segPos[3],  colors[2][i], colors[3][i]),
                             (segPos[4],  colors[3][i], colors[4][i]),
                               
                             (segPos[5],  colors[4][i], colors[4][i]) ]
    
# mode 2 smooth version   
for i in range (3):
    segmentdata[keys[i]] = [ (segPos[0],  colors[0][i], colors[0][i]),               
                             
                             (segPos[1],  colors[0][i], colors[1][i]),
                             (segPos[2],  colors[1][i], colors[1][i]),
                               
                             (segPos[3],  colors[3][i], colors[3][i]),
                             (segPos[4],  colors[3][i], colors[4][i]),
                               
                             (segPos[5],  colors[4][i], colors[4][i]) ]
    

# everything smoother
colors = arr([[200, 30, 32],
              [248,120, 64],
              [248,180, 64],
              [ 64,180,248],
              [ 64,120,248],
              [ 32, 30,200]]) / 255.0
   
for i in range (3):
    segmentdata[keys[i]] = [ (segPos[0],  colors[0][i], colors[0][i]),               
                             
                             (segPos[1],  colors[1][i], colors[1][i]),
                             (segPos[2],  colors[2][i], colors[2][i]),
                               
                             (segPos[3],  colors[3][i], colors[3][i]),
                             (segPos[4],  colors[4][i], colors[4][i]),
                               
                             (segPos[5],  colors[5][i], colors[5][i]) ]
    
#CMAP = LinearSegmentedColormap.from_list('custom',colors,N=9)
CMAP = LinearSegmentedColormap('custom', segmentdata)

#plt.pcolor(LambdaRef_list_centered*100.0,chi_list_centered*100.0,Type)
plt.register_cmap(cmap=CMAP)
plt.set_cmap("custom")