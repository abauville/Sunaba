#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:22:22 2018

@author: abauville
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:58 2018

@author: abauville
"""

import sys
import os
sys.path.insert(0, '../../../src/UserInput')
sys.path.insert(0, '../CriticalTaper')
sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
#import CritTaper_dataMaker
import Figz_Utils
#import CritTaper_Style
from numpy import array as arr
import OutputDef as Output
from PaperDecollement_Utils import getColormap, get_XYandPattern
from matplotlib.colors import LinearSegmentedColormap



#chi_list = [ 1, 10, 20, 30, 40, 50, 60, 70, 80]
chi_list = [40]
#chi_list = [15,20,25,30,35]
Lambda_list = arr([60])
tstep_list = np.arange(000,400,50)
tstep_list = np.concatenate([tstep_list,np.arange(400,550,5)])
tstep_list = np.concatenate([tstep_list,np.arange(550,751,50)])

#tstep_list = np.arange(000,801,50)


nC = len(chi_list)
nL = len(Lambda_list)
nStep = len(tstep_list)
nHor = nL
nVer = nStep


aspectRatio = 0.2
fig  = Figz_Utils.Figure(106,height=29.7,width=21.0,mode='draft')

#bigAxes = Figz_Utils.makeAxes(fig,1,1,aspectRatio=0.66,leftMarginPad=1.25,rightMarginPad=0.25,topMarginPad=1.5,bottomMarginPad = 0.0,xPad = 0.5,yPad=.25,setAspectRatioBasedOn='x')
#ax = plt.gca()
##plt.xlabel("$\\mathbf{\\lambda}$ [%]",weight='bold',verticalAlignment='center')
#ax.xaxis.tick_top()
#ax.xaxis.set_label_position('top')
#ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
##plt.axis([.0,100.0,100.0,.0])

Axes = Figz_Utils.makeAxes(fig,nVer,nHor,aspectRatio=aspectRatio,leftMarginPad=1.5,rightMarginPad=0.25,topMarginPad=.0,bottomMarginPad = 0.0,xPad = 0.5,yPad=.00,setAspectRatioBasedOn='x')


fig2  = Figz_Utils.Figure(107,height=21.0,width=21.0,mode='draft')
AxesPressure = Figz_Utils.makeAxes(fig2,1,1,bottomMarginPad=4.5)
AxesPsi = Figz_Utils.makeAxes(fig2,1,1,aspectRatio=.15,topMarginPad=15.0)
apectRatioMohr = 0.5
AxesMohr = Figz_Utils.makeAxes(fig2,4,1,aspectRatio=apectRatioMohr,topMarginPad=0.5,leftMarginPad=1.5,rightMarginPad=13.0)
for i in range(4):
    plt.sca(AxesMohr['%i1' % (i+1)])
    ax = plt.gca()
    plt.xticks([])
    plt.yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


superRootFolder = "/Users/abauville/Output/Paper_Decollement/Output_AllSteps/wWater/Beta00/"
superDirList = []
i = 0
for iC in range(nC):
    for iL in range(nL):
        superDirList.append("Weak%02d/Lambda%02d" % (chi_list[iC], Lambda_list[iL]))



ProductionMode = False
if ProductionMode:
#    sampleRate = 1
#    pointSize = 0.01
    sampleRate = 1
    pointSize = sampleRate/3.0
else:
    sampleRate = 300
    pointSize = sampleRate/20.0


nSim = nC*nL
iSim = 0


### Choose a point from final step
#tstep = 400 
#outFolder = "Out_%05d" % (tstep)
#dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
#xP = -1.8
#yP = 0.7
#lc = 2.0e3
#PartX       = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::sampleRate]/lc
#PartY       = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::sampleRate]/lc
#PartXIni    = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::sampleRate]/lc
#PartYIni    = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::sampleRate]/lc
#I = np.argmin(np.abs(PartX-xP)**2 + np.abs(PartY-yP)**2)
#xPIni = PartXIni[I]
#yPIni = PartYIni[I]
#plt.plot(xP,yP,'or')


# Choose a point arbitrarily
n = 1
#xPIni = np.linspace(-8.5,-7.5,n)
xPIni = np.linspace(-8.5,-5.5,n)
#xPIni = np.linspace(-12.25,-12.,n)
yPIni = arr([.75])

xPIni,yPIni = np.meshgrid(xPIni,yPIni)
xPIni = xPIni.flatten()
yPIni = yPIni.flatten()


nP = len(xPIni)

Pressure = np.zeros([nStep,nP])
PWaterColumn = np.zeros([nStep,nP])
Pl = np.zeros([nStep,nP])
Pf = np.zeros([nStep,nP])
Sxx = np.zeros([nStep,nP])
Sxy = np.zeros([nStep,nP])
SII = np.zeros([nStep,nP])
Ty = np.zeros([nStep,nP])
S1 = np.zeros([nStep,nP])
S3 = np.zeros([nStep,nP])
psi = np.zeros([nStep,nP])
depth = np.zeros([nStep,nP])
depth_w = np.zeros([nStep,nP])
xP  = np.zeros([nStep,nP])
yP = np.zeros([nStep,nP])
tstep_list_mat = np.zeros([nStep,nP])
#lc=2.0e3
rho = 2500.0
rho_w = 1000.0
g = 9.81
Lambda = Lambda_list[0]/100.0
for iStep in range(nStep):
    tstep = tstep_list[iStep]
    print("tstep = %i/%i" % (tstep,tstep_list[-1]))
    outFolder = "Out_%05d" % (tstep)
    dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
    
    
    if iStep == 0:
        Setup = Output.readInput(superRootFolder + superDirList[iSim] + '/Input/input.json')
        cohesion = Setup.MatProps['1'].cohesion
        frictionAngle = Setup.MatProps['1'].frictionAngle
        rho = Setup.MatProps['1'].rho0
        rho_w = Setup.MatProps['0'].rho0
        lc = Setup.Char.length
        raw = Output.getData(dataFolder + 'P.bin',True)
        
        xmin = raw.xmin; ymin = raw.ymin
        xmax = raw.xmax; ymax = raw.ymax
        x = np.linspace(raw.xmin,raw.xmax,raw.nx)
        y = np.linspace(raw.ymin,raw.ymax,raw.ny)
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        X, Y = np.meshgrid(x,y)
        X = X.T
        Y = Y.T
    
    
    PartXIni    = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data[0::1]/lc
    PartYIni    = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data[0::1]/lc
    PartX       = Output.getParticleData(dataFolder + 'particles_x.bin',True).data[0::1]/lc
    PartY       = Output.getParticleData(dataFolder + 'particles_y.bin',True).data[0::1]/lc
    phase       = Output.getData(dataFolder + 'phase.bin',False).data
    
    for iP in range(nP):
        tstep_list_mat[:,iP] = tstep_list
        I = np.argmin(np.abs(PartXIni-xPIni[iP])**2 + np.abs(PartYIni-yPIni[iP])**2)
        xP[iStep,iP] = PartX[I]
        yP[iStep,iP] = PartY[I]
        ix = np.argmin(np.abs(xP[iStep,iP]-x))
        iy = np.argmin(np.abs(yP[iStep,iP]-y))
        depth[iStep,iP] = np.sum(phase[ix,-1:iy-1:-1])*dy*lc
        depth_w[iStep,iP] =  np.sum(np.abs(phase[ix,-1:iy:-1]-1.0))*dy*lc
        
        
        Pl[iStep,iP] = depth[iStep,iP]*rho*g# + depth_w[iStep,iP]*rho_w*g*(1.0-Lambda)

        
#        print("ix=%i, iy=%i" % (ix,iy))
        Pressure[iStep,iP] = Output.getData(dataFolder + 'P.bin',True).data[ix,iy]
        Isurf = np.argmin(phase[ix,:])
        PWaterColumn[iStep,iP] = Output.getData(dataFolder + 'P.bin',True).data[ix,Isurf]
        
        Pressure[iStep,iP] = Pressure[iStep,iP]-PWaterColumn[iStep,iP]
        
        Sxx[iStep,iP] = Output.getData(dataFolder + 'sigma_xx.bin',True).data[ix,iy]
        Sxy[iStep,iP] = Output.getData(dataFolder + 'sigma_xy.bin',True).data[ix,iy]
        SII[iStep,iP] = np.sqrt(Sxx[iStep,iP]**2+Sxy[iStep,iP]**2)
        S1[iStep,iP] = Pressure[iStep,iP]*(1.0-Lambda)+SII[iStep,iP]
        S3[iStep,iP] = Pressure[iStep,iP]*(1.0-Lambda)-SII[iStep,iP]
        
        Ty[iStep,iP] = cohesion*np.cos(frictionAngle) + (1.0-Lambda)*Pressure[iStep,iP]*np.sin(frictionAngle)
        
        Tau = Sxx[iStep,iP] / Sxy[iStep,iP];

        thisPsi = np.zeros(Sxy[iStep,iP].shape)
        if Sxy[iStep,iP]<0.0:
            psi[iStep,iP] = np.arctan(-Tau+np.sqrt(Tau**2+1.0))
        else:
            psi[iStep,iP] = np.arctan(-Tau-np.sqrt(Tau**2+1.0))
        

#        Pf[iStep,iP] = (depth_w*rho_w + depth[iStep,iP]*rho*Lambda_list[0])*g
        Pf[iStep,iP] = depth_w[iStep,iP]*rho_w*g + Pressure[iStep,iP]*Lambda

plt.sca(AxesPressure['11'])
#plt.plot(tstep_list_mat,Pressure,'-')
#plt.plot(tstep_list_mat,Pressure,'-')
#Pl = depth*lc*2500.0*9.81

P_ocean = depth_w*rho_w*g
Pl=Pl
#    plt.plot(xP,yP,'or')

# find the max
ISteps = np.argmax(Pressure,axis=0)
IP = np.argmax(np.max(Pressure,axis=0))
IStep = ISteps[IP]

## P vs Pl
#plt.plot(arr([0.0,Pressure[IStep,IP]/1e6]),arr([0.0,Pressure[IStep,IP]/1e6]),'--')
#plt.plot(Pl/1e6,(S3)/1e6,'-ob',linewidth=.5,markersize=.75)
#plt.plot(Pl/1e6,(Pressure)/1e6,'-ok',linewidth=.5,markersize=.75)
#plt.plot(Pl/1e6,(Pressure*(1.0-Lambda))/1e6,'--ok',linewidth=.5,markersize=.75)
#plt.plot(Pl/1e6,(S1)/1e6,'-or',linewidth=.5,markersize=.75)

#P vs time
#plt.plot(tstep_list,(Pressure)/1e6,'-ok',linewidth=.5,markersize=.75)
plt.plot(tstep_list,(Pressure)*(1.0-Lambda)/1e6,'-ok',linewidth=.5,markersize=.75)
plt.fill(np.concatenate([tstep_list,np.flipud(tstep_list)]),np.concatenate([ ((Pressure)*(1.0-Lambda)+Ty)/1e6,np.flipud((Pressure)*(1.0-Lambda)-Ty)/1e6 ]),'k',alpha=.05)
plt.plot(tstep_list,(S1)/1e6,'-or',linewidth=.5,markersize=.75)
plt.plot(tstep_list,(S3)/1e6,'-ob',linewidth=.5,markersize=.75)
plt.plot(tstep_list,(Pl)*(1.0-Lambda)/1e6,'--b',linewidth=.5,markersize=.75)

plt.xlim(0,tstep_list[-1])

#tstep_list_plot = arr([100,200,300,400,500,600,700])
tstep_list_plot = arr([445,480,530,550])
nStepPlot = len(tstep_list_plot)
#tstep_list_plot = [400,425,450,475,500]
#for iStep in range(nStepPlot):
#    IStep = np.argmin(np.abs(tstep_list-np.abs(tstep)))
#    plt.text(tstep_list[IStep],1.0+S1[IStep]/1e6,tstep_list[IStep],verticalAlignment='baseline')

# Plot stress orientation
# ==========================================
plt.sca(AxesPsi['11'])
#plt.plot(arr([0.0,Pressure[IStep,IP]])/1e6,[90.0,90.0],':k')
#plt.plot(arr([0.0,Pressure[IStep,IP]])/1e6,[45,45],':y')
#plt.plot(arr([0.0,Pressure[IStep,IP]])/1e6,[0.0,0.0],':k')

plt.plot(arr([0.0,tstep_list[-1]]),[90.0,90.0],':k')
plt.plot(arr([0.0,tstep_list[-1]]),[45,45],':y')
plt.plot(arr([0.0,tstep_list[-1]]),[0.0,0.0],':k')

#plt.plot(Pl/1e6,(psi)*180.0/np.pi,'-om',linewidth=.5,markersize=.75)
plt.plot(tstep_list,(psi)*180.0/np.pi,'-om',linewidth=.5,markersize=.75)
plt.ylabel('$\\psi$ [°]')
plt.xlabel('timesteps')
plt.xlim(0,tstep_list[-1])



# Plot Mohr Diagram
# ==========================================
stepMohr = [445,480,530,550]

for iMohr in range (len(stepMohr)):
    I = np.argmin(np.abs(tstep_list-stepMohr[iMohr]))
    plt.sca(AxesPressure['11'])
    plt.plot(stepMohr[iMohr],Pressure[I,0]*(1.0-Lambda)/1e6,'or',markerFaceColor='None')
    plt.text(stepMohr[iMohr],1.0+Pressure[I,0]*(1.0-Lambda)/1e6,iMohr)
    
    plt.sca(AxesPsi['11'])
    plt.plot(stepMohr[iMohr],psi[I,0]*180.0/np.pi,'or',markerFaceColor='None')
    plt.text(stepMohr[iMohr],-30.0+psi[I,0]*180.0/np.pi,iMohr)
    
    plt.sca(AxesMohr['%i1' % (iMohr+1)])
    phi = np.linspace(0,np.pi,100)
    
    maxSigmaN = (Pressure[I,0]*(1.0-Lambda)+SII[I,0])*1.2
#    plt.plot([0.0,maxSigmaN],[0.0,cohesion*np.cos(frictionAngle)+maxSigmaN*np.sin(frictionAngle)],'-k')
    plt.plot([0.0,maxSigmaN],[0.0,cohesion+maxSigmaN*np.tan(frictionAngle)],'-k')
    
    plt.plot(Pressure[I,0]*(1.0-Lambda)+SII[I,0]*np.cos(phi),SII[I,0]*np.sin(phi),'-k')
    plt.plot(Pressure[I,0]*(1.0-Lambda),0.0,'|k',markerSize=5.0)
    plt.plot(Pl[I,0]*(1.0-Lambda),0.0,'|r',markerSize=5.0)
    plt.xlim(0.0,maxSigmaN)
    plt.ylim(0.0,maxSigmaN*apectRatioMohr)
    


#plt.plot(Pl/1e6,(S1+P_ocean)/1e6,'-',linewidth=.5)
#plt.plot(Pl/1e6,(Pressure+P_ocean)/1e6,'-x',linewidth=.5)
#plt.plot(Pl/1e6,(S3+P_ocean)/1e6,'-',linewidth=.5)
#
#plt.plot((S1+P_ocean)/1e6,'-',linewidth=.5)
#plt.plot((Pressure+P_ocean)/1e6,'-x',linewidth=.5)
#plt.plot((S3+P_ocean)/1e6,'-',linewidth=.5)

#plt.plot(Pl[IStep,IP]/1e6,Pressure[IStep,IP]/1e6,'ok')


#tstep_list_plot = [0]


Plot=True
#Plot=False
if Plot==True:
    iPlot = 0
    for iC in range(nC):
        for iL in range(nL):
            for iStep in range(nStepPlot):
                tstep = tstep_list_plot[iStep]
                
                IStep = np.argmin(np.abs(tstep_list-np.abs(tstep)))
                
                print("tstep = %i/%i" % (tstep,tstep_list[-1]))
                outFolder = "Out_%05d" % (tstep)
                dataFolder = superRootFolder + superDirList[iSim] + "/Output/" + outFolder + "/"
                Char = Output.readInput(superRootFolder + superDirList[iSim] + '/Input/input.json').Char
                timeSim = Output.readState(dataFolder + "modelState.json").time*Char.time
                
                PartX = []
                PartY = []
                PartPattern = []
                
    #            ax = plt.sca(Axes["%i%i" % (iC+1,iL+1)])
                ax = plt.sca(Axes["%i%i" % (iPlot+1,iL+1)])
                
        
                ymax = 3.5
                plt.axis([-1.0/aspectRatio*ymax,0.0,0.0,ymax])
                plt.axis("off")
                rx = 1
                ry = 1
                
                phase = Output.getData(dataFolder + 'phase.bin',False).data[::rx,::ry]
                mask = phase==0
                
                
                
                # Stress massaging
                # =============================================
                Sxx = Output.getData(dataFolder + 'sigma_xx.bin',True).data[::rx,::ry]
                Sxy = Output.getData(dataFolder + 'sigma_xy.bin',True).data[::rx,::ry]
                
                if iStep == 0:
                    raw = Output.getData(dataFolder + 'sigma_xx.bin',True)
                    xmin = raw.xmin; ymin = raw.ymin
                    xmax = raw.xmax; ymax = raw.ymax
                    x = np.linspace(raw.xmin,raw.xmax,raw.nx)
                    y = np.linspace(raw.ymin,raw.ymax,raw.ny)
                    dx = x[1]-x[0]
                    X, Y = np.meshgrid(x,y)
                    X = X.T
                    Y = Y.T
                Tau = Sxx / Sxy;
                SII = np.sqrt(Sxx**2+Sxy**2);
        
                psi = np.zeros(Sxy.shape)
                I = Sxy<0.0
        
                psi[I] = np.arctan(-Tau[I]+np.sqrt(Tau[I]**2+1.0))
                psi[~I] = np.arctan(-Tau[~I]-np.sqrt(Tau[~I]**2+1.0))
                
        #        psi = np.ma.masked_array(psi, mask_sub)
        
                length = 0.2 * phase[::rx,::ry].flatten()
        
                xLine = arr([-np.cos(psi.flatten())*length,np.cos(psi.flatten())*length])
                yLine = arr([-np.sin(psi.flatten())*length,np.sin(psi.flatten())*length])
                
                Svec_x = np.cos(psi)
                Svec_y = np.sin(psi)
                
                Svec_x *= phase
                Svec_y *= phase
                
                
                # Plot particles
                # =============================================
                
                PartX, PartY, PartPattern, nColors = get_XYandPattern(dataFolder, sampleRate=sampleRate, nLayersX=0, nLayersY=0.00,maxStrain=5.0)
                if iStep == 0:
                    CMAP = arr([ [0.4,0.5,0.8,1.0] ])
                    CMAP = getColormap(nColors,"myColorMap",CMAP=CMAP,shiftHLayerColors=False,strainDarknessFactor=0.0)
                    plt.register_cmap(cmap=CMAP)
                    plt.set_cmap("myColorMap")
                plt.scatter(PartX,PartY,c=PartPattern,s=pointSize,vmin=0.0,vmax=4*nColors-1,edgecolors='None')
                
                
                plt.plot(xP[IStep,: ],yP[IStep,: ],'or',markerFaceColor='None')
                plt.plot(xP[IStep,IP],yP[IStep,IP],'or',markerFaceColor='y')
               
                
                
                # Plot stresses
                # =============================================
                smin = 0.0*1e6
                smax = 10.0*1e6
                SII[SII>smax] = smax
                SII[SII<smin] = smin
                n = 100
                start_points = arr([np.linspace(xmin+dx,xmax-dx,n),
                                    np.linspace(dx,dx,n)]).T
                
                colors = arr([[  0,255,255],
                              [255,255,255],
                              [255,255,  0]]) / 255.0
                nTot = 256
                CMAP = LinearSegmentedColormap.from_list('custom2',colors,N=nTot)       
                plt.register_cmap(cmap=CMAP)
                plt.streamplot(x[::rx],y[::ry],Svec_x.T,Svec_y.T,color=SII.T, density=2.0,arrowsize=0.01,cmap='custom2')        
        
#                plt.text((xmin+xmax)/2.0,ymax-0.5,"tstep = %i", tstep,horizontalAlignment='center')
    
    
#                P = Output.getData(dataFolder + 'P.bin',True).data
#                plt.contourf(x,y,P.T/1e6,256,cmap='jet')
#                plt.colorbar()
                
                iPlot+=1
            iSim+=1
