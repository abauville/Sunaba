# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys
sys.path.insert(0, '/Users/abauville/Work/StokesFD/src/UserInput')

import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.pyplot as plt

import OutputDef as Output
import InputDef as Input
s       = 1.0
mn      = 60        * s
hour    = 60        * mn
day     = 24        * hour
yr      = 365       * day

# Create the colormap c-b-k-r-y
#condensedCMAP = np.array([[0.0,1.0,1.0,1.0], # cyan
#                          [0.0,0.0,1.0,1.0], # blue
#                          [0.0,0.0,0.0,1.0], # black
#                          [1.0,0.0,0.0,1.0], # red
#                          [1.0,1.0,0.0,1.0]]) #yellow

#condensedCMAP = np.array([
#                          [0.0,0.0,1.0,1.0], # blue
#                          [0.0,0.0,1.0,1.0], # blue
#                          [0.0,0.0,1.0,1.0], # blue
#                          [0.0,0.0,0.0,1.0], # black
#                          [1.0,0.0,0.0,1.0], # red
#                          [1.0,0.0,0.0,1.0], # red
#                          [1.0,0.0,0.0,1.0], # red
#                                           ]) #yellow

### Create the colormap passive
#shadeFac = 0.25
#shadeArray = np.array([shadeFac,shadeFac,shadeFac,1.0])
#condensedCMAP = np.array([np.array([0.0,0.8,0.5,1.0]),   
#                          np.array([0.0,0.8,0.5,1.0])*shadeArray, 
#                          np.array([1.0,0.9,0.6,1.0]),
#                          np.array([1.0,0.9,0.6,1.0])*shadeArray,]) 
#
#shadeFac = 1.0
#shadeArray = np.array([1.0,shadeFac,shadeFac,1.0])
#condensedCMAP = np.array([np.array([0.25,0.25,0.25,1.0]),   
#                          np.array([.25,0.25,0.25,1.0])*shadeArray, 
#                          np.array([0.8,0.8,0.8,1.0]),
#                          np.array([0.8,0.8,0.8,1.0])*shadeArray,]) 


# Create the colormap many random colors
nRandColors  = 24
condensedCMAPtemp = np.random.rand(nRandColors,4) * 0.5+.1
condensedCMAPtemp[:,0] += .4
condensedCMAPtemp[:,1] += .25
#condensedCMAPtemp[:,3] += .1
condensedCMAP = np.zeros((2*nRandColors,4))
condensedCMAP[0::2,:] = condensedCMAPtemp
condensedCMAP[1::2,:] = condensedCMAPtemp*0.25
#condensedCMAP = np.concatenate(condensedCMAP,0.5*condensedCMAP)
condensedCMAP[:,-1] = 1.0


CMAP = np.zeros((255,4))

n = condensedCMAP.shape[0]
for i in range(n-1):
    iS = int(round(255 *  i/(n-1) ))
    iE = int(round(255 * (i+1)/(n-1)))
    for j in range(4):
        CMAP[iS:iE,j] = np.linspace(condensedCMAP[i,j],condensedCMAP[i+1,j],int(iE-iS))

CMAP = (CMAP*255.0).astype("uint8")

#n = 100000
#x = np.zeros(n)
#y = np.random.rand(n)
#z = np.random.rand(n)
#s = np.random.rand(n)
iStep = 959
rootFolder = "/Users/abauville/Work/Output/Test/Output/"
outFolder = "Out_%05d" % iStep #DirList[iStep]
Setup = Output.readInput(rootFolder +  'Input/input.json')
#Char = Setup.Char
#CharExtra = Input.CharExtra(Char)

dataFolder = rootFolder + outFolder + "/"
State       = Output.readState(dataFolder + "modelState.json")
timeSim   = (State.time+ State.dt) * Setup.Char.time

#            phase = Output.getData(dataFolder + 'phase.bin').data
#            mask = phase == 0

#            strainRate  = Output.getData(dataFolder + 'strainRate.bin',True,mask).data

lc = 2.0 * 1000.0

print("load data")
PartX       = Output.getParticleData(dataFolder + 'particles_x.bin',True).data/lc
PartY       = Output.getParticleData(dataFolder + 'particles_y.bin',True).data/lc
PartXIni    = Output.getParticleData(dataFolder + 'particles_xIni.bin',True).data/lc
PartYIni    = Output.getParticleData(dataFolder + 'particles_yIni.bin',True).data/lc
PartStrain  = Output.getParticleData(dataFolder + 'particles_strain.bin',True).data
#PartTimeLastPlastic  = Output.getParticleData(dataFolder + 'particles_timeLastPlastic.bin',True).data / yr


# subsampling
sr = 1 # sample rate
PartX    = PartX[0::sr]
PartY    = PartY[0::sr]
PartXIni = PartXIni[0::sr]
PartYIni = PartYIni[0::sr]
PartStrain    = PartStrain[0::sr]

xmax = np.max(PartXIni)
xmin = np.min(PartXIni)
PartXPat = PartXIni
PartXPat -= xmin
PartXPat /= xmax-xmin

nLayers = 8
YPattern = np.cos(PartYIni*1.0*np.pi*nLayers)
nLayers = 1
#XPattern = np.sin(PartXIni*1.0*np.pi*nLayers)
XPattern = np.sin(PartXPat*0.5*np.pi*nLayers)


PartPattern = YPattern*XPattern
maxVal= np.max(PartPattern)
#PartPattern = np.round((PartPattern+maxVal)/(2.0*maxVal))

maxVal= np.max(XPattern)
#XPattern = np.round((XPattern+maxVal)/(2.0*maxVal))

maxVal= np.max(YPattern)
YPattern = np.round((YPattern+maxVal)/(2.0*maxVal))

#PartPattern+=2.0*YPattern
#PartPattern=2.0*YPattern
#PartPattern=2.0*XPattern

PartPattern= 2.0*np.round(XPattern * (nRandColors-1))
#PartPattern+= YPattern

#plt.scatter(PartX,PartY,c=PartPattern)



x = np.zeros(PartX.size)
print("render")
scene = mlab.figure(1)
mlab.clf()

#scene = engine.scenes[0]
scene.scene.x_plus_view()
scene.scene.parallel_projection = True

maxStrain = 5.0
minStrain = 0.0
PartStrain -= minStrain
PartStrain[PartStrain<0.0] = 0.0
PartStrain /= maxStrain-minStrain

PartStrainLogical = PartStrain>1.0
PartStrain[PartStrainLogical] = 1.0
#PartPattern += 1.0*PartStrainLogical

PartPattern += 1.0*PartStrain

#PartTimeLastPlastic[PartStrainLogical] = np.nan
print("render2")
glyph0 = mlab.points3d(x, PartX, PartY, PartPattern, mask_points=1,scale_mode="none",mode="cone",resolution=4)
glyph0.actor.property.lighting = False
glyph0.glyph.glyph_source.glyph_source.capping = False
glyph0.module_manager.scalar_lut_manager.lut.table = CMAP
#
module_manager = scene.children[0].children[0]
module_manager.scalar_lut_manager.show_legend = True


#mlab.draw()