#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
import os
sys.path.insert(0, '../../src/UserInput')
#import json
#from InputDef import *
import InputDef as input
import MaterialsDef as material
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos
print("\n"*5)

##             Units
## =====================================
m       = 1.0
s       = 1.0
K       = 1.0
kg      = 1.0

cm      = 0.01      * m
km      = 1000.0    * m

mn      = 60        * s
hour    = 60        * mn
day     = 24        * hour
yr      = 365       * day
Myr     = 1e6       * yr






##      Declare singleton objects
## =====================================
Setup = input.Setup(isDimensional=True)
Grid = Setup.Grid
Numerics = Setup.Numerics
Particles = Setup.Particles
Physics = Setup.Physics
Visu = Setup.Visu
Char = Setup.Char
BCStokes = Setup.BC.Stokes
BCThermal = Setup.BC.Thermal
ICThermal = Setup.IC.Thermal
ICDarcy = Setup.IC.Darcy
MatProps = Setup.MatProps
Geometry = Setup.Geometry

## Description
## =====================================
Setup.Description = ""



Numerics.phiMin = 1e-5
Numerics.phiMax = 0.8

Numerics.etaMin = 1e-4
Numerics.etaMax = 1e5

##          Material properties
## =====================================
StickyAir   = input.Material("StickyAir")
Sediment    = input.Material("Wet_Quartzite")
Mantle      = input.Material("Dry_Olivine")

#Phase0 = input.Material()
#Phase1 = input.Material()
Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Mantle}

PhaseRef = Mantle
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Mantle.name = "Mantle"
Sediment.name = "Sediment"

Mantle.cohesion = 50e6
#StickyAir.cohesion = 1e6
#Sediment.cohesion = 10e6


#Mantle.vDisl.n = 1.0

Mantle.vPei.isActive = False
#
#
#StickyAir.phiIni = 0.9
#
#Sediment.phiIni = 0.25
#Mantle.phiIni = Numerics.phiMin
#
#Mantle.perm0 = 1e-9
#Sediment.perm0 = 1e-9
#StickyAir.perm0 = 1e-9
#
#
##Mantle.cohesion = 1e100
#StickyAir.rho0 = 1000.0
##StickyAir.G = 1e100
##Mantle.G = 1e100
##Sediment.G = 1e100
#
##Mantle.cohesion = 1e100
##Sediment.cohesion = 1e100
#
##Backphi = 0.0001
##RefPerm = 1e-18
##StickyAir.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))
##Mantle.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))
##Sediment.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))
#
##StickyAir.vDiff = material.DiffusionCreep(eta0=1e17)
#Sediment.G = 5e8
#Mantle.G = 5e9

Backphi = 0.0001
RefPerm = StickyAir.perm0*(Backphi * Backphi * Backphi  *  (1.0-Backphi)*(1.0-Backphi))


##              Grid
## =====================================
#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3
Grid.xmin = 1*-270.0e3
Grid.xmax = 2* 270e3
Grid.ymin = 1*-250e3
Grid.ymax = 1* 20.0e3
Grid.nxC = 1/1*256#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 1/1*128#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)


Grid.fixedBox = True



##              Numerics
## =====================================
Numerics.nTimeSteps = 50000
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac_Stokes = 0.5
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 5
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 50

Numerics.absoluteTolerance = 5e-6

#Numerics.use_dtMaxwellLimit = False



Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.1


Numerics.dtMaxwellFac_EP_ov_E  = .7;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .3;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
#Numerics.use_dtMaxwellLimit = False

#Numerics.maxTime = 8e5*yr

VatBound = 5.0 * cm/yr
dx = (Grid.xmax-Grid.xmin)/Grid.nxC
#Numerics.dtVep = 1.0*Numerics.CFL_fac_Stokes*dx/abs(VatBound) 




#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##                 BC
## =====================================
BCStokes.SetupType = "CornerFlow"


BCStokes.refValue       = VatBound

BCThermal.TB = 1300.0 + 273.0
BCThermal.TT = 0.0    + 273.0
#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);


##                 BC
## =====================================
#BCStokes.SetupType = "CornerFlow"
BCStokes.refValue       = 5.0 * cm/yr




##                 IC
## =====================================
#Setup.IC.Thermal = input.IC_HSC(age=100*Myr)
ICThermal.age = 25*Myr

ICDarcy.background = 0.0#Numerics.phiMin
ICDarcy.Amp = 0.0
ICDarcy.xc = Grid.xmin+(Grid.xmax-Grid.xmin)/2
ICDarcy.yc = Grid.ymin+(Grid.ymax-Grid.ymin)/2
ICDarcy.wx = (Grid.xmax-Grid.xmin)/16.0
ICDarcy.wy = (Grid.xmax-Grid.xmin)/16.0


##              Non Dim
## =====================================

L = (Grid.xmax-Grid.xmin)/2.0
BCStokes.backStrainRate = - BCStokes.refValue / L

Char.set_based_on_corner_flow(PhaseRef,BCStokes,BCThermal,Physics,Grid,L)
#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)



##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
Hsed = -0.0e3
H = -5e3

DetHL = 0.25*H
DetHR = 0.15*H

InterH = 0.15*H
InterY = Grid.ymin+0.6*H-InterH/2


Xbitonio = Grid.xmin + (Grid.xmax-Grid.xmin)/3.0
Lbitonio = 25e3


i = 0
MantlePhase = 2
SedPhase = 1
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)

#Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.0,Hsed,"y","<",Grid.xmin,Xbitonio)

#
#i+=1
#Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.0,1*Hsed   ,"y","<",Xbitonio + Lbitonio, Grid.xmax)
i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,H   ,"y","<",Grid.xmin,Xbitonio+Lbitonio)

i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,H+.2*H   ,"y","<",Xbitonio+Lbitonio,Grid.xmax)

i+=1 
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,-0.4,Hsed ,"y",">",Xbitonio,Xbitonio + Lbitonio)
i+=1 
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.8,Hsed+ -0.25*Lbitonio,"y",">",Xbitonio + Lbitonio,Xbitonio + Lbitonio*1.5)


i+=1 
Geometry["%05d_line" % i] = input.Geom_Line(0,0.0,Hsed,"y",">",Grid.xmin,Grid.xmax)

#Geometry["%05d_sine" % i] = input.Geom_Sine(MantlePhase,H, 2e3, 0.0, W/64, "y","<",Grid.xmin,Grid.xmax)
#i+=1
#Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,-0.4,0.4*(Xbitonio)+H ,"y","<",Xbitonio,Xbitonio + Lbitonio)

#i+=1
#Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,1*H   ,"y","<",Xbitonio + Lbitonio, Grid.xmax)

#Geometry["%05d_circle" % i] = (input.Geom_Circle(phase,0.0,0.0,0.33/2.0))

Numerics.stickyAirSwitchingDepth = -25e3;
Numerics.stickyAirSwitchPhaseTo  = 2;
Numerics.stickyAirSwitchPassiveTo  = 0;
Numerics.stickyAirTimeSwitchPassive = 250e3 * yr


##            Visualization
## =====================================
Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1.0 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Porosity"
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest2/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_Corner/G%.1e_phi%.1e/" % (Sediment.G, Sediment.phiIni)
Visu.transparency = True

Visu.showGlyphs = True
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.5/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.51
#Visu.shiftFacZ = 0.1
Visu.shaderFolder = "../Shaders/CornerFlow" # Relative path from the running folder (of StokesFD)

print("\n"*5)
CharExtra = input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/10000.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 200.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Temperature.scale  = 1
Visu.colorMap.Temperature.center = 273.0/Char.temperature
Visu.colorMap.Temperature.max    = 1.0
Visu.colorMap.Porosity.log10on  = False
Visu.colorMap.Porosity.scale    = Sediment.phiIni/1.0
Visu.colorMap.Porosity.center    = 1.0
#Visu.colorMap.Porosity.center    = #0.1#Sediment.phiIni #ICDarcy.background
#Visu.colorMap.Porosity.max       = Sediment.phiIni+0.02 #Sediment.phiIni
#Visu.colorMap.Porosity.center = 0.0
Visu.colorMap.Porosity.max = 1.5


Visu.colorMap.Pressure.scale  = 1000e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00
Visu.colorMap.CompactionPressure.scale  = 200e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.00
Visu.colorMap.FluidPressure.scale  = 1000e6/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00

Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0

Visu.colorMap.Permeability.scale = RefPerm/Physics.eta_f / (Char.length*Char.length / (Char.mass/Char.length/Char.time) )
Visu.colorMap.Permeability.max = 10.0

Visu.colorMap.Vorticity.max = 0.00001/yr /  (1.0/Char.time) # in rad/yr
###                 Info
### =====================================




###          Write the input file
### =====================================
input.writeInputFile(Setup)

os.system("mkdir " + Visu.outputFolder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")

