#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
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

##          Material properties
## =====================================
StickyAir   = input.Material("StickyAir")
Mantle      = input.Material("Dry_Olivine")
Sediment    = input.Material("Wet_Quartzite")
#Phase0 = input.Material()
#Phase1 = input.Material()
Setup.MatProps = {"0":StickyAir,"1":Mantle,"2":Sediment}

PhaseRef = Mantle
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Mantle.name = "Mantle"
Sediment.name = "Sediment"

Mantle.cohesion = 50e6

#Mantle.vDisl.n = 1.0

Mantle.vPei.isActive = False


StickyAir.phiIni = 0.9
Mantle.phiIni = 0.00001
Sediment.phiIni = 0.1

Mantle.perm0 = 1e-6
Sediment.perm0 = 1e-6
StickyAir.perm0 = 1e-6


#Mantle.cohesion = 1e100
StickyAir.rho0 = 1000.0
#StickyAir.G = 1e100
#Mantle.G = 1e100
#Sediment.G = 1e100

#Mantle.cohesion = 1e100
#Sediment.cohesion = 1e100

Backphi = 0.0001
RefPerm = StickyAir.perm0*(Backphi * Backphi * Backphi  *  (1.0-Backphi)*(1.0-Backphi))
#RefPerm = 1e-18
#StickyAir.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))
#Mantle.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))
#Sediment.perm0 = RefPerm/(Backphi * Backphi * Backphi  /  (1.0-Backphi)*(1.0-Backphi))



##              Grid
## =====================================
#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3
Grid.xmin = 1*-125.0e3
Grid.xmax = 1* 125e3
Grid.ymin = 1*-140e3
Grid.ymax = 1* 10.0e3
Grid.nxC = 1/1*1024#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 1/1*512#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = False



##              Numerics
## =====================================
Numerics.nTimeSteps = 10
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac_Stokes = 0.1
Numerics.CFL_fac_Darcy = 0.1
Numerics.CFL_fac_Thermal = 1.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 100

Numerics.absoluteTolerance = 1e-5





Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.1






#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##                 BC
## =====================================
#BCStokes.SetupType = "CornerFlow"


BCStokes.refValue       = 1.0 * cm/yr


BCThermal.TB = 1300.0 + 273.0
BCThermal.TT = 0.0    + 273.0
#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);



##                 IC
## =====================================
#Setup.IC.Thermal = input.IC_HSC(age=100*Myr)
ICThermal.age = 180*Myr

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
H = -7e3

DetHL = 0.25*H
DetHR = 0.15*H

InterH = 0.15*H
InterY = Grid.ymin+0.6*H-InterH/2


XcIsland = 0.0 # Grid.xmin + (Grid.xmax-Grid.xmin)/3.0
LIsland = 60e3
HIsland = 5e3
Hroot   = 7e3
slope = HIsland/(LIsland/2)
slopeRoot = Hroot/(LIsland/2)

i = 0
MantlePhase = 1
SedPhase = 2
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)

#Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,0.0,Hsed,"y","<",Grid.xmin,Xbitonio)

#



i+=1
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,slope,-slope*(XcIsland- LIsland/2)    ,"y","<",XcIsland - LIsland/2, XcIsland)

i+=1
Geometry["%05d_line" % i] = input.Geom_Line(SedPhase,-slope,slope*(XcIsland)+HIsland    ,"y","<",XcIsland, XcIsland + LIsland/2)







i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,-slopeRoot,slopeRoot*(XcIsland- LIsland/2)+H    ,"y","<",XcIsland - LIsland/2, XcIsland)

i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,+slopeRoot,-slopeRoot*(XcIsland)+H-Hroot    ,"y","<",XcIsland, XcIsland + LIsland/2)


i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,H   ,"y","<",Grid.xmin,XcIsland - LIsland/2)

i+=1
Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,H   ,"y","<",XcIsland + LIsland/2,Grid.xmax)


#Geometry["%05d_sine" % i] = input.Geom_Sine(MantlePhase,H, 2e3, 0.0, W/64, "y","<",Grid.xmin,Grid.xmax)
#i+=1
#Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,-0.4,0.4*(Xbitonio)+H ,"y","<",Xbitonio,Xbitonio + Lbitonio)

#i+=1
#Geometry["%05d_line" % i] = input.Geom_Line(MantlePhase,0.0,1*H   ,"y","<",Xbitonio + Lbitonio, Grid.xmax)

#Geometry["%05d_circle" % i] = (input.Geom_Circle(phase,0.0,0.0,0.33/2.0))

Numerics.stickyAirSwitchingDepth = -50e3;
Numerics.stickyAirSwitchPhaseTo  = 2;
Numerics.stickyAirSwitchPassiveTo  = 0;
Numerics.stickyAirTimeSwitchPassive = 250e3 * yr


##            Visualization
## =====================================
Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/8
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1.0 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "StrainRate"
Visu.writeImages = False
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/OutputNew/"
Visu.transparency = True

Visu.showGlyphs = False
Visu.glyphType = "DarcyGradient"
Visu.glyphMeshType = "Triangle"
#Visu.glyphScale = 0.5/(BCStokes.refValue/(Char.length/Char.time))
Visu.glyphScale = 10.0/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 2 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.6


print("\n"*5)
CharExtra = input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1000.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 100.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.0
Visu.colorMap.Temperature.scale  = 1.0
Visu.colorMap.Temperature.center = 273.0/Char.temperature
Visu.colorMap.Temperature.max    = 1.0
Visu.colorMap.Porosity.log10on  = True
Visu.colorMap.Porosity.scale    = 0.001
#Visu.colorMap.Porosity.center    = #0.1#Sediment.phiIni #ICDarcy.background
#Visu.colorMap.Porosity.max       = Sediment.phiIni+0.02 #Sediment.phiIni
Visu.colorMap.Porosity.max = 2.0


Visu.colorMap.Pressure.scale  = 1000e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00
Visu.colorMap.CompactionPressure.scale  = 1000e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.00
Visu.colorMap.FluidPressure.scale  = 25e6/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00

Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 10.0
Visu.colorMap.Khib.max = 10.0

Visu.colorMap.Permeability.scale = RefPerm/Physics.eta_f / (Char.length*Char.length / (Char.mass/Char.length/Char.time) )
Visu.colorMap.Permeability.max = 10.0



###          Write the input file
### =====================================
input.writeInputFile(Setup)


