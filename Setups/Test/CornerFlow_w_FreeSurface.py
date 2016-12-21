#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
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
MatProps = Setup.MatProps
Geometry = Setup.Geometry

## Description
## =====================================
Setup.Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"


##          Material properties
## =====================================
StickyAir   = input.Material("StickyAir")
Mantle      = input.Material("Dry_Olivine")
#Phase0 = input.Material()
#Phase1 = input.Material()
Setup.MatProps = {"0":StickyAir,"1":Mantle}

PhaseRef = Mantle
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Mantle.name = "Mantle"
#Mantle.vDisl.n = 1.0

Mantle.vPei.isActive = False

#Mantle.cohesion = 1e100






Backphi = 0.0001
RefPerm = 1e-20
StickyAir.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))
Mantle.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))




##              Grid
## =====================================
Grid.xmin = -275.0e3
Grid.xmax =  1000e3
Grid.ymin = -500e3
Grid.ymax =  10.0e3
Grid.nxC = 3#257#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 8#128#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = True




##              Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac_Stokes = 0.5
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 10

Numerics.absoluteTolerance = 1e-5



Visu.showParticles = False


Particles.nPCX = 4
Particles.nPCY = 4


Particles.noiseFactor = 0.0

#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##                 BC
## =====================================
BCStokes.SetupType = "CornerFlow"
#BCStokes.SetupType = "PureShear"
#BCThermal.SetupType = "PureShear"
#BCStokes.SetupType = "SandBox"
#BCThermal.SetupType = "SandBox"
#BCThermal.SetupType = "TT_TBExternal_LRNoFlux"

BCStokes.refValue       =  10.0 * cm/yr


BCThermal.TB = 1300.0 + 273.0
BCThermal.TT = 0.0    + 273.0
#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);


##              Non Dim
## =====================================

L = (Grid.xmax-Grid.xmin)/2.0
BCStokes.backStrainRate = - BCStokes.refValue / L

Char.set_based_on_corner_flow(PhaseRef,BCStokes,BCThermal,Physics,Grid,L)
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)

##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = 0e3

DetHL = 0.25*H
DetHR = 0.15*H

InterH = 0.15*H
InterY = Grid.ymin+0.6*H-InterH/2


i = 0
phase = 1
Geometry["%05d_line" % i] = (input.Geom_Line(phase,0.0,H,"y","<",Grid.xmin,Grid.xmax))
#Geometry["%05d_circle" % i] = (input.Geom_Circle(phase,0.0,0.0,0.33/2.0))




##            Visualization
## =====================================
Visu.showParticles = False
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Temperature"
Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
Visu.transparency = True

Visu.showGlyphs = True
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.5 * 1.0/(BCStokes.refValue/(Char.length/Char.time))
Visu.glyphSamplingRateX = 8
Visu.glyphSamplingRateY = 8

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Viscosity"
#Visu.filter = "Linear"
Visu.filter = "Nearest"


print("\n"*5)
CharExtra = input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 1.0
Visu.colorMap.Stress.center = 1.0
Visu.colorMap.Stress.max    = 1.75
Visu.colorMap.Pressure.scale  = RefP/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.75
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.0
Visu.colorMap.Temperature.scale  = 1.0
Visu.colorMap.Temperature.center = 0.0
Visu.colorMap.Temperature.max    = 1.0

###          Write the input file
### =====================================
input.writeInputFile(Setup)


