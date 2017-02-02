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
import InputDef as Input
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

Pa      = 1.0
MPa     = 1e6

deg     = pi/180




##      Declare singleton objects
## =====================================
Setup = Input.Setup(isDimensional=True)
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
Output = Setup.Output

## Description
## =====================================
Setup.Description = "Angle of shear bands benchmark, based on Kaus, 2010 (doi:10.1016/j.tecto.2009.08.042)"



Numerics.phiMin = 1e-5
Numerics.phiMax = 0.9

Numerics.etaMin = 1e-6

##          Material properties
## =====================================
Matrix      = Input.Material("Sediments")
Inclusion   = Input.Material("Sediments")
Basement    = Input.Material("Sediments")

Setup.MatProps = {"0":Matrix,"1":Inclusion}

PhaseRef = Matrix
PhaseRef.isRef = True

Matrix.name     = "Sediment"
Inclusion.name  = "Inclusion"

Matrix.vDiff = material.DiffusionCreep       ("Off")
Inclusion.vDiff = material.DiffusionCreep       ("Off")

Matrix.vDisl    = material.DislocationCreep     (eta0=1E25, n=1)
Inclusion.vDisl = material.DislocationCreep     (eta0=1E20, n=1)

Matrix.rho0     = 2700  * kg/(m**3)
Inclusion.rho0  = 2700  * kg/(m**3)

Matrix.cohesion     = 40    * MPa
Inclusion.cohesion  = 40    * MPa

Matrix.frictionAngle    = 30 * deg
Inclusion.frictionAngle = 00 * deg

Matrix.G                = 5e10 * Pa
Inclusion.G             = 5e10 * Pa

##              Grid
## =====================================
#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3
RFac = 4; # Resolution Factor
HFac = 1/10

W = HFac * 40 * km
H = HFac * 10 * km
d = HFac*800*4          # inclusion size

Grid.xmin = -W/2
Grid.xmax = +W/2
Grid.ymin =  0.0
Grid.ymax =  H
Grid.nxC = RFac*128
Grid.nyC = RFac*32

Grid.fixedBox = False



##              Numerics
## =====================================
Numerics.nTimeSteps = 2
BCStokes.backStrainRate = -1.0e-15
Numerics.CFL_fac_Stokes = 0.5
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 5
Numerics.maxNonLinearIter = 250

Numerics.absoluteTolerance = 1e-5





Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.5


##              Output
## =====================================
Output.folder = "/Users/abauville/GoogleDrive/StokesFD_Output/OutputTest"
#Output.khi = True
#Output.strainRate = True
#Output.frequency = Numerics.nTimeSteps-1





##                 BC
## =====================================
#BCStokes.SetupType = "Sandbox"


#BCStokes.refValue       = 1.0 * cm/yr / 1.0


#BCThermal.TB = 30.0 + 273.0
#BCThermal.TT = 0.0    + 273.0


#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);

##                 IC
## =====================================
#Setup.IC.Thermal = Input.IC_HSC(age=100*Myr)
ICThermal.age = 100*Myr

ICDarcy.background = 0.0#Numerics.phiMin
ICDarcy.Amp = 0.0
ICDarcy.xc = Grid.xmin+(Grid.xmax-Grid.xmin)/2
ICDarcy.yc = Grid.ymin+(Grid.ymax-Grid.ymin)/2
ICDarcy.wx = (Grid.xmax-Grid.xmin)/16.0
ICDarcy.wy = (Grid.xmax-Grid.xmin)/16.0

##              Non Dim
## =====================================

#L = (Grid.xmax-Grid.xmin)/2.0
L = (Grid.ymax-Grid.ymin)/2.0
#BCStokes.backStrainRate = - BCStokes.refValue / L

#Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid)
Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)



##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin


inclusion_w = d
inclusion_h = d/2


slope = tan(0*pi/180)

i = 0
InclusionPhase = 1
Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,slope,inclusion_h,"y","<",Grid.xmin + W/2 - inclusion_w/2,Grid.xmin + W/2 + inclusion_w/2)


##            Visualization
## =====================================
Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = False
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC


Visu.type = "StrainRate"
Visu.writeImages = False
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_Sandbox/"
Visu.transparency = True

Visu.showGlyphs = False
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.1/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 0.5 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = 0.0


print("\n"*5)
CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 100.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.0


Visu.colorMap.Pressure.scale  = 250e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00


Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0




###          Write the Input file
### =====================================
Input.writeInputFile(Setup)


