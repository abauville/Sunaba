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
Matrix      = input.Material("Sediments")
Inclusion   = input.Material("Sediments")

Setup.MatProps = {"0":Matrix,"1":Inclusion}

PhaseRef = Inclusion
PhaseRef.isRef = True

Matrix.name = "Matrix"
Inclusion.name = "Inclusion"


Matrix.cohesion = 1e6 * 1e100
Inclusion.cohesion = 1e6 * 1e100

Matrix.rho0     = 1.0
Inclusion.rho0  = 1.0

Matrix.G    = 1e100
Inclusion.G = 1e100

Matrix.use_dtMaxwellLimit = True



Plitho = Matrix.rho0 * abs(Physics.gy) * 1.0*1e3 * 100.0
Sigma_y = Matrix.cohesion*cos(Matrix.frictionAngle) + sin(Matrix.frictionAngle)*1.0*Plitho



print("RefViscBrittle = %.2e Pa.s" % (Sigma_y/abs(BCStokes.backStrainRate)))
print("backStrainRate = %.2e, Sigma_y = %.2e MPa" % (BCStokes.backStrainRate, Sigma_y/1e6))

Matrix.use_dtMaxwellLimit = True


RefVisc =  1.0*(Sigma_y/abs(BCStokes.backStrainRate))



Physics.gy = 0.0

#Matrix.vDisl    = material.DislocationCreep    (eta0=RefVisc*10, n=1)
#Inclusion.vDisl = material.DislocationCreep    (eta0=RefVisc*1e-1, n=1)
#Matrix.vDiff    = material.DiffusionCreep   ("Off")
#Inclusion.vDiff = material.DiffusionCreep   ("Off")



VatBound = 10 * m/s
InclusionRadius = 0.5*m

Re = 52 # Reynolds number

RefVisc = VatBound*2*InclusionRadius/Re
Matrix.vDiff    = material.DiffusionCreep    (eta0=RefVisc*1)
Inclusion.vDiff = material.DiffusionCreep    (eta0=RefVisc*1)
Matrix.vDisl    = material.DislocationCreep   ("Off")
Inclusion.vDisl = material.DislocationCreep   ("Off")


##              Grid
## =====================================
Grid.xmin = -3.0*m
Grid.xmax = Grid.xmin + 12.0*m
Grid.ymin = -4.0*m
Grid.ymax = +4.0*m

Grid.nyC = 128
Grid.nxC = round(Grid.nyC * (Grid.xmax-Grid.xmin)/(Grid.ymax-Grid.ymin))

Grid.fixedBox = True





##              Numerics
## =====================================
Numerics.nTimeSteps = 2000

Numerics.CFL_fac_Stokes = 0.4
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1

Numerics.absoluteTolerance = 1e-7



Particles.nPCX = 3
Particles.nPCY = 3
Particles.noiseFactor = 0.9


Numerics.dtMaxwellFac_EP_ov_E  = .5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = False
Physics.gy = 0.0
Physics.gx = 0.0
#VatBound = 5.0 * cm/yr
#dx = (Grid.xmax-Grid.xmin)/Grid.nxC
#Numerics.dtVep = 1.0*Numerics.CFL_fac_Stokes*dx/abs(VatBound) 




#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##                 BC
## =====================================
BCStokes.SetupType = "WindTunnel"


BCStokes.refValue       = VatBound
BCStokes.backStrainRate = BCStokes.refValue/(Grid.xmax-Grid.xmin)

#BCThermal.TB = 1300.0 + 273.0
#BCThermal.TT = 0.0    + 273.0


#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);


#
###                 IC
### =====================================
##Setup.IC.Thermal = input.IC_HSC(age=100*Myr)
#ICThermal.age = 25*Myr
#
#ICDarcy.background = 0.0#Numerics.phiMin
#ICDarcy.Amp = 0.0
#ICDarcy.xc = Grid.xmin+(Grid.xmax-Grid.xmin)/2
#ICDarcy.yc = Grid.ymin+(Grid.ymax-Grid.ymin)/2
#ICDarcy.wx = (Grid.xmax-Grid.xmin)/16.0
#ICDarcy.wy = (Grid.xmax-Grid.xmin)/16.0


##              Non Dim
## =====================================
#L = (Grid.xmax-Grid.xmin)/2.0
#BCStokes.backStrainRate = - BCStokes.refValue / L

#Char.set_based_on_corner_flow(PhaseRef,BCStokes,BCThermal,Physics,Grid,L)
#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)

Char.length = 2.0*InclusionRadius
Char.time = Char.length/VatBound
Char.mass = 1.0
Char.temperature = 1.0


##              Geometry
## =====================================


W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin
cx = 0.0#.5*(Grid.xmax+Grid.xmin)
cy = 0.0#.5*(Grid.ymax+Grid.ymin)



i = 0
phase = 1
Geometry["%05d_circle" % i] = input.Geom_Circle(phase,cx,cy,InclusionRadius)



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

Visu.type = "Velocity"
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest2/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/FunOutput/"
Visu.transparency = False

#Visu.showGlyphs = True
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.5/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.0
Visu.shiftFacZ = 0.1
#Visu.shaderFolder = "../Shaders/CornerFlow" # Relative path from the running folder (of StokesFD)

print("\n"*5)
CharExtra = input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
SedVisc = Inclusion.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))



MatrixVisc = Matrix.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Inclusion Visc = %.2e" % SedVisc)
print("MatrixVisc = %.2e" % MatrixVisc)


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 1 #e6/CharExtra.stress
Visu.colorMap.Stress.center = 0 #*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5


Visu.colorMap.Pressure.scale  = 1000#1e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0


Visu.colorMap.POvPlitho.log10on = True
Visu.colorMap.POvPlitho.center = 0.0
Visu.colorMap.POvPlitho.max = 1.0

Visu.colorMap.Vorticity.max = 0.00001/yr /  (1.0/Char.time) # in rad/yr
Visu.colorMap.Velocity.max = 2.0
Visu.colorMap.Velocity.log10on = False
###                 Info
### =====================================




###          Write the input file
### =====================================
input.writeInputFile(Setup)

os.system("mkdir " + Visu.outputFolder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")

