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

PhaseRef = Matrix
PhaseRef.isRef = True

Matrix.name = "Matrix"
Inclusion.name = "Inclusion"


Matrix.cohesion = 1e6
Inclusion.cohesion = 1e6

Matrix.rho0     = 1000.0
Inclusion.rho0  = 1000.0

Matrix.G    = 5e8
Inclusion.G = 5e8




Plitho = Matrix.rho0 * abs(Physics.gy) * 1.0*1e3
Sigma_y = Matrix.cohesion*cos(Matrix.frictionAngle) + sin(Matrix.frictionAngle)*1.0*Plitho

BCStokes.backStrainRate = -1.0e-14

print("RefViscBrittle = %.2e Pa.s" % (Sigma_y/abs(BCStokes.backStrainRate)))
print("backStrainRate = %.2e, Sigma_y = %.2e MPa" % (BCStokes.backStrainRate, Sigma_y/1e6))


RefVisc =  (Sigma_y/abs(BCStokes.backStrainRate))




Matrix.vDisl = material.DislocationCreep    (eta0=RefVisc*100, n=1)
Inclusion.vDisl = material.DislocationCreep    (eta0=RefVisc*1e-2, n=1)

##              Grid
## =====================================
Grid.xmin = -1e3
Grid.xmax = +1e3
Grid.ymin = -1e3
Grid.ymax = +1e3
Grid.nxC = 128
Grid.nyC = 128

Grid.fixedBox = True



##              Numerics
## =====================================
Numerics.nTimeSteps = 50000

Numerics.CFL_fac_Stokes = 0.5
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 5
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 30

Numerics.absoluteTolerance = 5e-6



Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.1


Numerics.dtMaxwellFac_EP_ov_E  = .9;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .1;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress


#VatBound = 5.0 * cm/yr
#dx = (Grid.xmax-Grid.xmin)/Grid.nxC
#Numerics.dtVep = 1.0*Numerics.CFL_fac_Stokes*dx/abs(VatBound) 




#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##                 BC
## =====================================
#BCStokes.SetupType = "CornerFlow"


#BCStokes.refValue       = VatBound

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
Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)



##              Geometry
## =====================================


W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin
cx = .5*(Grid.xmax+Grid.xmin)
cy = .5*(Grid.ymax+Grid.ymin)
radius = W/16.0;


i = 0
phase = 1
Geometry["%05d_circle" % i] = input.Geom_Circle(phase,cx,cy,radius)



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

Visu.type = "StrainRate"
Visu.writeImages = False
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest2/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_Corner/G%.1e_phi%.1e/" % (Inclusion.G, Inclusion.phiIni)
Visu.transparency = True

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

Visu.colorMap.Stress.scale  = 200.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5


Visu.colorMap.Pressure.scale  = 10e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0


Visu.colorMap.Vorticity.max = 0.00001/yr /  (1.0/Char.time) # in rad/yr
###                 Info
### =====================================




###          Write the input file
### =====================================
input.writeInputFile(Setup)

os.system("mkdir " + Visu.outputFolder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")

