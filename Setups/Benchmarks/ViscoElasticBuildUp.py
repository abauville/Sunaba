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
Setup.Description = ""



Numerics.phiMin = 1e-5
Numerics.phiMax = 0.8

Numerics.etaMin = 1e-4
Numerics.etaMax = 1e5

##          Material properties
## =====================================
Matrix      = Input.Material("Sediments")
Inclusion   = Input.Material("Sediments")

Setup.MatProps = {"0":Matrix,"1":Inclusion}

PhaseRef = Matrix
PhaseRef.isRef = True

Matrix.name = "Matrix"
Inclusion.name = "Inclusion"


Matrix.rho0     = 1000.0
Inclusion.rho0  = 1000.0

Matrix.G    = 5e10
Inclusion.G = 5e10

Matrix.cohesion     = 1e100
Inclusion.cohesion  = 1e100

#Matrix.use_dtMaxwellLimit = True

BCStokes.backStrainRate = -1.0e-15


RefVisc =  1e20

Matrix.vDiff    = material.DiffusionCreep    (eta0=RefVisc*1)
Inclusion.vDiff = material.DiffusionCreep    (eta0=RefVisc*1e-3)
Matrix.vDisl    = material.DislocationCreep   ("Off")
Inclusion.vDisl = material.DislocationCreep   ("Off")

##              Grid
## =====================================
Grid.xmin = -1e0
Grid.xmax = +1e0
Grid.ymin = -1e0
Grid.ymax = +1e0
Grid.nxC = 32
Grid.nyC = 32

Grid.fixedBox = False



##              Numerics
## =====================================
Numerics.nTimeSteps = 20

Numerics.CFL_fac_Stokes = 0.5
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1

Numerics.absoluteTolerance = 1e-7


Particles.nPCX = 3
Particles.nPCY = 3
Particles.noiseFactor = 0.0


Numerics.dtMaxwellFac_EP_ov_E  = .5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = False

# Fix the timestep
Numerics.dtMin = 20 * yr
Numerics.dtMax = Numerics.dtMin


Physics.gy = 0.0
Physics.gx = 0.0

##                 BC
## =====================================
#BCStokes.SetupType = "CornerFlow"


#BCStokes.refValue       = VatBound

#BCThermal.TB = 1300.0 + 273.0
#BCThermal.TT = 0.0    + 273.0


#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);


##               Output
## =====================================
Output.folder = "/Users/abauville/StokesFD_Output/ViscoElasticBuildUp"
os.system("mkdir " + Output.folder)
Output.sigma_xx0 = True

##              Non Dim
## =====================================
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
#Geometry["%05d_circle" % i] = Input.Geom_Circle(phase,cx,cy,radius)



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
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest2/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/FunOutput/"
Visu.transparency = False

Visu.height = 1 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.0
Visu.shiftFacZ = 0.1

Visu.colorMap.Stress.scale  = 1.0
Visu.colorMap.Stress.center = 0
Visu.colorMap.Stress.max    = 1.0

Visu.closeAtTheEndOfSimulation = True


###                 Info
### =====================================




###          Write the Input file
### =====================================
Input.writeInputFile(Setup)

if (Visu.writeImages):
    os.system("mkdir " + Visu.outputFolder)
    
    
###                 Run
### =====================================
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ../Input.json")
#import ViscoElasticBuildUp_PostProc



###          Post-processing
### =====================================























































