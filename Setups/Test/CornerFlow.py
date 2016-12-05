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


##             Description
## =====================================
Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"

##      Declare singleton objects
## =====================================
Grid = input.Grid()
Numerics = input.Numerics()
Particles = input.Particles()
Physics = input.Physics(True)
Visu = input.Visu()
Char = input.Char()
BCStokes = input.BCStokes()
BCThermal = input.BCThermal()
Geometry = {}


##          Material properties
## =====================================
Phase0 = input.Material("Sediments")

#Phase1   = input.Material("Sediments")
Phase0.cohesion = 1e100
Phase0.n = 1.0;
Phase0.eta0 = 1e21
Phase0.G = 1e100


Backphi = 0.0001
RefPerm = 1e-20
Phase0.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))

Phase0.isRef = True
PhaseRef = Phase0;

MatProps = {'0': Phase0.__dict__}











##              Particles
## =====================================
Particles.noiseFactor = 0.95
Particles.nPCX = 3
Particles.nPCY = 3





##              Grid
## =====================================


Grid.xmin = -50.0e3
Grid.xmax =  80e3
Grid.ymin = -30e3
Grid.ymax = 0.0
Grid.nxC = 257#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 128#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = True


##                 BC
## =====================================
BCStokes.SetupType = "CornerFlow"
#BCStokes.SetupType = "PureShear"
#BCThermal.SetupType = "PureShear"
#BCStokes.SetupType = "SandBox"
#BCThermal.SetupType = "SandBox"

BCStokes.refValue       =  10.0 * cm/yr


BCThermal.TB = 1300.0
BCThermal.TT = 0.0


##              Non Dim
## =====================================
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
Char.set_based_on_corner_flow(PhaseRef,BCStokes,BCThermal,Physics,Grid)

BCStokes.backStrainRate = BCStokes.refValue / (Char.length)

##              Some info
## =====================================
Backphi = 0.001
RefPerm = Phase0.perm0*(Backphi * Backphi *Backphi  *  (1.0-Backphi)*(1.0-Backphi))
CompactionLength = sqrt(4/3*RefPerm/Physics.eta_f * (PhaseRef.eta0/Backphi))
DeltaRho = (Physics.rho_f-PhaseRef.rho0)
CompactionVelocity = abs(RefPerm/(Physics.eta_f*Backphi) * (1-Backphi) * DeltaRho*Physics.gy)
CompactionTime = CompactionLength/CompactionVelocity

print("CompactionLength: " + str(CompactionLength) + " m")
print("CompactionTime: " + str(CompactionTime/(3600*24*365*1e6)) + " Myr")


##            Visualization
## =====================================
Visu.showParticles = False
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Pressure"
Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
Visu.transparency = True

Visu.showGlyphs = True
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 1.0#BCStokes.refValue
Visu.glyphSamplingRateX = 8
Visu.glyphSamplingRateY = 8
Visu.showParticles = False




##              Numerics
## =====================================
Numerics.nTimeSteps = -1
Numerics.CFL_fac_Stokes = 0.75
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 25

Numerics.absoluteTolerance = 1e-5

Numerics.etaMin = 1e-5


#Visu.colorMap = vars(Visu.colorMap)




    

    
    

##          Write the input file
## =====================================
Visu.dictionarize()
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)


