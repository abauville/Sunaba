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
from math import pi, sqrt, tan, sin, cos, log10
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
Setup.Description = "Viagra test"



##          Material properties
## =====================================
Matrix    = Input.Material("Sediments")
Block    = Input.Material("Sediments")



Setup.MatProps = {"0":Matrix,"1":Block}




PhaseRef = Matrix
PhaseRef.isRef = True

Matrix.name = "Matrix"
Block.name = "Block"

Matrix.vDiff = material.DiffusionCreep       ("Off")
Block.vDiff = material.DiffusionCreep       ("Off")


## Main parameters for this setup
## =====================================

Matrix.frictionAngle  = 30/180*pi
Block.frictionAngle  = Matrix.frictionAngle


Matrix.cohesion =  1e100
Block.cohesion = 1e100

Matrix.rho0 = 1.0
Block.rho0 = 4000.0

Matrix.G    = 1e15
Block.G    = 1e10


Matrix.vDisl = material.DislocationCreep     (eta0=1e21, n=1)
Block.vDisl = material.DislocationCreep     (eta0=1e25, n=1)

#Matrix.vDisl = material.DislocationCreep     ('Off')
#Block.vDisl = material.DislocationCreep     ('Off')
#
#Matrix.vDiff = material.DiffusionCreep     (eta0=1e21)
#Block.vDiff = material.DiffusionCreep     (eta0=1e25)

Physics.gy = -10.0



Grid.xmin = 0.0
Grid.xmax = 1000.0e3
Grid.ymin = 0.0
Grid.ymax = 1000.0e3
Grid.nxC = 64
Grid.nyC = 64

Grid.fixedBox = False


Numerics.etaMax = 1e10;


##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin
i = 0
MatrixPhase = 0
BlockPhase = 1


Geometry["%05d_line" % i] = Input.Geom_Line(BlockPhase,0.0,800e3,"y","<",Grid.xmin,800e3)
i+=1
Geometry["%05d_line" % i] = Input.Geom_Line(MatrixPhase,0.0,200e3,"y","<",Grid.xmin,800e3)

##              Numerics
## =====================================
Numerics.nTimeSteps = 1000
Numerics.CFL_fac_Stokes = .1
Numerics.CFL_fac_Darcy = 1000.0
Numerics.CFL_fac_Thermal = 10000.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 1
Numerics.dtAlphaCorr = .3
Numerics.absoluteTolerance = 1e-14


Numerics.dtMaxwellFac_EP_ov_E  = .5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = False

#Numerics.maxTime = 8e5*yr
dx = (Grid.xmax-Grid.xmin)/Grid.nxC
Numerics.dtVep = 1.0*Numerics.CFL_fac_Stokes*dx/abs(.1*cm/yr) 


Numerics.dtMax = 10*yr
Numerics.dtMin = Numerics.dtMax

Particles.nPCX = 12
Particles.nPCY = 12
Particles.noiseFactor = 0.0


###              Output
### =====================================
#Output.folder = "/Users/abauville/StokesFD_Output/TestDarcy"
#Output.phase = True
#Output.strainRate = True
#Output.sigma_II = True
#Output.khi = True
#Output.particles_pos = True
#Output.particles_posIni = True
#Output.particles_phase = True
#Output.Pc = True
#Output.Pf = True
#Output.porosity = True
#Output.strainRate = True

#
#
#Output.frequency = 5




##                 BC
## =====================================
BCStokes.SetupType = "FixedLeftWall"

#BCStokes.Sandbox_NoSlipWall = True

#BCStokes.refValue       = 1.0 * cm/yr / 1.0


#BCThermal.TB = 30.0 + 273.0
#BCThermal.TT = 0.0  + 273.0


#Mantle.vDisl.E = 0.0
#Mantle.vDiff.E = 0.0
#BCThermal.DeltaL = 1000e3+(Grid.ymin);


##              Non Dim
## =====================================

#L = (Grid.xmax-Grid.xmin)/2.0
L = (Grid.ymax-Grid.ymin)/2.0
BCStokes.backStrainRate = 0.0

Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid,(Grid.ymax-Grid.ymin)/2.0)

#Char.length = Hsed/8.0
#
#Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
#CharVisc = RefVisc
#  
#CharStress  = PhaseRef.rho0*abs(Physics.gy)*Char.length
#
#
#Char.time   = CharVisc/CharStress
#Char.mass   = CharStress*Char.time*Char.time*Char.length

##            Visualization
## =====================================
Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

#Visu.shaderFolder = "../Shaders/Sandbox" # Relative path from the running folder (of StokesFD)


Visu.type = "Stress"
Visu.closeAtTheEndOfSimulation = False
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_SandboxNew/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_Test3/"
Visu.transparency = False

Visu.showGlyphs = False
Visu.glyphType = "DarcyGradient"
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.0001/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 0.75 * Visu.height
Visu.width = 0.75* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.0
Visu.shiftFacZ = 0.1


print("\n"*5)
CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
SedVisc = Matrix.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
BaseVisc = Block.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

#Visu.colorMap.Stress.scale  = 10000e6/CharExtra.stress
#Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
#Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(1e-13/(1.0/Char.time))#abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Temperature.scale  = 1.0
Visu.colorMap.Temperature.center = 273.0/Char.temperature
Visu.colorMap.Temperature.max    = 1.0
Visu.colorMap.Porosity.log10on  = False
Visu.colorMap.Porosity.scale    = Matrix.phiIni/1.0
#Visu.colorMap.Porosity.center    = Matrix.phiIni/2.0
#Visu.colorMap.Porosity.center    = #0.1#Matrix.phiIni #ICDarcy.background
#Visu.colorMap.Porosity.max       = Matrix.phiIni+0.02 #Matrix.phiIni
#Visu.colorMap.Porosity.center = 0.0
Visu.colorMap.Porosity.max = 1.0

Visu.colorMap.Pressure.scale  = 10000e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00
Visu.colorMap.CompactionPressure.scale  = 5e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.0
Visu.colorMap.FluidPressure.scale  = 50e6/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00




Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0

Visu.colorMap.Velocity.scale = 5.0 * (cm/yr) / (Char.length/Char.time)

Visu.colorMap.Vorticity.max = 0.00005/yr /  (1.0/Char.time) # in rad/yr

Visu.colorMap.ShearModulus.center = log10(Block.G/CharExtra.stress)
Visu.colorMap.ShearModulus.max =  2.5*Visu.colorMap.ShearModulus.center


Visu.colorMap.EffectiveViscosity.scale = RefVisc/1e2 / CharExtra.visc
Visu.colorMap.EffectiveViscosity.max = 2.0

Visu.colorMap.Velocity.log10on = True
Visu.colorMap.Velocity.scale = 100



Visu.colorMap.Stress.scale  = 100e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0

##              Some info
## ======================================
print("Lc = " + str(  (Matrix.cohesion*cos(Matrix.frictionAngle)) / (Matrix.rho0 * abs(Physics.gy) * sin(Matrix.frictionAngle))  ))
#CompactionLength = sqrt(4.0/3.0*Matrix.perm0/Physics.eta_f * (Physics->eta[iCell]/phi));



###          Write the Input file
### =====================================
Input.writeInputFile(Setup)


