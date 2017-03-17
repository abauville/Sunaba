#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD

import os
import sys
sys.path.insert(0, '../../src/UserInput')
#import json
#from InputDef import *
import InputDef as Input
import MaterialsDef as material
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos, atan
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
Setup.Description = "Setup to check the angle of decollement"



Numerics.phiMin = 1e-5
Numerics.phiMax = 0.9

Numerics.etaMin = 1e-4
Numerics.etaMax = 1e4

##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Sediment    = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Sediment}

PhaseRef = Sediment
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"

Sediment.vDiff = material.DiffusionCreep       ("Off")


Sediment.vDisl = material.DislocationCreep     (eta0=5E19, n=1)

StickyAir.rho0 = 0.0
#StickyAir.rho0 = 0000.00


Sediment.G = 1e10
StickyAir.G = 1e12

StickyAir.cohesion = .02e6/1.0#1.0*Sediment.cohesion
StickyAir.vDiff = material.DiffusionCreep(eta0=1E15)

## Main parameters for this setup
## =====================================

mu = 0.1
Sediment.frictionAngle = atan(mu)
slope = tan(0*pi/180)


HFac = 1.0
Hsed = HFac*1.0e3
Hnd = 1.0
Sediment.cohesion = Hsed / Hnd * (Sediment.rho0 * abs(Physics.gy))
StickyAir.cohesion = Sediment.cohesion/10.0#1.0*Sediment.cohesion





##              Grid

#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3

LWRatio = 5

Grid.xmin = HFac* -4.0e3*LWRatio
Grid.xmax = HFac*  0.0e3
Grid.ymin = HFac* 0.0e3
Grid.ymax = HFac* 4.0e3
Grid.nxC = 1/1*((128+64)*LWRatio) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 1/1*((128+64))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = True


##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin

i = 0
SedPhase = 1



Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W,"y","<",Grid.xmin,Grid.xmax)

#BCStokes.Sandbox_TopSeg00 = 0.525e3*HFac
#BCStokes.Sandbox_TopSeg01 = 0.475e3*HFac
BCStokes.Sandbox_TopSeg00 = 0.475e3*HFac
BCStokes.Sandbox_TopSeg01 = 0.525e3*HFac

##              Numerics
## =====================================
Numerics.nTimeSteps = 5000
Vboundary = 5 * cm/yr
BCStokes.backStrainRate = Vboundary/(Grid.xmin)
Numerics.CFL_fac_Stokes = 0.4
Numerics.CFL_fac_Darcy = 0.1
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 20
Numerics.maxNonLinearIter = 20

Numerics.absoluteTolerance = 5e-6





Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.9


##              Some info
## ======================================
Lc2 =  (Sediment.cohesion) / (Sediment.rho0 * abs(Physics.gy) * tan(Sediment.frictionAngle))  
Lc =  (Sediment.cohesion) / (Sediment.rho0 * abs(Physics.gy))  

print("Lc = " + str( Lc))
print("Lc2 = " + str( Lc2))
print("HnonDim = " + str(  Hsed/Lc ))
print("HnonDim2 = " + str( Hsed/Lc2))
print("Cohesion = " + str( Sediment.cohesion/1e6) + "MPa")



##              Output
## =====================================
Output.folder = "/Users/abauville/StokesFD_Output/WedgeSystematics_H_vs_mu/" + "H" +  ("%07d" % round(Hsed/Lc*1000)) + "_Mu" + ("%03d" % round(mu*100))
Output.phase            = True
Output.strainRate       = True
Output.P                = True
Output.khi              = True
Output.sigma_xx         = True
Output.sigma_xy         = True
Output.sigma_II         = True
Output.particles_pos    = True
Output.particles_posIni = True
Output.particles_phase  = True


#Vboundary = (abs(BCStokes.backStrainRate) * (Grid.xmax-Grid.xmin))
Output.timeFrequency    = (0.05*Hsed) / Vboundary
Numerics.maxTime        = 2*(Grid.xmax-Grid.xmin) / Vboundary




##                 BC
## =====================================
BCStokes.SetupType = "Sandbox"



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





##            Visualization
## =====================================
Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC


Visu.type = "StrainRate"
Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = Output.folder + "/00_Out_Visu/"
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

Visu.shiftFacY = -0.51


print("\n"*5)
CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

#StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1000.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))
print("nx = " + str(Grid.nxC) + ", ny = " + str(Grid.nyC))

print(Output.folder)

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 100.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Temperature.scale  = 1.0
Visu.colorMap.Temperature.center = 273.0/Char.temperature
Visu.colorMap.Temperature.max    = 1.0
Visu.colorMap.Porosity.log10on  = False
Visu.colorMap.Porosity.scale    = Sediment.phiIni/1.0
#Visu.colorMap.Porosity.center    = Sediment.phiIni/2.0
#Visu.colorMap.Porosity.center    = #0.1#Sediment.phiIni #ICDarcy.background
#Visu.colorMap.Porosity.max       = Sediment.phiIni+0.02 #Sediment.phiIni
#Visu.colorMap.Porosity.center = 0.0
Visu.colorMap.Porosity.max = 1.0

Visu.colorMap.Pressure.scale  = 50e6/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00
Visu.colorMap.CompactionPressure.scale  = 10e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.0
Visu.colorMap.FluidPressure.scale  = 50e6/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00




Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0





###          Write the Input file
### =====================================
Input.writeInputFile(Setup)
os.system("mkdir " + Output.folder)
os.system("mkdir " + Visu.outputFolder)

