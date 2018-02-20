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
MPa     = 1e6       * Pa



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

Numerics.etaMin = 1e-6

##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Sediment    = Input.Material("Sand")
Basement    = Input.Material("Sand")

Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Basement}

PhaseRef = Sediment
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"
Basement.name = "Basement"

#StickyAir.rho0 = 1000.0
StickyAir.rho0 = 0000.00

Sediment.vDiff = material.DiffusionCreep     (eta0=1e7)
StickyAir.vDiff = material.DiffusionCreep     (eta0=1e3)
Basement.vDiff = material.DiffusionCreep     (eta0=1e7)
#Sediment.vDiff = DiffusionCreep     (eta0=1e10)


#Sediment.frictionAngle = 30/180*pi
Basement.frictionAngle = 10/180*pi

Sediment.perm0 = 1e-8

Sediment.rho0 = 1600
Basement.rho0 = 1600

Sediment.G = 1e10
Basement.G = 1e10
StickyAir.G = 1e10

StickyAir.cohesion = 1.0*Sediment.cohesion/1.0


##              Grid
## =====================================
#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3
HFac = 1.0;
LWRatio = 2

#Grid.xmin = HFac* -12.0 * cm
#Grid.xmax = HFac* 0.0
Grid.xmin = HFac* 0.0  * cm
Grid.xmax = LWRatio*HFac* 12.0 * cm
Grid.ymin = HFac* 0.0e3
Grid.ymax = HFac* 4 * cm
Grid.nxC = LWRatio*2/1*(64*LWRatio) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nyC = 2/1*(64)#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

#Grid.fixedBox = True


##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin
Hsed = HFac*3.0 * cm
Hbase = HFac*0.1 * cm

#slope = tan(20*pi/180)
BoxTilt = -5.0 * pi/180
slope = -BoxTilt #tan(0*pi/180)

Physics.gx = -9.81*sin(BoxTilt);
Physics.gy = -9.81*cos(BoxTilt);

i = 0
SedPhase = 1
BasementPhase = 2
Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W,"y","<",Grid.xmin,Grid.xmax)
#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,0.0,Hbase,"y","<",Grid.xmax-8.24*cm,Grid.xmax)
#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,0.0,Grid.xmax-Hbase,"x",">",Grid.ymin,Hsed)


BCStokes.Sandbox_TopSeg00 = 0.0e3*HFac
BCStokes.Sandbox_TopSeg01 = 0.0e3*HFac

##              Numerics
## =====================================
Numerics.nTimeSteps = -200
Numerics.CFL_fac_Stokes = 0.5
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 5
Numerics.maxNonLinearIter = 10

Numerics.absoluteTolerance = 1e-8
Numerics.relativeTolerance = 1e-5


Numerics.dtMin = 15.0*s
Numerics.dtMax = 15.0*s


Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.5

Numerics.yieldComputationType = 0
Numerics.stressSubGridDiffFac = 0.0

##              Output
## =====================================
Output.folder = "/Users/abauville/GoogleDrive/StokesFD_Output/OutputTest"
#Output.khi = True
#Output.strainRate = True
#Output.frequency = Numerics.nTimeSteps-1





##                 BC
## =====================================
BCStokes.SetupType = "Sandbox"



BCStokes.refValue       = -2.5 * cm/hour
BCStokes.backStrainRate = BCStokes.refValue/(Grid.xmax-Grid.xmin)
#BCStokes.Sandbox_NoSlipWall = True




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

Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid)
#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)





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
#Visu.transparency = True

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

Visu.shiftFacY = -0.0


print("\n"*5)
CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
BaseVisc = Basement.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

#StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/10000.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)
print("BaseVisc = %.2e" %  BaseVisc)


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

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
Visu.colorMap.Porosity.center    = Sediment.phiIni/2.0
#Visu.colorMap.Porosity.center    = #0.1#Sediment.phiIni #ICDarcy.background
#Visu.colorMap.Porosity.max       = Sediment.phiIni+0.02 #Sediment.phiIni
Visu.colorMap.Porosity.max = 1.0


Visu.colorMap.Pressure.scale  = 600 * Pa/CharExtra.stress
Visu.colorMap.Pressure.center = 0.0
Visu.colorMap.Pressure.max    = 1.00
Visu.colorMap.CompactionPressure.scale  = 2e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.0
Visu.colorMap.FluidPressure.scale  = 10e6/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00




Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0

Visu.colorMap.SIIOvYield.max = 1.05


###          Write the Input file
### =====================================
Input.writeInputFile(Setup)


