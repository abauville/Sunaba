#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:24:44 2016

@author: abauville
"""

# Input Test for Stokes FD
import sys
import os
sys.path.insert(0, '../../../src/UserInput')
#import json
#from InputDef import *
import InputDef as Input
import MaterialsDef as material
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *
from math import pi, sqrt, tan, sin, cos, log10, log2
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

Pa      = kg/m/s/s
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

ProductionMode = False
Numerics.phiCrit = 1e-3
Numerics.phiMin = 1e-4
Numerics.phiMax = 0.9

Numerics.etaMin = 1e-8
Numerics.etaMax = 1e8

##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Sediment    = Input.Material("Sediments")
Basement    = Input.Material("Sediments")
WeakLayer   = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Basement, "3":WeakLayer}




PhaseRef = Sediment
#PhaseRef = StickyAir
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"
Basement.name = "Basement"
WeakLayer.name = "WeakLayer"

Sediment.vDiff = material.DiffusionCreep       ("Off")
Basement.vDiff = material.DiffusionCreep       ("Off")
WeakLayer.vDiff = material.DiffusionCreep       ("Off")
#Basement.vDiff = material.DiffusionCreep       (eta0 = 1e23)

#Sediment.vDisl = material.DislocationCreep     (eta0=1E90, n=10)
#Basement.vDisl = material.DislocationCreep     (eta0=1E150, n=10)



#StickyAir.rho0 = 1.0
StickyAir.rho0 = 0000.00


StickyAir.phiIni = Numerics.phiMin
Sediment.phiIni = Numerics.phiMin
WeakLayer.phiIni = 0.6
Basement.phiIni = Numerics.phiMin

StickyAir.perm0 = 1e-8
WeakLayer.perm0 = 1e-8
Sediment.perm0 = 1e-8
Basement.perm0 = 1e-12



Sediment.G  = 1e9
WeakLayer.G = 1e9

Basement.G  = Sediment.G*100.0
StickyAir.G = Sediment.G*1.0
StickyAir.cohesion = 1e6/1.0#1.0*Sediment.cohesion

Sediment.use_dtMaxwellLimit = True


## Main parameters for this setup
## =====================================

Sediment.frictionAngle  = 30/180*pi
WeakLayer.frictionAngle = 30/180*pi
Basement.frictionAngle  = Sediment.frictionAngle



WeakLayer.cohesion = 10e6
Sediment.cohesion =  10e6
Basement.cohesion = 50*1e6

Numerics.deltaSigmaMin = 10.0 * MPa

HFac = 1.0


LWRatio = 1.5
Hsed = HFac*1.0e3


Grid.xmin = -1.5*Hsed*LWRatio
Grid.xmax = 0.0e3
Grid.ymin = 0.0e3
Grid.ymax = 1.5*Hsed
if ProductionMode:
    Grid.nxC = round(1/1*((64+64+128)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(1/1*((64+64+128)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)
else:
    Grid.nxC = round(1/1*((64)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(1/1*((32)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.fixedBox = True

print("Grid.nxC = %i, Grid.nyC = %i" % (Grid.nxC, Grid.nyC))



VatBound = - 10 * cm/yr
dx = (Grid.xmax-Grid.xmin)/Grid.nxC
dy = (Grid.ymax-Grid.ymin)/Grid.nyC
BCStokes.backStrainRate = VatBound / (Grid.xmax-Grid.xmin)

Plitho = Sediment.rho0 * abs(Physics.gy) * 1.0*Hsed
Sigma_y = Sediment.cohesion*cos(Sediment.frictionAngle) + sin(Sediment.frictionAngle)*1.0*Plitho
print("RefViscBrittle = %.2e Pa.s" % (Sigma_y/abs(BCStokes.backStrainRate)))
print("backStrainRate = %.2e, Sigma_y = %.2e MPa" % (BCStokes.backStrainRate, Sigma_y/1e6))


RefVisc =  (Sigma_y/abs(BCStokes.backStrainRate))


RefVisc *= 1
StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1000)
Sediment.vDisl = material.DislocationCreep     (eta0=RefVisc*1, n=1)
WeakLayer.vDisl = material.DislocationCreep    (eta0=RefVisc*1, n=1)
Basement.vDisl = material.DislocationCreep     (eta0=RefVisc*10000, n=1)



BoxTilt = 0 * pi/180
slope = -BoxTilt #tan(0*pi/180)

Physics.gx = -9.81*sin(BoxTilt);
Physics.gy = -9.81*cos(BoxTilt);




##              Grid

#Grid.xmin = -1000.0e3
#Grid.xmax =  1000e3
#Grid.ymin = -380e3
#Grid.ymax =  20.0e3



##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin

Hbase = HFac*0.1e3

Wseamount = .15e3*HFac
xseamount = Grid.xmin + 1e3

i = 0
SedPhase = 1
BasementPhase = 2
WeakPhase = 3

Lweak = Grid.xmax-Grid.xmin
Hweak = .24e3*HFac
ThickWeak = .05e3*HFac



Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W,"y","<",Grid.xmin,Grid.xmax)

## Weak Layer
#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(WeakPhase,slope,Hweak - slope*W,"y","<",Grid.xmin,Grid.xmin+Lweak)
#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hweak - ThickWeak - slope*W,"y","<",Grid.xmin,Grid.xmin+Lweak)
#
#

#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,0.0,Hbase,"y","<",Grid.xmin,Grid.xmax)



HSFac = 8
#BCStokes.Sandbox_TopSeg00 = 0.395e3*HFac
BCStokes.Sandbox_TopSeg00 = Hbase + 0*Hbase + dy + 0*HSFac*dy
BCStokes.Sandbox_TopSeg01 = BCStokes.Sandbox_TopSeg00+HSFac*dy#0.405e3*HFac

#
#i+=1
#BackStopSlope = BoxTilt#tan(-10*pi/180)
#Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,BackStopSlope,Grid.xmax-Hbase,"x",">",BCStokes.Sandbox_TopSeg01,Grid.ymax-3*Hbase)





##              Numerics
## =====================================
Numerics.nTimeSteps = 50000
Numerics.CFL_fac_Stokes = .5
Numerics.CFL_fac_Darcy = 1000.0
Numerics.CFL_fac_Thermal = 10000.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 2
if ProductionMode:
    Numerics.maxNonLinearIter = 15
else:
    Numerics.maxNonLinearIter = 20
Numerics.dtAlphaCorr = .3
Numerics.absoluteTolerance = 1e-6


Numerics.dtMaxwellFac_EP_ov_E  = .5   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = True


Numerics.maxTime = 1e6*yr

timeFac = 4
Numerics.dtMin = 1.0*s #50/4*yr
Numerics.dtMax = 1e4 * yr#Numerics.dtMin




if (ProductionMode):
    Particles.nPCX = 4
    Particles.nPCY = 4
    Particles.noiseFactor = 0.75
#    Particles.minPartPerCellFactor = 0.5
else:
    Particles.nPCX = 4
    Particles.nPCY = 32
    Particles.noiseFactor = 0.5
#    Particles.minPartPerCellFactor = 0.5
    


###              Output
### =====================================
#Output.folder = "/Users/abauville/Output_Paper_DynDecollement/DynStress/nx_%i_ny_%i_G_%.2e_C_%.2e_fric_%.2e_Hsed_%.2e" % (Grid.nxC, Grid.nyC, Sediment.G, Sediment.cohesion, Sediment.frictionAngle*180/pi, Hsed)
#Output.strainRate = True
#Output.sigma_II = True
#Output.khi = True
#Output.P = True
#
#Output.frequency = timeFac





##                 BC
## =====================================
BCStokes.SetupType = "Sandbox"

BCStokes.Sandbox_NoSlipWall = True


##                 IC
## =====================================
ICThermal.age = 100*Myr

ICDarcy.background = 0.0#Numerics.phiMin
ICDarcy.Amp = 0.0
ICDarcy.xc = Grid.xmin+(Grid.xmax-Grid.xmin)/2
ICDarcy.yc = Grid.ymin+(Grid.ymax-Grid.ymin)/2
ICDarcy.wx = (Grid.xmax-Grid.xmin)/24.0
ICDarcy.wy = (Grid.xmax-Grid.xmin)/24.0

##              Non Dim
## =====================================
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
BaseVisc = Basement.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

L = (Grid.ymax-Grid.ymin)/2.0

Char.length = Hsed/1.0

Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
Char.time   = 100*yr
CharVisc = 1.0/(1.0/SedVisc + 1.0/(Sediment.G*Char.time))#RefVisc
CharStress  = PhaseRef.rho0*abs(Physics.gy)*Char.length
Char.mass   = CharStress*Char.time*Char.time*Char.length


##              Info
## =====================================
print("\n"*5)
CharExtra = Input.CharExtra(Char)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)
print("BaseVisc = %.2e" %  BaseVisc)

print("Lc = " + str(  (Sediment.cohesion*cos(Sediment.frictionAngle)) / (Sediment.rho0 * abs(Physics.gy) * sin(Sediment.frictionAngle))  ))




##            Visualization
## =====================================

Particles.passiveGeom = "Grid_w_Layers"

Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/32
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.shaderFolder = "../Shaders/Sandbox_w_Layers" # Relative path from the running folder (of StokesFD)

Visu.type = "StrainRate"
#if ProductionMode:
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/StokesFD_Outputs/Test_Sandbox_ObliqueBackStop"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output"
Visu.transparency = False

Visu.glyphMeshType = "TensorCross"
Visu.glyphType = "DeviatoricStressTensor"
#Visu.showGlyphs = True
#Visu.glyphScale = 8.0/(abs(VatBound)/(Char.length/Char.time))
Visu.glyphScale = 0.2
glyphSpacing = (Grid.ymax-Grid.ymin)/32 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1.0 * Visu.height
Visu.width = 1.0 * Visu.width

Visu.filter = "Nearest"

Visu.shiftFacY = -0.5
Visu.shiftFacZ = 0.1

print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = .5*Plitho/CharExtra.stress
Visu.colorMap.Stress.center = 0.0
Visu.colorMap.Stress.max    = 2.00
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Porosity.max = 1.0

Visu.colorMap.Pressure.scale  = .5*Plitho/CharExtra.stress
Visu.colorMap.Pressure.center = 2.0
Visu.colorMap.Pressure.max    = 4.00


Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0

Visu.colorMap.Velocity.log10on = True
Visu.colorMap.Velocity.scale = (10.0*cm/yr) / (Char.length/Char.time)#abs(VatBound) / (Char.length/Char.time)
Visu.colorMap.Velocity.center = 0
Visu.colorMap.Velocity.max = 2.0

Visu.colorMap.Vorticity.max = 0.0005/yr /  (1.0/Char.time) # in rad/yr


Visu.colorMap.POvPlitho.log10on = True
Visu.colorMap.POvPlitho.center = 0.0
Visu.colorMap.POvPlitho.max = log10(2.0)

Visu.closeAtTheEndOfSimulation = False

###          Write the Input file
### =====================================
Input.writeInputFile(Setup)

os.system("mkdir " + Visu.outputFolder)
os.system("mkdir " + Output.folder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")
