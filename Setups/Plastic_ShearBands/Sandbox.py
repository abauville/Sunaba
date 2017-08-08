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

Numerics.etaMin = 1e-5
Numerics.etaMax = 1e5

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
StickyAir.G = Sediment.G*100.0
StickyAir.cohesion = 1e6/1.0#1.0*Sediment.cohesion

Sediment.use_dtMaxwellLimit = True


## Main parameters for this setup
## =====================================

Sediment.frictionAngle  = 30/180*pi
WeakLayer.frictionAngle = 30/180*pi
Basement.frictionAngle  = Sediment.frictionAngle
slope = tan(0*pi/180)


WeakLayer.cohesion = 10e6
Sediment.cohesion =  10e6
Basement.cohesion = 50*1e6

HFac = 1.0


LWRatio = 1.0
Hsed = HFac*1.5e3


Grid.xmin = -3.0*Hsed*LWRatio
Grid.xmax = 0.0e3
Grid.ymin = 0.0e3
Grid.ymax = 3.0*Hsed
if ProductionMode:
    Grid.nxC = round(1/1*((64+64+128)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(1/1*((64+64+128)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)
else:
    Grid.nxC = round(1/1*((64+16)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(1/1*((64+16)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

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

#Numerics.etaMax = 1e10;



#WeakLayer.cohesion = 1e30
#WeakLayer.cohesion = 1e30
#Sediment.cohesion = 1e30
#Basement.cohesion = 1e30


BoxTilt = -00 * pi/180
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

Hbase = HFac*0.15e3

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
i+=1
Geometry["%05d_line" % i] = Input.Geom_Line(BasementPhase,slope,Hbase - slope*W,"y","<",Grid.xmin,Grid.xmax)
i+=1
#Geometry["%05d_sine" % i] = Input.Geom_Sine(BasementPhase,Hbase - slope*W,3*Hbase,0,Wseamount*2,"y","<",xseamount-Wseamount/2,xseamount+Wseamount/2)
Geometry["%05d_sine" % i] = Input.Geom_Sine(BasementPhase,Hbase - slope*W,0*0.25*Hbase,pi/16,Wseamount*2/3,"y","<",Grid.xmin,Grid.xmax)
i+=1
Geometry["%05d_sine" % i] = Input.Geom_Sine(BasementPhase,Hbase - slope*W,0*0.25*Hbase,pi+pi/16,Wseamount*2/3,"y","<",Grid.xmin,Grid.xmax)


HSFac = 2
#BCStokes.Sandbox_TopSeg00 = 0.395e3*HFac
BCStokes.Sandbox_TopSeg00 = Hbase + HSFac*dy
BCStokes.Sandbox_TopSeg01 = BCStokes.Sandbox_TopSeg00+HSFac*dy#0.405e3*HFac

##              Numerics
## =====================================
Numerics.nTimeSteps = 1000
Numerics.CFL_fac_Stokes = .25
Numerics.CFL_fac_Darcy = 1000.0
Numerics.CFL_fac_Thermal = 10000.0
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 2
if ProductionMode:
    Numerics.maxNonLinearIter = 15
else:
    Numerics.maxNonLinearIter = 2
Numerics.dtAlphaCorr = .3
Numerics.absoluteTolerance = 1e-6


Numerics.dtMaxwellFac_EP_ov_E  = .5   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = True

Numerics.maxTime = (Grid.xmax-Grid.xmin)/abs(VatBound)

Numerics.dtMin = 100*yr
Numerics.dtMax = Numerics.dtMin


#Numerics.dtVep = 1.0*Numerics.CFL_fac_Stokes*dx/abs(VatBound) 


if (ProductionMode):
    Particles.nPCX = 4
    Particles.nPCY = 4
    Particles.noiseFactor = 0.75
#    Particles.minPartPerCellFactor = 0.5
else:
    Particles.nPCX = 4
    Particles.nPCY = 4
    Particles.noiseFactor = 0.0
#    Particles.minPartPerCellFactor = 0.5
    


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
BCStokes.SetupType = "Sandbox"

BCStokes.Sandbox_NoSlipWall = True

#BCStokes.refValue       = 1.0 * cm/yr / 1.0


#BCThermal.TB = 30.0 + 273.0
#BCThermal.TT = 0.0  + 273.0


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
SedVisc = Sediment.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
BaseVisc = Basement.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

#L = (Grid.xmax-Grid.xmin)/2.0
L = (Grid.ymax-Grid.ymin)/2.0
#BCStokes.backStrainRate = - BCStokes.refValue / L

#Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid,Hsed/8.0)
Char.length = Hsed/8.0

Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
CharVisc = 1.0/(1.0/SedVisc + 1.0/(Sediment.G*Numerics.dtMin))#RefVisc
  
CharStress  = PhaseRef.rho0*abs(Physics.gy)*Char.length


Char.time   = CharVisc/CharStress
Char.mass   = CharStress*Char.time*Char.time*Char.length


##              Info
## =====================================
print("\n"*5)
CharExtra = Input.CharExtra(Char)
#RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


#StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/1000.0)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)
print("BaseVisc = %.2e" %  BaseVisc)



##            Visualization
## =====================================

Particles.passiveGeom = "Grid_w_Layers"

Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/48
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.shaderFolder = "../Shaders/Sandbox_w_Layers" # Relative path from the running folder (of StokesFD)
#Visu.shaderFolder = "../Shaders/CornerFlow"

Visu.type = "StrainRate"
#if ProductionMode:
Visu.writeImages = False
    
#Visu.outputFolder = "/Users/abauville/StokesFD_Output/EffectiveStrainRateFormulationTest"
Visu.outputFolder = "/Users/abauville/GoogleDrive/StokesFD_Output/Test_Sandbox_New"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_SandboxNew/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Sandbox_Outputs/PfHydro_dt99_01_G5e8/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Seismic_Sandbox_Outputs/nx%i_ny%i_G%.e_D%.f_C%.1e_fric%.f_MethodAv_HSFac%i_dtMaxwell_08_02_ManyIter/" % (Grid.nxC, Grid.nyC, Sediment.G, HFac, Sediment.cohesion, Sediment.frictionAngle*180/pi, HSFac)
#Visu.outputFolder = "/Users/abauville/StokesFD_Outputs/Seismic_Sandbox_Outputs/nx%i_ny%i_G%.e_D%.f_C%.1e_fric%.f_HSFac%i_dtMaxwell_05_05_New/" % (Grid.nxC, Grid.nyC, Sediment.G, HFac, Sediment.cohesion, Sediment.frictionAngle*180/pi, HSFac)

Visu.transparency = False

#Visu.showGlyphs = True
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 8.0/(abs(VatBound)/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1.0 * Visu.height
Visu.width = 1.0 * Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.5
Visu.shiftFacZ = 0.1





print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = .5*Plitho/CharExtra.stress
Visu.colorMap.Stress.center = 1.0
Visu.colorMap.Stress.max    = 2.00
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

Visu.colorMap.Pressure.scale  = .5*Plitho/CharExtra.stress
Visu.colorMap.Pressure.center = 2.0
Visu.colorMap.Pressure.max    = 4.00
Visu.colorMap.CompactionPressure.scale  = 5e6/CharExtra.stress
Visu.colorMap.CompactionPressure.center = 0.0
Visu.colorMap.CompactionPressure.max    = 1.0
Visu.colorMap.FluidPressure.scale  = 1*Plitho/CharExtra.stress
Visu.colorMap.FluidPressure.center = 0.0
Visu.colorMap.FluidPressure.max    = 1.00




Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0



Visu.colorMap.Velocity.log10on = True
Visu.colorMap.Velocity.scale = (10.0*cm/yr) / (Char.length/Char.time)#abs(VatBound) / (Char.length/Char.time)
Visu.colorMap.Velocity.center = 0
Visu.colorMap.Velocity.max = 2.0

Visu.colorMap.Vorticity.max = 0.0005/yr /  (1.0/Char.time) # in rad/yr


Visu.colorMap.POvPlitho.log10on = True
Visu.colorMap.POvPlitho.center = 0.0
Visu.colorMap.POvPlitho.max = log10(2.0)


##              Some info
## ======================================
print("Lc = " + str(  (Sediment.cohesion*cos(Sediment.frictionAngle)) / (Sediment.rho0 * abs(Physics.gy) * sin(Sediment.frictionAngle))  ))
#CompactionLength = sqrt(4.0/3.0*Sediment.perm0/Physics.eta_f * (Physics->eta[iCell]/phi));



###          Write the Input file
### =====================================
Input.writeInputFile(Setup)

os.system("mkdir " + Visu.outputFolder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")
