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
from math import pi, sqrt, tan, sin, cos, log10, log2, log
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
Setup.Description = "phib 28, cohesionWeak = 0.5, fricWeak = 0.75"

ProductionMode = True
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
Backstop   = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Sediment,"2":Basement, "3":Backstop}

Numerics.stickyAirSwitchPhaseTo = 1
Numerics.stickyAirSwitchPassiveTo = 4


PhaseRef = Sediment
#PhaseRef = StickyAir
PhaseRef.isRef = True

StickyAir.name = "StickyAir"
Sediment.name = "Sediment"
Basement.name = "Basement"
Backstop.name = "Backstop"

Sediment.vDiff = material.DiffusionCreep       ("Off")
Basement.vDiff = material.DiffusionCreep       ("Off")
Backstop.vDiff = material.DiffusionCreep       ("Off")

StickyAir.rho0 = 0000.00


StickyAir.phiIni = Numerics.phiMin
Sediment.phiIni = Numerics.phiMin
Backstop.phiIni = 0.6
Basement.phiIni = Numerics.phiMin

StickyAir.perm0 = 1e-8
Backstop.perm0 = 1e-8
Sediment.perm0 = 1e-8
Basement.perm0 = 1e-12

Sediment.G  = 1e9
Backstop.G = 1e9

Basement.G  = Sediment.G*10.0
StickyAir.G = Sediment.G/2.0


Sediment.use_dtMaxwellLimit = False

Numerics.invariantComputationType = 0

##              Numerics
## =====================================
Numerics.nTimeSteps = 10000000
Numerics.CFL_fac_Stokes = .25
Numerics.CFL_fac_Darcy = 1000.0
Numerics.CFL_fac_Thermal = 10000.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 5
if ProductionMode:
    Numerics.maxNonLinearIter = 30
else:
    Numerics.maxNonLinearIter = 30
#    Numerics.maxNonLinearIter = 10
    Numerics.dtAlphaCorr = .3
Numerics.absoluteTolerance = 1e-6
Numerics.relativeTolerance  = 1e-3


Numerics.dtMaxwellFac_EP_ov_E  = .5   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = .0   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = .5   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = True



Numerics.dt_stressFac = 0.5 # between 0 and 1; dt = Fac*time_needed_to_reach_yield # i.e. see RefTime in this file
Numerics.dt_plasticFac = 0.5 # between 0 and 1; 0 = EP/E limit; 1 = VP/EP limit


Numerics.stressSubGridDiffFac = 1.0



#Numerics.dtMin = 1e-2   *yr #0.1*Char.time #50/4*yr
#Numerics.dtMax = 1e3   *yr#50.0*Char.time#Numerics.dtMin


if (ProductionMode):
    Particles.nPCX = 4
    Particles.nPCY = 4
    Particles.noiseFactor = 0.0
#    Particles.minPartPerCellFactor = 0.5
else:
    Particles.nPCX = 4
    Particles.nPCY = 4
    Particles.noiseFactor = 0.00
#    Particles.minPartPerCellFactor = 0.5
    

Numerics.yieldComputationType = 1


## Main parameters for this setup
## =====================================

Sediment.frictionAngle  = 30/180*pi
Backstop.frictionAngle = 30/180*pi
Basement.frictionAngle  = Sediment.frictionAngle



Backstop.cohesion = 1.0e6
Sediment.cohesion =  5.0*1.0e6# * 20.0
Basement.cohesion = 50*1e6
StickyAir.cohesion = 1.0*Sediment.cohesion

HFac        = 1.0
Hsed        = HFac*1.0e3


if (ProductionMode):
    ResFac      = 2.0
    LWRatio     = 3.5
else:
    ResFac      = 0.75
    LWRatio     = 2.5

timeFac = 5
Numerics.maxTime = 40*12800*yr * HFac
#Numerics.dtMin = 2**timeFac   *yr * HFac #0.1*Char.time #50/4*yr
#Numerics.dtMax = 2**timeFac   *yr * HFac#50.0*Char.time#Numerics.dtMin


Numerics.dtMin = 2**1   *yr * HFac #0.1*Char.time #50/4*yr
Numerics.dtMax = 2**10   *yr * HFac#50.0*Char.time#Numerics.dtMin



if (ProductionMode):
    Grid.xmin = -6.5*Hsed*LWRatio
    Grid.ymax = 6.5*Hsed
else:
    Grid.xmin = -3.5*Hsed*LWRatio
    Grid.ymax = 3.5*Hsed

Grid.xmax = 0.0e3
Grid.ymin = 0.0e3


if ProductionMode:
    Grid.nxC = round(2**ResFac*((48)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(2**ResFac*((48)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)
else:
    Grid.nxC = round(2**ResFac*((48)*LWRatio)) #round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
    Grid.nyC = round(2**ResFac*((48)))#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

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


RefVisc =  10.0*(Sigma_y/abs(BCStokes.backStrainRate))


RefVisc *= 1
StickyAir.vDiff = material.DiffusionCreep(eta0=RefVisc/100000)
Sediment.vDisl = material.DislocationCreep     (eta0=RefVisc*100, n=1)
Backstop.vDisl = material.DislocationCreep    (eta0=RefVisc*1, n=1)
Basement.vDisl = material.DislocationCreep     (eta0=RefVisc*100, n=1)



BoxTilt = 0 * pi/180
slope = -BoxTilt #tan(0*pi/180)

Physics.gx = -9.81*sin(BoxTilt);
Physics.gy = -9.81*cos(BoxTilt);








##              Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = Grid.ymax-Grid.ymin

Hbase = 0*HFac*0.2e3

Wseamount = .15e3*HFac
xseamount = Grid.xmin + 1e3

i = 0
AirPhase = 0
SedPhase = 1
BasementPhase = 2
BackPhase = 3

Lweak = Grid.xmax-Grid.xmin
Hweak = .24e3*HFac
ThickWeak = .05e3*HFac



Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W,"y","<",Grid.xmin,Grid.xmax)
#i+=1
#Geometry["%05d_line" % i] = Input.Geom_Line(BackPhase,slope,Hsed - slope*W,"y","<",Grid.xmin+W/2,Grid.xmax)



#i+=1
#slope = 15.0*pi/180.0
#Geometry["%05d_line" % i] = Input.Geom_Line(SedPhase,slope,Hsed - slope*W/2.0,"y","<",Grid.xmin,Grid.xmax)


HSFac = 1
BCStokes.Sandbox_TopSeg00 = Hbase + 0*Hbase + 0*dy + 0*HSFac*dy
BCStokes.Sandbox_TopSeg01 = BCStokes.Sandbox_TopSeg00+HSFac*dy#0.405e3*HFac





###              Output
### =====================================
#baseFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/CorotationalNewInvType1/FixedDt_Method%i/" % Numerics.yieldComputationType
#baseFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/Test3b/"
baseFolder = "/Users/abauville/Output/Paper_AccretionVsSubduction/HFac%.0f_G%.0e_CWeak/" % (HFac,Sediment.G)
#baseFolder = "/Users/abauville/Output/EGU2018_PosterDecollement/StrucStyle/Test/"
##baseFolder = "/Users/abauville/Output/EGU2018_PosterFail/dxdtSensitivity3/AdaptativeDt_UpperConvected_Method0/"


Output.folder = (baseFolder + "Output/dxFac%i_dtFac%i" % (ResFac, timeFac) )
Output.strainRate = True
Output.strain     = True
Output.sigma_II = True
Output.khi = True
Output.P = True
Output.phase = True
Output.sigma_xx = True
Output.sigma_xy = True
Output.Vx = True
Output.Vy = True

Output.frequency = 4#round(128*yr/Numerics.dtMin)
#Output.timeFrequency = 50*yr




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
#
#L = (Grid.ymax-Grid.ymin)/2.0
#
#Char.length = Hsed/1.0
#
#Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
#Char.time   = 100*yr
#CharVisc = 1.0/(1.0/SedVisc + 1.0/(Sediment.G*Char.time))#RefVisc
#CharStress  = PhaseRef.rho0*abs(Physics.gy)*Char.length
#Char.mass   = CharStress*Char.time*Char.time*Char.length


##              Non Dim
## =====================================
Char.length =  (Grid.xmax-Grid.xmin)/2
Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
#    CharStress =    Matrix.cohesion*cos(Matrix.frictionAngle) + Physics.Pback *sin((Matrix.frictionAngle))

#if ProductionMode:
#    a_f = 100.0
#else:    
#    a_f = 100.0
#    




#timeFac = -1
#DeltaSigma = CharStress*dt_stressFac ;
G = Sediment.G
EII = abs(BCStokes.backStrainRate)
eta = Sediment.getRefVisc(0.0,Char.temperature,EII)

S3 = Setup.Physics.Pback
C = Sediment.cohesion
phi = Sediment.frictionAngle
S1 = 1.0/(1.0-sin(phi)) * (  2*C*cos(phi) + S3*(1+sin(phi))  )
Sy_back = (S1-S3)/2.0
P_Lim = (S1+S3)/2.0
#    P = Setup.Physics.Pback
#    Sy_back = C*cos(phi) + P*sin(phi)
RefTime  = eta/G * log(2.0*eta*EII / (2.0*eta*EII - Sy_back )); # time at which stress has built up to the 
#Char.time = timeFac*RefTime*Numerics.dt_stressFac
Char.time = Numerics.dtMin




CharVisc = 1.0/(1.0/eta+1.0/(G*Char.time))
CharStress = CharVisc/Char.time

Char.mass   = CharStress*Char.time*Char.time*Char.length





















##              Info
## =====================================
print("\n"*5)
CharExtra = Input.CharExtra(Char)

StickyAirVisc = StickyAir.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))

CharVisc2 = 1.0/(1.0/StickyAirVisc+1.0/(StickyAir.G*Char.time))


print("RefVisc = %.2e" % RefVisc)
print("Sediment Visc = %.2e" % SedVisc)
print("StickyAirVisc = %.2e" % StickyAirVisc)
print("BaseVisc = %.2e" %  BaseVisc)

print("Lc = " + str(  (Sediment.cohesion*cos(Sediment.frictionAngle)) / (Sediment.rho0 * abs(Physics.gy) * sin(Sediment.frictionAngle))  ))




##            Visualization
## =====================================

Particles.passiveGeom = "Grid_w_Layers"

Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/10 / (Grid.ymax/Hsed)
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.shaderFolder = "../Shaders/Sandbox_w_Layers_Backstop" # Relative path from the running folder (of StokesFD)

#Visu.type = "StrainRate"
Visu.type = "StrainRate"
#if ProductionMode:
Visu.renderFrequency = 4#round(128*yr/Numerics.dtMin)
#Visu.renderTimeFrequency = 128*yr
Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/StokesFD_Output/Test_NewRotation"
#Visu.outputFolder = ("/Users/abauville/Output/Sandbox_NumericalConvergenceTest_NewRHS/dt_%.0fyr/ResFac_%.1f" % (Numerics.dtMin/yr, ResFac) )
Visu.outputFolder = (baseFolder + "Visu/dxFac%i_dtFac%i" % (ResFac, timeFac) )
Visu.transparency = True

Visu.glyphMeshType = "TensorCross"
Visu.glyphType = "DeviatoricStressTensor"
#Visu.showGlyphs = True
#Visu.glyphScale = 8.0/(abs(VatBound)/(Char.length/Char.time))
Visu.glyphScale = 0.2
#glyphSpacing = (Grid.ymax-Grid.ymin)/64 #50 * km
#Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
#Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1.0 * Visu.height
Visu.width = 1.5 * Visu.width


Visu.filter = "Nearest"

Visu.shiftFacY = -0.5
Visu.shiftFacZ = -0.1

print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0

Visu.colorMap.Stress.scale  = 0.25*Plitho/CharExtra.stress
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

Visu.colorMap.Khi.max = 10.0

Visu.colorMap.Velocity.log10on = True
Visu.colorMap.Velocity.scale = (10.0*cm/yr) / (Char.length/Char.time)#abs(VatBound) / (Char.length/Char.time)
Visu.colorMap.Velocity.center = 0
Visu.colorMap.Velocity.max = 2.0

Visu.colorMap.Vorticity.max = 0.0005/yr /  (1.0/Char.time) # in rad/yr





Visu.colorMap.POvPlitho.log10on = True
Visu.colorMap.POvPlitho.center = 0.0
Visu.colorMap.POvPlitho.max = log10(2.0)

#Visu.closeAtTheEndOfSimulation = False

Visu.colorMap.VxRes.scale = 1e-6
Visu.colorMap.VyRes.scale = 1e-6
Visu.colorMap.PRes.scale = 1e-6

###          Write the Input file
### =====================================
Input.writeInputFile(Setup)

#
os.system("mkdir " + Visu.outputFolder)
os.system("mkdir " + Output.folder)
os.system("/Users/abauville/JAMSTEC/StokesFD/Debug/StokesFD ./input.json")
