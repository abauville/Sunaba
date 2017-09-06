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
from math import pi, sqrt, tan, sin, cos, exp
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
MPa     = 1e6
GPa     = 1e9

deg     = pi/180




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
Setup.Description = "Angle of shear bands benchmark, based on Kaus, 2010 (doi:10.1016/j.tecto.2009.08.042)"



Numerics.phiMin = 1e-5
Numerics.phiMax = 0.9

Numerics.etaMin = 1e-5
Numerics.etaMax = 1e+5

##          Material properties
## =====================================
StickyAir   = Input.Material("StickyAir")
Matrix      = Input.Material("Sediments")
Inclusion   = Input.Material("Sediments")


Setup.MatProps = {"0":StickyAir,"1":Matrix,"2":Inclusion}

PhaseRef = Matrix#Inclusion
PhaseRef.isRef = True

Matrix.name     = "Sediment"
Inclusion.name  = "Inclusion"

Matrix.vDiff    = material.DiffusionCreep       ("Off")
Inclusion.vDiff = material.DiffusionCreep       ("Off")
StickyAir.vDiff = material.DiffusionCreep       (eta0=1E16)

Matrix.use_dtMaxwellLimit = True

Matrix.vDisl    = material.DislocationCreep     (eta0=1E24, n=1)
Inclusion.vDisl = material.DislocationCreep     (eta0=1E20, n=1)



Matrix.rho0     = 0.0*2700  * kg/(m**3)
Inclusion.rho0  = 0.0*2700  * kg/(m**3)

Matrix.cohesion     = 50    * MPa
Inclusion.cohesion  = 50    * MPa

Matrix.frictionAngle    = 30 * deg
Inclusion.frictionAngle = 30 * deg

Matrix.G                = 1.0 * GPa
Inclusion.G             = 1.0 * GPa
StickyAir.G             = Matrix.G
StickyAir.cohesion      = .1 * MPa

#StickyAir.cohesion = Matrix.cohesion



##              Grid
## =====================================
RFac = 1; # Resolution Factor
HFac = 1.0


H = HFac * 1 * km
HStickyAir = H/5.0
Grid.ymin =  0.0
Grid.ymax =  H + HStickyAir
Grid.nyC = 128
dy = (Grid.ymax-Grid.ymin)/(Grid.nyC+1)



r =4*dy# H/8.0         # inclusion radius
d = 2.0*r
theta = 33/180*pi # effective shear zone angle
W = r*cos(45/180*pi) + (H-r*sin(45/180*pi))/tan(theta) # takes into account that the shear zone starts at 45 degree on the inclusion perimeter
#HStickyAir = 0.0
W = W*1.5



Grid.xmin = 0.0
Grid.xmax = W
#Grid.xmin = -W
#Grid.xmax = 0.0


Grid.nxC = round(128*W/(H+HStickyAir))

dx = (Grid.xmax-Grid.xmin)/(Grid.nxC+1)




Grid.fixedBox = True

Physics.Pback = 100 * MPa


##              Numerics
## =====================================
Numerics.nTimeSteps = -100
BCStokes.backStrainRate = -1.0e-15
Numerics.CFL_fac_Stokes = 0.25
Numerics.CFL_fac_Darcy = 0.8
Numerics.CFL_fac_Thermal = 10.0
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 3
Numerics.maxNonLinearIter = 20

Numerics.absoluteTolerance = 1e-5


Numerics.dtMaxwellFac_EP_ov_E  = 0.5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = 0.0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = 0.5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress
Numerics.use_dtMaxwellLimit = False



Particles.nPCX = 4
Particles.nPCY = 4
Particles.noiseFactor = 0.0



##              Non Dim
## =====================================

#L = (Grid.xmax-Grid.xmin)/2.0
#L = (Grid.ymax-Grid.ymin)/2.0
#BCStokes.backStrainRate = - BCStokes.refValue / L

#Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid,Length=H/2.0)
#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)

#Char.length =  (Grid.xmax-Grid.xmin)/2
#Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
#RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
##Char.time   = 1.0/abs(BCStokes.backStrainRate)/1.0
#Char.time   = RefVisc/PhaseRef.G*0.01#(1.0/abs(BCStokes.backStrainRate)/1.0*1e-6)
#CharVisc = 1.0/(1.0/RefVisc + 1.0/(PhaseRef.G*Char.time))#RefVisc
#CharStress  = 2.0*CharVisc*1.0/Char.time#PhaseRef.rho0*abs(Physics.gy)*Char.length
#Char.mass   = CharStress*Char.time*Char.time*Char.length


Char.length =  (Grid.xmax-Grid.xmin)/2
Char.temperature = (BCThermal.TB + BCThermal.TT)/2.0
CharStress =    Matrix.cohesion*cos(Matrix.frictionAngle) + Physics.Pback *sin((Matrix.frictionAngle))
n = 50.0
DeltaSigma = CharStress/n;
G = Matrix.G
EII = abs(BCStokes.backStrainRate)
eta = Matrix.getRefVisc(0.0,Char.temperature,EII)
t = 0.0
Char.time = DeltaSigma / (2*G*EII * exp(-G/eta*t));
Char.mass   = CharStress*Char.time*Char.time*Char.length


#
#
#CharStress = Physics.Pback



Numerics.dtMin = Char.time
Numerics.dtMax = Numerics.dtMin

#Numerics.maxTime = 500000*yr

##                 BC
## =====================================
#BCStokes.SetupType = "Sandbox"
#BCStokes.Sandbox_NoSlipWall = True


##              Geometry
## =====================================
#W = Grid.xmax-Grid.xmin
#H = Grid.ymax-Grid.ymin


inclusion_w = d/2
inclusion_h = d


slope = tan(0*pi/180)

i = 0
MatrixPhase = 1
Geometry["%05d_line" % i] = Input.Geom_Line(MatrixPhase,slope,H,"y","<",Grid.xmin,Grid.xmax)
InclusionPhase = 2
i+=1
Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w,"y","<",Grid.xmin,Grid.xmin+inclusion_w)
#Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w*4,"y","<",Grid.xmin,Grid.xmin+inclusion_w)

#Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w,"y","<",Grid.xmin+W/3.0,Grid.xmin+W/3.0+inclusion_w)

#Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w,"y","<",Grid.xmin+2.0*W/3.0,Grid.xmin+2.0*W/3.0+2.0*inclusion_w)
#Geometry["%05d_sine" % i] = Input.Geom_Sine(InclusionPhase,inclusion_w,inclusion_w,0.0,W/8.0,"y","<",Grid.xmin,Grid.xmax)
#Geometry["%05d_line" % i] = Input.Geom_Line(InclusionPhase,0.0,inclusion_w,"y","<",Grid.xmax-inclusion_w,Grid.xmax)


#Geometry["%05d_circle" % i] = Input.Geom_Circle(InclusionPhase, Grid.xmin, Grid.ymin, inclusion_w)




CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


###              Output
### =====================================
Output.folder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/nx_%i_ny_%i_G_%.2e_C_%.2e_fric_%.2e_Pref_%.2e" % (Grid.nxC, Grid.nyC, Matrix.G, Matrix.cohesion, Matrix.frictionAngle*180/pi, Physics.Pback)
Output.folder = "/Users/abauville/Output_Paper_DynDecollement/DynStress_PureShear/Test"
Output.strainRate = True
Output.sigma_II = True
Output.sigma_xx = True
Output.sigma_xy = True
Output.khi = True
Output.P = True

Output.frequency = 1#timeFac



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
Visu.transparency = False

Visu.showGlyphs =  False
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 250.0
glyphSpacing = (Grid.ymax-Grid.ymin)/24 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 1.0 * Visu.height
Visu.width = 1* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = 0.0


print("\n"*5)
CharExtra = Input.CharExtra(Char)
RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))


print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*H/2.0

Visu.colorMap.Stress.scale  = 100.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.0


Visu.colorMap.Pressure.scale  = Physics.Pback/CharExtra.stress
Visu.colorMap.Pressure.center = 1.0
Visu.colorMap.Pressure.max    = 2.0


Visu.colorMap.VelocityDiv.scale = 1e-1

Visu.colorMap.Khi.max = 5.0
Visu.colorMap.Khib.max = 5.0




###          Write the Input file
### =====================================
os.system("mkdir " + Output.folder)
os.system("cd " + Output.folder + "\n rm -r *")
if (Visu.writeImages):
    os.system("mkdir " + Visu.outputFolder)
Input.writeInputFile(Setup)


