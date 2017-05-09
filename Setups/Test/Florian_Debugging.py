# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json

import InputDef as Input
import MaterialsDef as material
from math import pi
#from GeometryGraphical import *

# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import * 

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

## Description
## =====================================
Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"


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



##       Modify Material properties
## =====================================
Phase0 = Input.Material("StickyAir")
Phase1 = Input.Material("Sediments")
Phase2 = Input.Material("Sediments")

PhaseRef = Input.Material("Sediments")



PhaseRef.name = "Reference"
PhaseRef.vDisl = material.DislocationCreep    (eta0=1e20, n=1)
PhaseRef.rho0 = 2500
PhaseRef.cohesion = 10*1e6
PhaseRef.frictionAngle = 30/180*pi
PhaseRef.G    = 1e8
PhaseRef.isRef = True







Phase0.name = "StickyAir"
Phase0.vDisl = material.DislocationCreep    (eta0=1e16, n=1)
Phase0.rho0 = 10
Phase0.cohesion = .1e6
Phase0.frictionAngle = 30/180*pi
Phase0.G    = 1e11

Phase1.name = "Sediments"
Phase1.vDisl = material.DislocationCreep    (eta0=1E25, n=1)
Phase1.rho0 = 2500
Phase1.cohesion = 1e6
Phase1.frictionAngle = 30/180*pi
Phase1.G    = 1e8

Phase2.name = "Detachment"
Phase2.vDisl = material.DislocationCreep    (eta0=1E20, n=1)
Phase2.rho0 = 2500
Phase2.cohesion = 1e6
Phase2.frictionAngle = 1/180*pi
Phase2.G    = 1e8


Setup.MatProps = {'0': Phase0,'1': Phase1,'2': Phase2}



#BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1e-13
Numerics.CFL_fac_Stokes = .5
Numerics.nLineSearch = 4
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 10

Numerics.absoluteTolerance = 1e-20

Numerics.use_dtMaxwellLimit = False
#Numerics.dtVep = 100000.0*Numerics.CFL_fac_Stokes*BCStokes.backStrainRate


Grid.nyC = 128
Grid.nxC = 256

Grid.xmin = -50.0e3
Grid.xmax =  50.0e3
Grid.ymax =  10.0e3
Grid.ymin = 0

#BCStokes.SetupType = "PureShear"
#BCStokes.SetupType = "Sandbox"
#BCThermal.SetupType = "Sandbox"

Particles.nPCX = 4
Particles.nPCY = 4

#Physics.gy = 0.
#Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)
Char.set_based_on_lithostatic_pressure(PhaseRef,BCStokes,BCThermal,Physics,Grid)


Numerics.dtMaxwellFac_EP_ov_E  = 0.5;   # lowest,       ElastoPlasticVisc   /   G
Numerics.dtMaxwellFac_VP_ov_E  = 0.0;   # intermediate, ViscoPlasticVisc    /   G
Numerics.dtMaxwellFac_VP_ov_EP = 0.5;   # highest,      ViscoPlasticVisc    /   ElastoPlasticStress

Numerics.use_dtMaxwellLimit = True


##            Define Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = 0.6*(Grid.ymax-Grid.ymin)

DetHL = 0.25*H
DetHR = 0.15*H

InterH = 0.055*H
InterY = Grid.ymin+0.6*H-InterH/2


i = 0
phase = 1
Geometry["%05d_line" % i] = (Input.Geom_Line(phase,0.0,H,"y","<",Grid.xmin,Grid.xmax))

i+=1
phase = 2
Geometry["%05d_rect" % i] = (Input.Geom_Rect(phase,Grid.xmin,Grid.ymin,W/2,DetHL))

i+=1
phase = 2
Geometry["%05d_rect" % i] = (Input.Geom_Rect(phase,Grid.xmin+W/2-.0001,Grid.ymin,W/2,DetHR))

i+=1
phase = 2
Geometry["%05d_rect" % i] = (Input.Geom_Rect(phase,Grid.xmin,InterY,W,InterH))



##
##for key in Geometry:
##    Geometry[key].plot()
##
##plt.axis([Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax])
##plt.show()









Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.0*(Grid.xmax-Grid.xmin)/Char.length/Grid.nxC



Particles.noiseFactor = 0.95




##            Visualization
## =====================================
Visu.shaderFolder = "../Shaders/Sandbox" # Relative path from the running folder (of StokesFD)

Particles.passiveDy = (Grid.ymax-Grid.ymin)*1/16
Particles.passiveDx = Particles.passiveDy

Visu.showParticles = True
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC


Visu.type = "Vorticity"
#Visu.writeImages = True
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
#Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_SandboxNew/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output_SandboxNew5/"
Visu.transparency = True

Visu.showGlyphs = False
Visu.glyphMeshType = "Triangle"
Visu.glyphScale = 0.1/(BCStokes.refValue/(Char.length/Char.time))
glyphSpacing = (Grid.ymax-Grid.ymin)/8 #50 * km
Visu.glyphSamplingRateX = round(Grid.nxC/((Grid.xmax-Grid.xmin)/glyphSpacing))
Visu.glyphSamplingRateY = round(Grid.nyC/((Grid.ymax-Grid.ymin)/glyphSpacing))

Visu.height = 0.75 * Visu.height
Visu.width = 1.0* Visu.width

#Visu.filter = "Linear"
Visu.filter = "Nearest"

Visu.shiftFacY = -0.0
Visu.shiftFacZ = 0.1

RefVisc = PhaseRef.getRefVisc(0.0,Char.temperature,abs(BCStokes.backStrainRate))
print("dx = " + str((Grid.xmax-Grid.xmin)/Grid.nxC) + ", dy = " + str((Grid.ymax-Grid.ymin)/Grid.nyC))

RefP = PhaseRef.rho0*abs(Physics.gy)*(-Grid.ymin)/2.0
CharExtra = Input.CharExtra(Char)
Visu.colorMap.Stress.scale  = 20.0e6/CharExtra.stress
Visu.colorMap.Stress.center = 0*200.0e6/CharExtra.stress
Visu.colorMap.Stress.max    = 1.0
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 4.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.5
Visu.colorMap.Temperature.scale  = 1.0
Visu.colorMap.Temperature.center = 273.0/Char.temperature
Visu.colorMap.Temperature.max    = 1.0


Visu.colorMap.Pressure.scale  = 50e6/CharExtra.stress
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

Visu.colorMap.Vorticity.max = 0.0001/yr /  (1.0/Char.time) # in rad/yr



##          Write the input file
## =====================================
Input.writeInputFile(Setup)

