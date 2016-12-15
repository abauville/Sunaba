# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import InputDef as Input
#from GeometryGraphical import *

# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import * 

print("\n"*5)

Setup = Input.Setup(isDimensional=True)

## Description
## =====================================
Setup.Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"


##      Declare singleton objects
## =====================================

Grid = Setup.Grid
Numerics = Setup.Numerics
Particles = Setup.Particles
Physics = Setup.Physics
Visu = Setup.Visu
Char = Setup.Char
BCStokes = Setup.BC.Stokes
BCThermal = Setup.BC.Thermal
MatProps = Setup.MatProps
Geometry = Setup.Geometry



##       Modify Material properties
## =====================================
#Phase0 = Input.Material("Wet_Olivine")
Phase0 = Input.Material()
Phase1 = Input.Material()
Setup.MatProps = {"0":Phase0, "1":Phase1}

PhaseRef = Phase0
PhaseRef.isRef = True

PhaseRef.name = "Reference"

Phase0.name = "Matrix"

Phase1.name = "Inclusion"
Phase1.vDisl.B = 1.0/(1.0/1000.0)/2.0
Phase1.vDisl.E = 0.
Phase1.vDisl.V = 0.
#Phase0.vDisl.B = 1.0/2.0
#Phase0.vDisl.E = 0.
#Phase0.vDisl.V = 0.
#Phase0.vDisl.n    = 15.0
Phase1.vDisl.n    = 25.0

Phase0.vDiff.isActive = False
Phase0.vPei.isActive = False


##            Define Numerics
## =====================================
Numerics.nTimeSteps = 2
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac_Stokes = 0.5
Numerics.nLineSearch = 3 
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1000

Numerics.absoluteTolerance = 1e-6

Grid.nyC = 128
Grid.nxC = Grid.nyC

Visu.showParticles = False


Particles.nPCX = 4
Particles.nPCY = 4


Char.set_based_on_strainrate(PhaseRef,BCStokes,BCThermal,Grid)


##            Define Geometry
## =====================================

W = Grid.xmax-Grid.xmin
H = 0.6*(Grid.ymax-Grid.ymin)

DetHL = 0.25*H
DetHR = 0.15*H

InterH = 0.15*H
InterY = Grid.ymin+0.6*H-InterH/2


i = 0
phase = 1
#Geometry["%05d_line" % i] = (Geom_Line(phase,0.0,H,"y","<",Grid.xmin,Grid.xmax))
Geometry["%05d_circle" % i] = (Input.Geom_Circle(phase,0.0,0.0,0.33/2.0))


Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.0*(Grid.xmax-Grid.xmin)/Char.length/Grid.nxC



Particles.noiseFactor = 0.0

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "StrainRate"
#Visu.filter = "Linear"
Visu.filter = "Nearest"


if PhaseRef.vDisl.isActive:
    RefVisc = 1.0/(2.0*pow(PhaseRef.vDisl.B,1.0/PhaseRef.vDisl.n)*pow(abs(BCStokes.backStrainRate),(-1.0/PhaseRef.vDisl.n+1.0)))
elif PhaseRef.vDiff.isActive:
    RefVisc = PhaseRef.vDiff.B

CharExtra = Input.CharExtra(Char)
Visu.colorMap.Stress.scale  = 1.0
Visu.colorMap.Stress.center = 1.0
Visu.colorMap.Stress.max    = 1.75
Visu.colorMap.Viscosity.max = 0.5
Visu.colorMap.Viscosity.scale = RefVisc/CharExtra.visc
Visu.colorMap.Viscosity.max = 3.0
Visu.colorMap.StrainRate.scale = abs(BCStokes.backStrainRate/(1.0/Char.time))
Visu.colorMap.StrainRate.max = 1.0


###          Write the input file
### =====================================
Input.writeInputFile(Setup)


