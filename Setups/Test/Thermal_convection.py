# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *

# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *

print("\n"*5)

## Description
## =====================================
Description = "This is a test input file. Which defines to materials: a matrix and an inclusion 100 times stronger in a square box in pure shear"


##      Declare singleton objects
## =====================================
Grid = Grid()
Numerics = Numerics()
Particles = Particles()
Physics = Physics(False)
Visu = Visu()
Char = Char()
BCStokes = BCStokes()
BCThermal = BCThermal()
Geometry = {}



##       Modify Material properties
## =====================================
Phase0 = Material()
Phase1 = Material()
Phase2 = Material()
Phase3 = Material()
Phase4 = Material()


Phase0.name = "Matrix"

Phase0.alpha = 0.2
Phase0.beta  = 0.0
Phase0.k     = 8e-8

Phase0.eta0 = 1.

Phase0.perm0 = 1.



#Phase1.rho0 = 1.5*Phase0.rho0

MatProps = {'0': Phase0.__dict__}



BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -0.
Numerics.CFL_fac = 1.0
Numerics.nLineSearch = 1
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1

Numerics.absoluteTolerance = 1e-20

Numerics.dtMax = 50.0


Grid.nyC = [64]
Grid.nxC = [128]


Grid.xmin = -2.0
Grid.xmax =  2.0

Visu.showParticles = False
#BCStokes.SetupType = "PureShear"
BCStokes.SetupType = "SimpleShearPeriodic"
BCThermal.SetupType = "SimpleShearPeriodic"

Particles.nPCX = 3
Particles.nPCY = 3

Visu.filter = "Nearest"

#Physics.gy = 0.
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
Char.set_based_on_lithostatic_pressure(Phase0,BCThermal,Physics,Grid)

##            Define Geometry
## =====================================
##i = 0
##phase = 2
##Geometry["%05d_line" % i] = vars(Geom_Line(phase,0.2,0,"y",">",Grid.xmin,Grid.xmax))
##i+=1
##phase = 1
##Geometry["%05d_rect" % i] = vars(Geom_Rect(phase,.5,.5,.2,.2))
##i+=1
##phase = 3
##Geometry["%05d_sine" % i] = vars(Geom_Sine(phase,-.2,0.2,0.,1.,"y","<",Grid.xmin,Grid.xmax))
##i+=1
##phase = 4
##Geometry["%05d_circle" % i] = vars(Geom_Circle(phase,-.5,-.5,0.2))
##

i=0
phase = 1
#Geometry["%05d_circle" % i] = vars(Geom_Circle(phase,0.,0.,0.2))

Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC[0]



Particles.noiseFactor = 0.9

Visu.height = 1/2 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Temperature"




##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

