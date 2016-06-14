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
physics = Physics()
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
Phase1.name = "Inclusion"
Phase1.eta0 = 1./100*Phase0.eta0
#Phase1.G = 100;

Phase1.rho0 = 1*Phase0.rho0
Phase2.eta0 = 1./100*Phase0.eta0
Phase0.n = 3.0
Phase1.n = 3.0

MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__, '2': Phase2.__dict__,'3': Phase3.__dict__,'4': Phase4.__dict__}



##            Define Numerics
## =====================================
Numerics.nTimeSteps = 1
BCStokes.backStrainRate = -1.0E-1
Numerics.CFL_fac = 0.5
Numerics.nLineSearch = 1
Numerics.maxNonLinearIter = 20

Numerics.absoluteTolerance = 1e-20

Grid.nxC = 256
Grid.nyC = 256

Visu.showParticles = False



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
Geometry["%05d_circle" % i] = vars(Geom_Circle(phase,0.,0.,0.2))

Visu.particleMeshRes = 6
Visu.particleMeshSize = 0.05

Particles.noiseFactor = 0.5



Visu.type = "Stress"


##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

