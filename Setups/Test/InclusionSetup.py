# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *
#from GeometryGraphical import *

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
Physics = Physics(True)
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

PhaseRef = Material()



PhaseRef.name = "Reference"

Phase0.name = "StickyAir"

Phase1.name = "Sediments"
Phase1.eta0 = 100.


MatProps = {'0': Phase0.__dict__,'1': Phase1.__dict__}



#BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = 1
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac = 1.0
Numerics.nLineSearch = 1
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 2

Numerics.absoluteTolerance = 5e-3

Grid.nyC = 16
Grid.nxC = 32

#Grid.xmin = -25.0e3
#Grid.xmax =  25.0e3
#Grid.ymax =  10.0e3
#Grid.ymin = 0

Visu.showParticles = False
#BCStokes.SetupType = "PureShear"
#BCStokes.SetupType = "Sandbox"
#BCThermal.SetupType = "Sandbox"

Particles.nPCX = 5
Particles.nPCY = 5

Visu.filter = "Nearest"

#Physics.gy = 0.
Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
#Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid)

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
Geometry["%05d_circle" % i] = (Geom_Circle(phase,0.0,0.0,0.25))




##for key in Geometry:
##    Geometry[key].plot()
##
##plt.axis([Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax])
##plt.show()




#make dict of geometry
for key in Geometry:
   Geometry[key] = vars(Geometry[key])









Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.0*(Grid.xmax-Grid.xmin)/Char.length/Grid.nxC



Particles.noiseFactor = 0.95

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "Viscosity"




##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

