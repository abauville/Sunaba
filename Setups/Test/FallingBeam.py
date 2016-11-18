# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *
# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
#from GeometryGraphical import *

print("\n"*5)

##             Description
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


##          Material properties
## =====================================
Phase0 = Material("Default")
Phase1   = Material("Default")
PhaseRef = Phase0


Phase0.n = 5.0
Phase1.eta0 = 100
Phase1.rho0 = 3


Phase0.isRef = True

MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__}




##                 BC
## =====================================
BCStokes.SetupType = "PureShear"




##              Particles
## =====================================
Particles.noiseFactor = 0.95
Particles.nPCX = 3
Particles.nPCY = 3





##              Grid
## =====================================
#Grid.xmin = -32*CompactionLength
#Grid.xmax =  32*CompactionLength
#Grid.ymin =  -1.0*(Grid.xmax-Grid.xmin)
#Grid.ymax =  0.0*(Grid.xmax-Grid.xmin)

#RefinementFac = 2.0
#Grid.nyC = round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
#Grid.nxC = round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

Grid.xmin = -10.0
Grid.xmax =  10.0
Grid.ymin = -10.0
Grid.ymax =  10.0;
Grid.nxC = 256
Grid.nyC = 256

Grid.fixedBox = False

##              Non Dim
## =====================================
#Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid)






##              Geometry
## =====================================

H = 10.0
L = 2.0
A = 0.5

i=0
Geometry["%05d_sine" % i] = Geom_Sine(1,0+L/2.0,A,-pi/2.0,H/5.0,"x","<",-H/2.0,H/2.0)
i+=1
Geometry["%05d_sine" % i] = Geom_Sine(0,0-L/2.0,A,pi/2,H/5.0,"x","<",-H/2.0,H/2.0)


#plt.axis([Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax])

#for key in Geometry:
#    Geometry[key].plot()

#plt.show()

#make dict of geometry
for key in Geometry:
   Geometry[key] = vars(Geometry[key])




##            Visualization
## =====================================
Visu.showParticles = False
Visu.filter = "Nearest"
Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC

Visu.height = 1 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "StrainRate"
Visu.writeImages = False
#Visu.outputFolder = "/Users/abauville/JAMSTEC/StokesFD_OutputTest/"
Visu.outputFolder = "/Users/abauville/GoogleDrive/Output/"
Visu.transparency = True



##              Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = 0.0e-15
Numerics.CFL_fac = 0.05
Numerics.nLineSearch = 10
Numerics.maxCorrection  = 1.0
Numerics.minNonLinearIter = 1
Numerics.maxNonLinearIter = 15

Numerics.absoluteTolerance = 1e-5

Numerics.etaMin = 1e-5



##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);
outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)


