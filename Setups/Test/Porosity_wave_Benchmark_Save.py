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
Physics = Physics(True)
Visu = Visu()
Char = Char()
BCStokes = BCStokes()
BCThermal = BCThermal()
Geometry = {}



##       Modify Material properties
## =====================================
Phase0 = Material("Sediments")
Phase1 = Material()
Phase2 = Material()
Phase3 = Material()
Phase4 = Material()


Phase0.name = "Matrix"
Phase0.eta0 = 1e20


MatProps = {'0': Phase0.__dict__}



#BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -0.
Numerics.CFL_fac = 0.5
Numerics.nLineSearch = 3
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 10

Numerics.absoluteTolerance = 1e-4

#Numerics.dtMax = 20000000000.0

#Grid.nyC = [64]
#Grid.nxC = [32]

# Characteristic length and time for the porosity wave
Backphi = 0.001
Aphi = 0.01 # peak amplitude of the gaussian


RefPerm = 5e-18 ##Phase0.perm0# * Aphi*Aphi*Aphi  *  (1.0-Aphi)*(1.0-Aphi)
Phase0.perm0 = 5e-18/(Backphi * Backphi *Backphi  *  (1.0-Backphi)*(1.0-Backphi))
CompactionLength = sqrt(4/3*RefPerm/Physics.eta_f * (Phase0.eta0/Backphi))
DeltaRho = (1000-Phase0.rho0)
#CompactionVelocity = (DeltaRho * Physics.gy * CompactionLength*CompactionLength)/(Phase0.eta0/Backphi)
#PercolationVelocity = RefPerm*DeltaRho*Physics.gy
#C = (2*Aphi+1)*PercolationVelocity


Grid.xmin = -30*CompactionLength
Grid.xmax =  30*CompactionLength
Grid.ymin =  5*Grid.xmin
Grid.ymax =  5*Grid.xmax

RefinementFac = 1.0


Grid.nyC = round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nxC = round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)





Visu.showParticles = False
BCStokes.SetupType = "PureShear"
#BCStokes.SetupType = "SimpleShearPeriodic"
#BCThermal.SetupType = "SimpleShearPeriodic"

Particles.nPCX = 3
Particles.nPCY = 3

Visu.filter = "Nearest"

#Physics.gy = 0.
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)
Char.set_based_on_lithostatic_pressure(Phase0,BCThermal,Physics,Grid)

Char.length = CompactionLength
#Char.mass = Phase0.rho0*Char.length*Char.length*Char.length
CharStress =Phase0.rho0 *abs(Physics.gy)*Char.length
Char.time = Phase0.eta0/CharStress
Char.mass   = CharStress*Char.time*Char.time*Char.length






#Numerics.dtMax = 1/1000 * (1./RefinementFac   *  CompactionLength/C )/Char.time
#Numerics.dtMin = 1/1000 * (1./RefinementFac   *  CompactionLength/C )/Char.time

C2 = abs(RefPerm/(Physics.eta_f*Backphi) * (1-Backphi) * DeltaRho*Physics.gy)
#C3 = abs(RefPerm/(Physics.eta_f*Aphi) * (1-Aphi) * DeltaRho*Physics.gy)
C = (2*Aphi+1)
Numerics.dtMax = 1/10/RefinementFac*1/C
Numerics.dtMin = 1/10/RefinementFac*1/C


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
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC



Particles.noiseFactor = 0.95

Visu.height = 1/2 * Visu.height
Visu.width = 1 * Visu.width

Visu.type = "CompactionPressure"




##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

