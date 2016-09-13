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
PhaseRef = Material("Sediments","Reference")
Phase0   = Material("Sediments")
#Phase1   = Material("StickyWater")
Phase1   = Material("Sediments")

Phase0.eta0 = 1e19
Backphi = 0.01
RefPerm = 5e-20
Phase0.perm0 = RefPerm/(Backphi * Backphi *Backphi  /  (1.0-Backphi)*(1.0-Backphi))
Phase1.perm0 = Phase0.perm0


MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__}#, '2': Phase2.__dict__}


#BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac = 0.5
Numerics.nLineSearch = 1
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1

Numerics.absoluteTolerance = 1e-10

Numerics.dtMax = 20000000000.0

#Grid.nyC = [64]
#Grid.nxC = [128]

#Grid.ymin =  0
#Grid.ymax =  30.0E3
#Grid.xmin =  0.
#Grid.xmax =  30.0E3

Visu.showParticles = True
BCStokes.SetupType = "PureShear"
#BCStokes.SetupType = "SimpleShearPeriodic"
#BCThermal.SetupType = "SimpleShearPeriodic"

Particles.nPCX = 3
Particles.nPCY = 3

Visu.filter = "Nearest"

#Physics.gy = 0.
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)









#CompactionLength = sqrt(RefPerm * (Phase1.eta0/Backphi))
CompactionLength = sqrt(4/3* RefPerm/Physics.eta_f * (Phase1.eta0/Backphi))

#Char.length = CompactionLength
#CharStress =PhaseRef.rho0 *abs(Physics.gy)*Char.length
#Char.time = PhaseRef.eta0/CharStress
#Char.mass   = CharStress*Char.time*Char.time*Char.length


#Grid.xmin = -2*CompactionLength
#Grid.xmax =  2*CompactionLength
#Grid.ymin =  1*Grid.xmin
#Grid.ymax =  1*Grid.xmax

Grid.xmin = -100e3
Grid.xmax =  0.0
Grid.ymin =  0.0
Grid.ymax = 50e3;

Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid,0)

RefinementFac = 0.25


Grid.nyC = 256#round( RefinementFac*(Grid.ymax-Grid.ymin)/ CompactionLength)
Grid.nxC = 64#round( RefinementFac*(Grid.xmax-Grid.xmin)/ CompactionLength)

print("nxC = " + str(Grid.nxC))
print("nyC = " + str(Grid.nyC))




print("Char length: " + str(Char.length))
print("Phase Ref Perm0: " + str(PhaseRef.perm0))


Visu.particleMeshRes = 6
Visu.particleMeshSize = 1.5*(Grid.xmax-Grid.xmin)/Grid.nxC



Particles.noiseFactor = 0.95

Visu.height = round((Grid.ymax-Grid.ymin)/(Grid.xmax-Grid.xmin) * Visu.height)
Visu.width = 1 * Visu.width

Visu.type = "CompactionPressure"


##            Define Geometry
## =====================================
H = Grid.ymax-Grid.ymin
L = Grid.xmax-Grid.xmin
Hsed = H/2.0
DepthWater = H/2.0
TopWater = Hsed+DepthWater

air = 0
#water = 1
sediments = 1

i = 0
#Geometry["%05d_rect" % i] = Geom_Rect(sediments,0.0,Hsed/2.0,L,Hsed/2.0)
#i+=1
#Geometry["%05d_line" % i] = Geom_Line(sediments,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)
#Geometry["%05d_sine" % i] = Geom_Sine(sediments,Hsed,Hsed/20.0,0.,L/5.0,"y","<",Grid.xmin,Grid.xmax)

i+=1
#Geometry["%05d_line" % i] = Geom_Line(air,0.0,Hsed/1.1,"y","<",Grid.xmin,Grid.xmax)
#Geometry["%05d_sine" % i] = Geom_Sine(air,Hsed-Hsed/5.0,Hsed/20.0,0.,L/5.0,"y","<",Grid.xmin,Grid.xmax)



##i+=1
##phase = 1
##Geometry["%05d_rect" % i] = vars(Geom_Rect(phase,.5,.5,.2,.2))
##i+=1
##phase = 2
##Geometry["%05d_sine" % i] = vars(Geom_Sine(phase,-.2,0.2,0.,1.,"y","<",Grid.xmin,Grid.xmax))
##i+=1
##phase = 4
##Geometry["%05d_circle" % i] = vars(Geom_Circle(phase,-.5,-.5,0.2))
##



#plt.axis([Grid.xmin, Grid.xmax, Grid.ymin, Grid.ymax])

#for key in Geometry:
#    Geometry[key].plot()

#plt.show()




#make dict of geometry
for key in Geometry:
   Geometry[key] = vars(Geometry[key])


##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

