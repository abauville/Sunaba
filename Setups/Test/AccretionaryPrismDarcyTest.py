# Input Test for Stokes FD
import sys
sys.path.insert(0, '../../src/UserInput')
import json
from InputDef import *

# Optional: uncomment the next line to activate the plotting methods to the Geometry objects, requires numpy and matplotlib
from GeometryGraphical import *

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
Phase0   = Material("StickyAir")
Phase1   = Material("StickyWater")
Phase2   = Material("Sediments")





MatProps = {'0': Phase0.__dict__, '1': Phase1.__dict__, '2': Phase2.__dict__}



#BCThermal.TT = 0.



##            Define Numerics
## =====================================
Numerics.nTimeSteps = -1
BCStokes.backStrainRate = -1.0
Numerics.CFL_fac = 0.1
Numerics.nLineSearch = 1
Numerics.maxCorrection  = 1.0
Numerics.maxNonLinearIter = 1

Numerics.absoluteTolerance = 1e-10

Numerics.dtMax = 20000000000.0

Grid.nyC = [16]
Grid.nxC = [16]

Grid.ymin =  0
Grid.ymax =  30.0E3
Grid.xmin =  0.
Grid.xmax =  30.0E3

Visu.showParticles = True
BCStokes.SetupType = "PureShear"
#BCStokes.SetupType = "SimpleShearPeriodic"
#BCThermal.SetupType = "SimpleShearPeriodic"

Particles.nPCX = 3
Particles.nPCY = 3

Visu.filter = "Nearest"

#Physics.gy = 0.
#Char.set_based_on_strainrate(Phase0,BCStokes,BCThermal,Grid)


##            Define Geometry
## =====================================
H = Grid.ymax-Grid.ymin
L = Grid.xmax-Grid.xmin
Hsed = 10.0E3
DepthWater = 10.0E3
TopWater = Hsed+DepthWater

air = 0
water = 1
sediments = 2

i = 0
Geometry["%05d_rect" % i] = Geom_Rect(water,0.0,Hsed,L,DepthWater)
i+=1
Geometry["%05d_line" % i] = Geom_Line(sediments,0.0,Hsed,"y","<",Grid.xmin,Grid.xmax)

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





Char.set_based_on_lithostatic_pressure(PhaseRef,BCThermal,Physics,Grid,3*Hsed)




print(Char.length)
print(PhaseRef.perm0)
print((PhaseRef.perm0/Char.length/Char.length) / (Physics.eta_f/(Char.mass/Char.length/Char.time/Char.time)))
print(Physics.eta_f)


Visu.particleMeshRes = 6
Visu.particleMeshSize = 0.4*(Grid.xmax-Grid.xmin)/Grid.nxC[0]



Particles.noiseFactor = 0.95

Visu.height = round((Grid.ymax-Grid.ymin)/(Grid.xmax-Grid.xmin) * Visu.height)
Visu.width = 1 * Visu.width

Visu.type = "CompactionPressure"




##          Write the input file
## =====================================
myJsonFile = dict(Description = Description, Grid = Grid.__dict__, Numerics = Numerics.__dict__, Particles = Particles.__dict__, Physics = Physics.__dict__, Visu = Visu.__dict__, MatProps = MatProps, Char = Char.__dict__, BCStokes = BCStokes.__dict__, BCThermal = BCThermal.__dict__, Geometry = Geometry);

outFile = open('input.json', 'w')
json.dump(myJsonFile, open('input.json', 'w') , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)

